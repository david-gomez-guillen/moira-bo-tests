# A Markov cohort model of a cancer that grows out of a precancerous lesion.
#
#   Cycle:    one year
#   Horizon:  from age 20 to age 85
#   Strategies: natural history only, cancers are diagnosed when they give
#               symptoms, which the calibration is run on; or examinations at
#               fixed ages, which remove lesions and find preclinical cancers
#
# The cohort starts at age 20 free of lesions. Each year a healthy person can
# develop a lesion, a lesion can regress or become a preclinical cancer, and a
# preclinical cancer is localized first and advanced later, and is diagnosed
# when it gives symptoms at the stage it has reached. Everyone alive runs the
# same risk of dying of other causes.
#
# The file reads top to bottom: the states and age groups, the yearly
# transitions, simulate.cohort(), which runs them, and the outputs the
# calibration compares with data.


# ---- States and age groups ----------------------------------------------------

MODEL.AGE.START <- 20
MODEL.AGE.END <- 85
MODEL.AGES <- MODEL.AGE.START:(MODEL.AGE.END - 1)

MODEL.STATES <- c('healthy', 'lesion', 'pre.early', 'pre.late', 'diagnosed', 'dead.other')

# Everyone never diagnosed with cancer: incidence and lesion prevalence are
# measured over them.
AT.RISK <- c('healthy', 'lesion', 'pre.early', 'pre.late')

# The age strata of the age-dependent parameters, as the years of age each holds.
MODEL.STRATA.AGES <- list(`20-34` = 20:34, `35-49` = 35:49, `50-64` = 50:64, `65-84` = 65:84)
MODEL.STRATA <- names(MODEL.STRATA.AGES)

# The age groups and ages the outputs are reported at.
INCIDENCE.GROUPS <- list(`30-39` = 30:39, `40-49` = 40:49, `50-54` = 50:54, `55-59` = 55:59,
                         `60-64` = 60:64, `65-69` = 65:69, `70-74` = 70:74, `75-84` = 75:84)
PREVALENCE.GROUPS <- MODEL.STRATA.AGES
STAGE.GROUPS <- list(`30-54` = 30:54, `55-84` = 55:84)

# The ages at which each strategy examines the cohort.
SCREENING.AGES <- list(
  no_intervention = integer(0),
  screening_5y = seq(50, 70, by = 5)
)


# ---- One year of natural history ----------------------------------------------

# Annual probability of dying of something other than this cancer: a Gompertz
# curve, taken from life tables rather than calibrated.
other.cause.mortality <- function(age, base, rate) {
  pmin(1, base * exp(rate * (age - MODEL.AGE.START)))
}

# A parameter as one value per year of age. A stratified parameter holds a
# value per stratum, and one given as a single value applies at every age.
by.age <- function(value) {
  value <- unlist(value)
  if (length(value) == 1) return(rep(as.numeric(value), length(MODEL.AGES)))
  ages <- numeric(length(MODEL.AGES))
  for (stratum in MODEL.STRATA)
    ages[MODEL.AGES %in% MODEL.STRATA.AGES[[stratum]]] <- value[[stratum]]
  ages
}

# The yearly probabilities of every transition, one value per year of age.
transition.probabilities <- function(p) {
  list(
    onset = by.age(p$p.lesion.onset),
    regress = by.age(p$p.lesion.regress),
    progress = by.age(p$p.lesion.progress),
    stage = by.age(p$p.stage.progress),
    symptoms.early = by.age(p$p.symptoms.early),
    symptoms.late = by.age(p$p.symptoms.late),
    death.other = other.cause.mortality(MODEL.AGES, unlist(p$mortality.other.base),
                                        unlist(p$mortality.other.rate))
  )
}


# ---- The cohort ---------------------------------------------------------------

# Runs the cohort from age 20 to 85. Returns, for each year of age, the share
# of the cohort in each state at the start of the year and the cancers
# diagnosed during it, at the localized and at the advanced stage.
#
# At each age of `screening.ages` the year starts with an examination. Of those
# who attend, it finds a lesion with one sensitivity and a preclinical cancer
# with another. A lesion it finds is removed, which leaves the person healthy;
# a cancer it finds is diagnosed at the stage it has reached, and counts among
# the diagnoses of that year. Returns too, per year, the share of the cohort
# examined and the lesions removed.
simulate.cohort <- function(p, screening.ages = integer(0)) {
  tp <- transition.probabilities(p)
  if (any(!is.finite(unlist(tp))) || any(unlist(tp) < 0) || any(unlist(tp) > 1))
    stop('every transition probability must lie in [0, 1]')
  n <- length(MODEL.AGES)
  trace <- matrix(0, nrow = n, ncol = length(MODEL.STATES),
                  dimnames = list(MODEL.AGES, MODEL.STATES))
  diagnosed.early <- numeric(n)
  diagnosed.late <- numeric(n)
  screened <- removed <- screen.early <- screen.late <- numeric(n)

  healthy <- 1; lesion <- 0; pre.early <- 0; pre.late <- 0; diagnosed <- 0; dead <- 0
  for (i in seq_len(n)) {
    if (MODEL.AGES[i] %in% screening.ages) {
      attending <- unlist(p$screening.adherence)
      screened[i] <- attending * (healthy + lesion + pre.early + pre.late)
      removed[i] <- lesion * attending * unlist(p$screening.sens.lesion)
      screen.early[i] <- pre.early * attending * unlist(p$screening.sens.cancer)
      screen.late[i] <- pre.late * attending * unlist(p$screening.sens.cancer)
      healthy <- healthy + removed[i]
      lesion <- lesion - removed[i]
      pre.early <- pre.early - screen.early[i]
      pre.late <- pre.late - screen.late[i]
      diagnosed <- diagnosed + screen.early[i] + screen.late[i]
    }
    trace[i, ] <- c(healthy, lesion, pre.early, pre.late, diagnosed, dead)
    alive <- 1 - tp$death.other[i]

    # Those who survive the year make at most one transition out of their
    # state. The ways out of a lesion compete, and so do those out of a
    # localized cancer: when they add up to more than 1 they are scaled down.
    leave.lesion <- max(1, tp$regress[i] + tp$progress[i])
    leave.early <- max(1, tp$stage[i] + tp$symptoms.early[i])

    onset <- healthy * tp$onset[i]
    regress <- lesion * tp$regress[i] / leave.lesion
    progress <- lesion * tp$progress[i] / leave.lesion
    stage <- pre.early * tp$stage[i] / leave.early
    found.early <- pre.early * tp$symptoms.early[i] / leave.early
    found.late <- pre.late * tp$symptoms.late[i]

    diagnosed.early[i] <- alive * found.early
    diagnosed.late[i] <- alive * found.late

    living <- c(healthy + regress - onset,
                lesion + onset - regress - progress,
                pre.early + progress - stage - found.early,
                pre.late + stage - found.late,
                diagnosed + found.early + found.late)
    dead <- dead + (1 - alive) * sum(living)
    living <- alive * living
    healthy <- living[1]; lesion <- living[2]; pre.early <- living[3]
    pre.late <- living[4]; diagnosed <- living[5]
  }
  list(trace = trace,
       diagnosed.early = diagnosed.early + screen.early,
       diagnosed.late = diagnosed.late + screen.late,
       screened = screened, removed = removed)
}


# ---- Outputs ------------------------------------------------------------------

# What the calibration compares with data:
#   Cancer incidence:   diagnoses per 100,000 person-years among those never
#                       diagnosed, by age group: after symptoms, and under a
#                       screening strategy at an examination too
#   Lesion prevalence:  share of those never diagnosed who carry a lesion, by
#                       age group, as a study examining a sample of every age
#                       for lesions would find
#   Late-stage share:   share of the diagnoses of each age group made at the
#                       advanced stage
cohort.outputs <- function(cohort) {
  ages <- MODEL.AGES
  at.risk <- rowSums(cohort$trace[, AT.RISK, drop = FALSE])
  diagnosed <- cohort$diagnosed.early + cohort$diagnosed.late

  incidence <- sapply(INCIDENCE.GROUPS, function(group)
    1e5 * sum(diagnosed[ages %in% group]) / sum(at.risk[ages %in% group]))
  prevalence <- sapply(PREVALENCE.GROUPS, function(group)
    sum(cohort$trace[ages %in% group, 'lesion']) / sum(at.risk[ages %in% group]))
  late.share <- sapply(STAGE.GROUPS, function(group)
    sum(cohort$diagnosed.late[ages %in% group]) / sum(diagnosed[ages %in% group]))

  list(`Cancer incidence` = incidence,
       `Lesion prevalence` = prevalence,
       `Late-stage share` = late.share)
}

# Costs and QALYs of the cohort, discounted, for the cost-effectiveness tabs: a
# cost for every diagnosis and a lower utility for every year lived after one,
# whatever the stage, and the cost of the examinations and of the lesions they
# remove.
cohort.summary <- function(cohort, p) {
  discount <- (1 + unlist(p$discount))^-(MODEL.AGES - MODEL.AGE.START)
  diagnosed <- cohort$diagnosed.early + cohort$diagnosed.late
  undiagnosed <- rowSums(cohort$trace[, AT.RISK, drop = FALSE])
  qalys <- undiagnosed + unlist(p$utility.after.diagnosis) * cohort$trace[, 'diagnosed']
  screening <- cohort$screened * unlist(p$cost.screening) +
    cohort$removed * unlist(p$cost.lesion.removal)
  c(C = sum(discount * diagnosed * unlist(p$cost.diagnosis)) + sum(discount * screening),
    E = sum(discount * qalys))
}

# The plots of the Base case tab for a cohort, each under the title of its tab:
# a numeric vector named by age, or a data frame of the ages and a line per
# column. Rates and shares are among those never diagnosed, as the outputs are.
cohort.plots <- function(cohort) {
  at.risk <- rowSums(cohort$trace[, AT.RISK, drop = FALSE])
  by.stage <- function(early, late)
    data.frame(Age = MODEL.AGES, Localized = early, Advanced = late)
  list(
    `Cancer incidence per 100,000` =
      by.stage(cohort$diagnosed.early / at.risk * 1e5, cohort$diagnosed.late / at.risk * 1e5),
    `Lesion prevalence` = setNames(cohort$trace[, 'lesion'] / at.risk, MODEL.AGES),
    `Undiagnosed cancer prevalence per 100,000` =
      by.stage(cohort$trace[, 'pre.early'] / at.risk * 1e5,
               cohort$trace[, 'pre.late'] / at.risk * 1e5),
    `Cumulative diagnoses per 100,000` =
      by.stage(cumsum(cohort$diagnosed.early) * 1e5, cumsum(cohort$diagnosed.late) * 1e5)
  )
}

# Runs the cohort once for each strategy. `cohort` is the cohort of each. The
# plots are drawn only when `make.plots` asks for them, so that every other
# run, the calibration's among them, returns what it did without them.
simulate <- function(strategies, p, make.plots = FALSE) {
  summary <- data.frame(strategy = character(0), C = numeric(0), E = numeric(0))
  outputs <- cohorts <- plots <- list()
  for (strategy in strategies) {
    ages <- SCREENING.AGES[[strategy]]
    if (is.null(ages)) stop(sprintf("Unknown strategy '%s'", strategy))
    cohort <- simulate.cohort(p, ages)
    totals <- cohort.summary(cohort, p)
    summary <- rbind(summary, data.frame(strategy = strategy, C = unname(totals[['C']]),
                                         E = unname(totals[['E']])))
    outputs[[strategy]] <- cohort.outputs(cohort)
    cohorts[[strategy]] <- cohort
    if (make.plots) plots[[strategy]] <- cohort.plots(cohort)
  }
  results <- list(summary = summary, outputs = outputs, cohort = cohorts)
  if (make.plots) results$plots <- plots
  results
}
