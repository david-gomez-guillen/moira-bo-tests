# A Markov cohort model of a cancer that grows out of a precancerous lesion.
#
#   Cycle:    one year
#   Horizon:  from age 20 to age 85
#   Strategy: natural history only, cancers are diagnosed when they give symptoms
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
simulate.cohort <- function(p) {
  tp <- transition.probabilities(p)
  if (any(!is.finite(unlist(tp))) || any(unlist(tp) < 0) || any(unlist(tp) > 1))
    stop('every transition probability must lie in [0, 1]')
  n <- length(MODEL.AGES)
  trace <- matrix(0, nrow = n, ncol = length(MODEL.STATES),
                  dimnames = list(MODEL.AGES, MODEL.STATES))
  diagnosed.early <- numeric(n)
  diagnosed.late <- numeric(n)

  healthy <- 1; lesion <- 0; pre.early <- 0; pre.late <- 0; diagnosed <- 0; dead <- 0
  for (i in seq_len(n)) {
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
  list(trace = trace, diagnosed.early = diagnosed.early, diagnosed.late = diagnosed.late)
}


# ---- Outputs ------------------------------------------------------------------

# What the calibration compares with data:
#   Cancer incidence:   symptomatic diagnoses per 100,000 person-years among
#                       those never diagnosed, by age group
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
# cost for every diagnosis and a lower utility for every year lived after one.
cohort.summary <- function(cohort, p) {
  discount <- (1 + unlist(p$discount))^-(MODEL.AGES - MODEL.AGE.START)
  diagnosed <- cohort$diagnosed.early + cohort$diagnosed.late
  undiagnosed <- rowSums(cohort$trace[, AT.RISK, drop = FALSE])
  qalys <- undiagnosed + unlist(p$utility.after.diagnosis) * cohort$trace[, 'diagnosed']
  c(C = sum(discount * diagnosed * unlist(p$cost.diagnosis)), E = sum(discount * qalys))
}

# Runs the cohort once for each strategy. There is only natural history to
# simulate, so every strategy gives the same results.
simulate <- function(strategies, p) {
  cohort <- simulate.cohort(p)
  outputs <- cohort.outputs(cohort)
  totals <- cohort.summary(cohort, p)
  summary <- data.frame(strategy = strategies, C = unname(totals[['C']]),
                        E = unname(totals[['E']]))
  list(summary = summary,
       outputs = setNames(rep(list(outputs), length(strategies)), strategies),
       cohort = cohort)
}
