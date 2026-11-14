# A cost-effectiveness model of screening for a hypothetical cancer.
#
# It is a Markov cohort model. Rather than following people one by one, it keeps
# track of which share of a group of people, the cohort, is in each health state,
# and moves those shares from state to state once a year.
#
#   Cycle:       one year
#   Horizon:     from age 20 to age 89
#   Strategies:  no screening, or screening rounds at fixed ages, done with a
#                test or with a procedure
#
# Cancer grows out of a lesion, along one of two pathways:
#
#   slow:  normal -> low-grade lesion -> high-grade lesion -> cancer
#   fast:  normal -> fast-pathway lesion -> cancer
#
# A cancer starts out preclinical, which means nobody knows it is there. It is
# localized first and advanced later, and it is diagnosed when it gives symptoms
# or when a screening round finds it. Screening pays off in two ways: it removes
# lesions before they turn into cancer, and it finds cancers while they are
# still localized and easier to cure.
#
# The file reads top to bottom: the states, the yearly transitions, what a
# screening round does, and simulate(), which puts them together. Bookkeeping
# that is not needed to follow the model, such as reading a parameter's value at
# an age or summarizing results per age group, is in model_helpers.R.

source('model_helpers.R')


# ---- States ---------------------------------------------------------------------

MODEL.AGE.START <- 20
MODEL.AGE.END <- 89

# Each pathway has preclinical states of its own because its cancers behave
# differently: those of the fast pathway advance sooner and take longer to give
# symptoms, so more of them are diagnosed late.
MODEL.STATES <- c(
  'normal',          # no lesion
  'lgl',             # low-grade lesion
  'hgl',             # high-grade lesion
  'fpl',             # fast-pathway lesion
  'pre.early.slow',  # preclinical localized cancer, slow pathway
  'pre.late.slow',   # preclinical advanced cancer, slow pathway
  'pre.early.fast',  # preclinical localized cancer, fast pathway
  'pre.late.fast',   # preclinical advanced cancer, fast pathway
  'clin.early',      # diagnosed localized cancer, under treatment
  'clin.late',       # diagnosed advanced cancer, under treatment
  'survivor',        # alive after treatment
  'dead.cancer',     # death from cancer
  'dead.other'       # death from other causes
)

# Groups of states the model refers to more than once.
LESION <- c('lgl', 'hgl', 'fpl')
PRECLINICAL <- c('pre.early.slow', 'pre.late.slow', 'pre.early.fast', 'pre.late.fast')
DEAD <- c('dead.cancer', 'dead.other')
ALIVE <- setdiff(MODEL.STATES, DEAD)

# Everyone who has never been diagnosed with cancer. Incidence and lesion
# prevalence are measured over them.
AT.RISK <- c('normal', LESION, PRECLINICAL)

# The ages at which each strategy offers a screening round, and what with.
SCREENING.SCHEDULES <- list(
  no_screening = list(modality = NULL, ages = integer(0)),
  test_biennial = list(modality = 'test', ages = seq(50, 74, by = 2)),
  procedure_10y = list(modality = 'procedure', ages = seq(50, 79, by = 10)),
  procedure_45 = list(modality = 'procedure', ages = seq(45, 79, by = 10))
)

# Seconds simulate() waits before it starts. The model runs in milliseconds, and
# the wait makes it behave like the slow models it stands in for. 0 disables it.
MODEL.SIMULATION.DELAY <- 3


# ---- One year of natural history ------------------------------------------------

# Annual probability of dying of something other than this cancer. It grows
# exponentially with age, which is known as a Gompertz curve.
other.cause.mortality <- function(age, base, rate) {
  min(1, base * exp(rate * (age - MODEL.AGE.START)))
}

# The transition matrix of one year at a given age.
# Row: the state a person is in now. Column: the state a year later.
transition.matrix <- function(age, p) {
  tp <- matrix(0, nrow = length(MODEL.STATES), ncol = length(MODEL.STATES),
               dimnames = list(MODEL.STATES, MODEL.STATES))

  # A lesion appears, on one pathway or the other.
  tp['normal', 'lgl'] <- p$p.lgl.onset
  tp['normal', 'fpl'] <- p$p.fpl.onset

  # A lesion gets worse, or regresses to what it came from.
  tp['lgl', 'hgl'] <- p$p.lgl.progress
  tp['lgl', 'normal'] <- p$p.lgl.regress
  tp['hgl', 'pre.early.slow'] <- p$p.hgl.progress
  tp['hgl', 'lgl'] <- p$p.hgl.regress
  tp['fpl', 'pre.early.fast'] <- p$p.fpl.progress
  tp['fpl', 'normal'] <- p$p.fpl.regress

  # A preclinical cancer gives symptoms and is diagnosed at the stage it is in,
  # or advances from localized to advanced without being noticed.
  tp['pre.early.slow', 'clin.early'] <- p$p.symptomatic.early.slow
  tp['pre.early.slow', 'pre.late.slow'] <- p$p.stage.progress.slow
  tp['pre.late.slow', 'clin.late'] <- p$p.symptomatic.late.slow
  tp['pre.early.fast', 'clin.early'] <- p$p.symptomatic.early.fast
  tp['pre.early.fast', 'pre.late.fast'] <- p$p.stage.progress.fast
  tp['pre.late.fast', 'clin.late'] <- p$p.symptomatic.late.fast

  # A diagnosed cancer is cured or kills, and a survivor can still die of it.
  tp['clin.early', 'survivor'] <- p$p.cure.early
  tp['clin.early', 'dead.cancer'] <- p$p.cancer.death.early
  tp['clin.late', 'survivor'] <- p$p.cure.late
  tp['clin.late', 'dead.cancer'] <- p$p.cancer.death.late
  tp['survivor', 'dead.cancer'] <- p$p.survivor.death

  # A probability cannot be negative, whatever value a parameter is given.
  tp <- pmax(tp, 0)

  # Everyone alive runs the same risk of dying of something else.
  tp[ALIVE, 'dead.other'] <- other.cause.mortality(age, p$mortality.other.base,
                                                   p$mortality.other.rate)

  # Each row has to add up to 1, since everyone ends up somewhere. A row whose
  # transitions add up to more is scaled down, so that a calibration trying
  # extreme values still gets a valid matrix. What the transitions leave is the
  # probability of staying, on the diagonal.
  leaving <- rowSums(tp)
  excess <- leaving > 1
  if (any(excess)) tp[excess, ] <- tp[excess, ] / leaving[excess]
  diag(tp) <- 1 - rowSums(tp)

  # The dead stay dead: these are the absorbing states.
  tp['dead.cancer', 'dead.cancer'] <- 1
  tp['dead.other', 'dead.other'] <- 1

  tp
}


# ---- A screening round ------------------------------------------------------------

# Sensitivity: the probability that an examination finds what a person in each
# state is carrying. An advanced cancer is somewhat easier to find than a
# localized one (`late.boost`), and those of the fast pathway are harder to find
# than those of the slow one (`rr.fast`).
sensitivity.by.state <- function(lgl, hgl, fpl, cancer, late.boost, rr.fast) {
  sensitivity <- per.state()
  sensitivity['lgl'] <- lgl
  sensitivity['hgl'] <- hgl
  sensitivity['fpl'] <- fpl
  sensitivity[c('pre.early.slow', 'pre.early.fast')] <- cancer
  sensitivity[c('pre.late.slow', 'pre.late.fast')] <- min(1, cancer * late.boost)

  fast <- c('pre.early.fast', 'pre.late.fast')
  sensitivity[fast] <- sensitivity[fast] * rr.fast
  sensitivity
}

# For a person in each state, the probability that a round finds their lesion or
# cancer (`detect`) and the probability that they undergo a procedure
# (`procedures`).
screening.detection <- function(modality, p) {
  detect <- per.state()
  procedures <- per.state()

  sens.procedure <- sensitivity.by.state(
    p$sens.procedure.lgl, p$sens.procedure.hgl, p$sens.procedure.fpl,
    p$sens.procedure.cancer, late.boost = 1.05, rr.fast = p$rr.detection.fast)

  if (identical(modality, 'test')) {
    # The test comes first, and only those who test positive get the procedure
    # that confirms and removes. Finding something takes both to succeed.
    positive <- sensitivity.by.state(
      p$sens.test.lgl, p$sens.test.hgl, p$sens.test.fpl,
      p$sens.test.cancer, late.boost = 1.15, rr.fast = p$rr.detection.fast)

    # Specificity is the probability that someone with nothing tests negative,
    # so the rest of them are false positives and get a procedure for nothing.
    positive[c('normal', 'survivor')] <- 1 - p$spec.test

    detect <- p$adherence.test * positive * sens.procedure
    procedures <- p$adherence.test * positive
  } else if (identical(modality, 'procedure')) {
    # Everyone who attends gets the procedure straight away.
    detect <- p$adherence.procedure * sens.procedure
    procedures[setdiff(ALIVE, c('clin.early', 'clin.late'))] <- p$adherence.procedure
  }

  list(detect = pmin(detect, 1), procedures = pmin(procedures, 1))
}

# What a screening round does to the cohort: the lesions it finds are removed,
# which sends those people back to normal, and the preclinical cancers it finds
# are diagnosed at the stage they are in.
apply.screening <- function(cohort, modality, p) {
  round <- screening.detection(modality, p)
  detect <- round$detect

  removed <- cohort[LESION] * detect[LESION]
  found.early <- cohort[['pre.early.slow']] * detect[['pre.early.slow']] +
                 cohort[['pre.early.fast']] * detect[['pre.early.fast']]
  found.late <- cohort[['pre.late.slow']] * detect[['pre.late.slow']] +
                cohort[['pre.late.fast']] * detect[['pre.late.fast']]

  cohort[LESION] <- cohort[LESION] - removed
  cohort['normal'] <- cohort['normal'] + sum(removed)
  cohort[PRECLINICAL] <- cohort[PRECLINICAL] * (1 - detect[PRECLINICAL])
  cohort['clin.early'] <- cohort['clin.early'] + found.early
  cohort['clin.late'] <- cohort['clin.late'] + found.late

  # Whoever is alive and not being treated for cancer is invited to the round.
  n.tests <- if (identical(modality, 'test')) {
    p$adherence.test * sum(cohort[c(AT.RISK, 'survivor')])
  } else 0

  list(cohort = cohort,
       n.tests = n.tests,
       n.procedures = sum(cohort * round$procedures),
       n.removals = sum(removed),
       found.early = found.early,
       found.late = found.late)
}


# ---- The simulation -----------------------------------------------------------------

simulate <- function(strategies, pars, delay = MODEL.SIMULATION.DELAY) {
  if (delay > 0) Sys.sleep(delay)

  ages <- MODEL.AGE.START:MODEL.AGE.END

  # One discount rate for the whole horizon, the one of the starting age.
  discount <- par.at.age(pars$discount, MODEL.AGE.START)

  summary <- data.frame()
  outputs <- list()
  cohort.info <- list()

  for (strategy in strategies) {
    schedule <- SCREENING.SCHEDULES[[strategy]]
    if (is.null(schedule)) stop(sprintf("Unknown strategy '%s'", strategy))

    # The share of the cohort in each state. Everyone starts with no lesion,
    # and the shares always add up to 1.
    cohort <- per.state()
    cohort['normal'] <- 1

    # One value per year, filled in as the cohort ages. The trace keeps the
    # cohort of every year, a row per age.
    n.years <- length(ages)
    costs <- numeric(n.years)
    qalys <- numeric(n.years)
    at.risk <- numeric(n.years)
    with.lesion <- numeric(n.years)
    new.early <- numeric(n.years)
    new.late <- numeric(n.years)
    trace <- matrix(NA_real_, nrow = n.years, ncol = length(MODEL.STATES),
                    dimnames = list(as.character(ages), MODEL.STATES))

    for (i in seq_along(ages)) {
      age <- ages[i]
      p <- parameters.at.age(pars, age)

      # ---- 1. Screening round, in the years the strategy has one ------------
      # It happens at the start of the year. It costs the tests, the procedures
      # and the removals, and a procedure can cause a complication, which costs
      # money and takes away some quality of life.
      screening.cost <- 0
      screening.disutility <- 0
      found.early <- 0
      found.late <- 0

      if (age %in% schedule$ages) {
        round <- apply.screening(cohort, schedule$modality, p)
        cohort <- round$cohort
        found.early <- round$found.early
        found.late <- round$found.late

        complications <- round$n.procedures * p$p.procedure.complication
        screening.cost <- round$n.tests * p$cost.test +
          round$n.procedures * p$cost.procedure +
          round$n.removals * p$cost.removal +
          complications * p$cost.complication
        screening.disutility <- complications * p$disutility.complication
      }

      trace[i, ] <- cohort

      # ---- 2. Costs and health of this year ---------------------------------
      # What a person in each state costs in a year, and the quality of life of
      # that year: 1 is a year in full health and 0 is being dead. A lesion or
      # a preclinical cancer is not felt, so it counts as full health.
      state.costs <- per.state()
      state.costs['clin.early'] <- p$cost.treatment.early
      state.costs['clin.late'] <- p$cost.treatment.late
      state.costs['survivor'] <- p$cost.followup

      state.utilities <- per.state()
      state.utilities[AT.RISK] <- 1
      state.utilities['clin.early'] <- p$utility.cancer.early
      state.utilities['clin.late'] <- p$utility.cancer.late
      state.utilities['survivor'] <- p$utility.survivor

      # Each state contributes in proportion to the share of the cohort in it.
      # Money and health count for less the further in the future they are,
      # which is what discounting expresses: at a rate of 3%, what happens a
      # year from now is worth 1 / 1.03 of the same today.
      discount.factor <- 1 / (1 + discount)^(age - MODEL.AGE.START)
      costs[i] <- (sum(state.costs * cohort) + screening.cost) * discount.factor
      qalys[i] <- (sum(state.utilities * cohort) - screening.disutility) * discount.factor

      # ---- 3. What is measured this year ------------------------------------
      # New diagnoses, as a cancer registry counts them: those that surface
      # through symptoms during the year plus those the screening round found.
      tp <- transition.matrix(age, p)
      new.early[i] <- cohort[['pre.early.slow']] * tp['pre.early.slow', 'clin.early'] +
                      cohort[['pre.early.fast']] * tp['pre.early.fast', 'clin.early'] +
                      found.early
      new.late[i] <- cohort[['pre.late.slow']] * tp['pre.late.slow', 'clin.late'] +
                     cohort[['pre.late.fast']] * tp['pre.late.fast', 'clin.late'] +
                     found.late
      at.risk[i] <- sum(cohort[AT.RISK])
      with.lesion[i] <- sum(cohort[LESION])

      # ---- 4. Move the cohort one year forward ------------------------------
      # Multiplying the shares by the matrix sends each share to where its row
      # says, which gives the shares at the start of next year.
      cohort <- drop(cohort %*% tp)
    }

    # ---- Results of the strategy ----------------------------------------------
    # C is the discounted cost and E the discounted quality-adjusted life years
    # (QALYs), both per person and added up over the whole horizon.
    summary <- rbind(summary, data.frame(strategy = strategy,
                                         C = sum(costs),
                                         E = sum(qalys)))

    # The outputs a calibration compares with observed data, per age group.
    outputs[[strategy]] <- calibration.outputs(ages, new.early, new.late,
                                               at.risk, with.lesion)
    cohort.info[[strategy]] <- trace
  }

  list(summary = summary,
       outputs = outputs,
       incidence = lapply(outputs, function(o) o$`Cancer incidence`),
       cohort.info = cohort.info)
}
