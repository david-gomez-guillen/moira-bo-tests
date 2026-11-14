# Individual-level microsimulation of a progressive disease: a birth cohort
# followed in annual cycles from age 20 to 85. Each person is in one of five
# states, and every year draws one uniform number per person that decides both
# other-cause death and the disease transition of that year.

START.AGE <- 20
END.AGE <- 85
AGES <- START.AGE:(END.AGE - 1)

# People simulated per evaluation. Outputs are Monte Carlo estimates, so two
# runs with the same parameters differ unless the random state is the same.
N.POPULATION <- 100000

HEALTHY <- 1L
LESION <- 2L
PRECLINICAL <- 3L
CLINICAL <- 4L

# Annual probability of death from other causes, by five-year age band.
LIFE.TABLE <- data.frame(
  age.from = seq(20, 80, by = 5),
  q = c(0.0004, 0.0005, 0.0006, 0.0009, 0.0014, 0.0022, 0.0034,
        0.0052, 0.0080, 0.0125, 0.0200, 0.0340, 0.0600)
)

# The strata: every natural-history probability takes one value per age group.
MODEL.STRATA <- c('20-39', '40-54', '55-69', '70-84')
STRATUM.START <- c(20, 40, 55, 70)
STRATUM.OF.AGE <- findInterval(AGES, STRATUM.START)

# What the model reports, by stratum and averaged over its years of age:
# clinical incidence per 100,000 person-years, and the share of the living with
# a lesion and with an undiagnosed cancer. Also by stratum, the persistence of
# lesions: of the people with a lesion at the follow-up age of the stratum, the
# share of those still alive five years later who have one then.
FOLLOW.UP.AGES <- c(25, 40, 55, 70)
FOLLOW.UP.YEARS <- 5

other.mortality <- function(ages) {
  LIFE.TABLE$q[findInterval(ages, LIFE.TABLE$age.from)]
}

# The value of a stratified parameter at every age of the simulation.
by.age <- function(value) {
  as.numeric(unlist(value)[MODEL.STRATA])[STRATUM.OF.AGE]
}

# The age-specific transition probabilities of every cycle.
transition.schedule <- function(pars) {
  list(
    q = other.mortality(AGES),
    onset = by.age(pars$p.lesion.onset),
    regress = by.age(pars$p.lesion.regress),
    progress = by.age(pars$p.lesion.progress),
    symptoms = by.age(pars$p.symptoms)
  )
}

# Runs the cohort and returns, per cycle, the people alive, with a lesion,
# with a preclinical cancer and diagnosed at its start, and the clinical
# diagnoses during it; and, per follow-up, the people with a lesion at its
# start and, among them, those alive and with a lesion at its end.
#
# The draw u of a person decides death first (u < q); otherwise v, u rescaled
# to [0, 1), decides the disease: a lesion regresses when v < p.lesion.regress
# and progresses when it falls in the next p.lesion.progress, and healthy and
# preclinical people move on when v < their own probability.
simulate.cohort <- function(pars, n = N.POPULATION) {
  s <- transition.schedule(pars)
  state <- rep.int(HEALTHY, n)

  followed <- NULL
  n.cycles <- length(AGES)
  alive <- lesions <- preclinical <- cases <- diagnosed <- numeric(n.cycles)
  follow.alive <- follow.lesion <- numeric(length(FOLLOW.UP.AGES))
  for (k in seq_len(n.cycles)) {
    f <- match(AGES[k] - FOLLOW.UP.YEARS, FOLLOW.UP.AGES)
    if (!is.na(f)) {
      follow.alive[f] <- sum(followed)
      follow.lesion[f] <- sum(followed & state == LESION)
    }
    # Who is followed: those with a lesion from the start of a follow-up to its
    # end, and no one in between follow-ups.
    followed <- if (AGES[k] %in% FOLLOW.UP.AGES) state == LESION
                else if (is.na(f)) followed

    counts <- tabulate(state, nbins = 4L)
    alive[k] <- length(state)
    lesions[k] <- counts[LESION]
    preclinical[k] <- counts[PRECLINICAL]
    diagnosed[k] <- counts[CLINICAL]

    u <- runif(length(state))
    v <- (u - s$q[k]) / (1 - s$q[k])
    down <- c(0, s$regress[k], 0, 0)[state]
    up <- c(s$onset[k], s$progress[k], s$symptoms[k], 0)[state]
    next.state <- state + (v >= down & v < down + up) - (v >= 0 & v < down)

    # Only the living are carried into the next cycle.
    survives <- u >= s$q[k]
    cases[k] <- sum(survives & state == PRECLINICAL & next.state == CLINICAL)
    state <- next.state[survives]
    if (!is.null(followed)) followed <- followed[survives]
  }
  list(alive = alive, lesions = lesions, preclinical = preclinical,
       cases = cases, diagnosed = diagnosed,
       follow.alive = follow.alive, follow.lesion = follow.lesion)
}

# The reported outputs of a cohort trace, in the layout of the targets.
cohort.outputs <- function(trace) {
  over <- function(count, scale = 1) {
    sapply(seq_along(MODEL.STRATA), function(i) {
      k <- STRATUM.OF.AGE == i
      scale * sum(count[k]) / sum(trace$alive[k])
    })
  }
  list(
    `Clinical incidence` = setNames(over(trace$cases, 1e5), MODEL.STRATA),
    `Lesion prevalence` = setNames(over(trace$lesions), MODEL.STRATA),
    `Preclinical cancer prevalence` = setNames(over(trace$preclinical), MODEL.STRATA),
    `Lesion persistence` = setNames(trace$follow.lesion / trace$follow.alive, MODEL.STRATA)
  )
}

# Costs and QALYs per person, discounted to age 20: a one-off cost at diagnosis,
# and a lower utility for every year lived after it.
cohort.summary <- function(trace, pars, n = N.POPULATION) {
  discount <- (1 + pars$discount)^-(AGES - START.AGE)
  undiagnosed <- trace$alive - trace$diagnosed
  qalys <- sum(discount * (undiagnosed + pars$utility.after.diagnosis * trace$diagnosed)) / n
  cost <- sum(discount * trace$cases * pars$cost.diagnosis) / n
  c(C = cost, E = qalys)
}

simulate <- function(strategies, pars) {
  outputs <- list()
  summary <- data.frame(strategy = character(0), C = numeric(0), E = numeric(0))
  for (strategy in strategies) {
    trace <- simulate.cohort(pars)
    outputs[[strategy]] <- cohort.outputs(trace)
    ce <- cohort.summary(trace, pars)
    summary <- rbind(summary, data.frame(strategy = strategy, C = ce[['C']], E = ce[['E']]))
  }
  list(summary = summary, outputs = outputs)
}
