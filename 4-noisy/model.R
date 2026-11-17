# Individual-level microsimulation of a progressive disease: a birth cohort
# followed in annual cycles from age 20 to 85. Each person is in one of five
# states, and every year draws one uniform number per person that decides both
# other-cause death and the disease transition of that year.
#
# Strategies: no intervention, which the calibration is run on, or an
# examination at fixed ages, which removes the lesions it finds and diagnoses
# the preclinical cancers it finds.

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

# The ages at which each strategy examines the cohort.
SCREENING.AGES <- list(
  no_intervention = integer(0),
  screening_5y = seq(50, 70, by = 5)
)

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
#
# At each age of `screening.ages` the year starts with an examination. Each
# person has a draw w for it, and with a lesion or a preclinical cancer is
# found when w is below the adherence times the sensitivity to what they have:
# a lesion found is removed, which leaves the person healthy, and a cancer
# found is diagnosed and counted among the cases of that year. The draws w come
# from a stream of their own, so the draws u are the same whether the cohort is
# examined or not, and so are the deaths from other causes.
simulate.cohort <- function(pars, n = N.POPULATION, screening.ages = integer(0)) {
  s <- transition.schedule(pars)
  state <- rep.int(HEALTHY, n)

  if (length(screening.ages)) {
    if (!exists('.Random.seed', envir = globalenv(), inherits = FALSE)) runif(1)
    main <- get('.Random.seed', envir = globalenv())
    set.seed(sample.int(.Machine$integer.max, 1))
    screening.draws <- matrix(runif(n * length(screening.ages)), n)
    assign('.Random.seed', main, envir = globalenv())
    person <- seq_len(n)
  }
  screened <- removed <- screen.found <- numeric(length(AGES))

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

    round <- match(AGES[k], screening.ages)
    if (!is.na(round)) {
      w <- screening.draws[person, round]
      lesion.found <- state == LESION & w < pars$screening.adherence * pars$screening.sens.lesion
      cancer.found <- state == PRECLINICAL &
        w < pars$screening.adherence * pars$screening.sens.cancer
      screened[k] <- sum(state != CLINICAL & w < pars$screening.adherence)
      removed[k] <- sum(lesion.found)
      screen.found[k] <- sum(cancer.found)
      state[lesion.found] <- HEALTHY
      state[cancer.found] <- CLINICAL
    }

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
    if (length(screening.ages)) person <- person[survives]
  }
  list(alive = alive, lesions = lesions, preclinical = preclinical,
       cases = cases + screen.found, diagnosed = diagnosed,
       follow.alive = follow.alive, follow.lesion = follow.lesion,
       screened = screened, removed = removed)
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
# after symptoms or at an examination, a lower utility for every year lived
# after it, and the cost of the examinations and of the lesions they remove.
cohort.summary <- function(trace, pars, n = N.POPULATION) {
  discount <- (1 + pars$discount)^-(AGES - START.AGE)
  undiagnosed <- trace$alive - trace$diagnosed
  qalys <- sum(discount * (undiagnosed + pars$utility.after.diagnosis * trace$diagnosed)) / n
  cost <- sum(discount * trace$cases * pars$cost.diagnosis) / n +
    sum(discount * (trace$screened * pars$cost.screening +
                      trace$removed * pars$cost.lesion.removal)) / n
  c(C = cost, E = qalys)
}

# The plots of the Base case tab for a cohort trace, each under the title of
# its tab: a numeric vector named by age, or a data frame of the ages and a line
# per column. Rates and shares are among the living, as the outputs are, and
# carry the Monte Carlo noise of the outputs. Drawing them takes no random
# numbers.
cohort.plots <- function(trace, n = N.POPULATION) {
  list(
    `Clinical incidence per 100,000` = setNames(trace$cases / trace$alive * 1e5, AGES),
    `Lesion prevalence` = setNames(trace$lesions / trace$alive, AGES),
    `Preclinical cancer prevalence per 100,000` =
      setNames(trace$preclinical / trace$alive * 1e5, AGES),
    `Cumulative diagnoses per 100,000` = setNames(cumsum(trace$cases) / n * 1e5, AGES)
  )
}

# Every strategy is simulated from the random state the first one starts from:
# with common random numbers, the difference between two strategies is that of
# the strategies rather than of the draws. The plots are drawn only when
# `make.plots` asks for them, so that every other run, the calibration's among
# them, returns what it did without them.
simulate <- function(strategies, pars, make.plots = FALSE) {
  outputs <- plots <- list()
  summary <- data.frame(strategy = character(0), C = numeric(0), E = numeric(0))
  if (length(strategies) > 1) {
    if (!exists('.Random.seed', envir = globalenv(), inherits = FALSE)) runif(1)
    start <- get('.Random.seed', envir = globalenv())
  }
  for (strategy in strategies) {
    ages <- SCREENING.AGES[[strategy]]
    if (is.null(ages)) stop(sprintf("Unknown strategy '%s'", strategy))
    if (strategy != strategies[1]) assign('.Random.seed', start, envir = globalenv())
    trace <- simulate.cohort(pars, screening.ages = ages)
    outputs[[strategy]] <- cohort.outputs(trace)
    ce <- cohort.summary(trace, pars)
    summary <- rbind(summary, data.frame(strategy = strategy, C = ce[['C']], E = ce[['E']]))
    if (make.plots) plots[[strategy]] <- cohort.plots(trace)
  }
  results <- list(summary = summary, outputs = outputs)
  if (make.plots) results$plots <- plots
  results
}
