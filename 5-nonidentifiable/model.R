# Markov cohort model of a cancer with two pathways, in annual cycles from age
# 30 to 89. The slow pathway runs through a low-grade and a high-grade lesion;
# the fast one goes straight from healthy to preclinical cancer. A preclinical
# cancer becomes clinical when it causes symptoms. Once diagnosed, a cancer has
# the same prognosis whichever pathway it came from: it is either cured or
# fatal.
#
# Strategies: no intervention, which the calibration is run on, or an
# examination at fixed ages, which removes the lesions it finds and treats the
# preclinical cancers it finds. A fast-pathway cancer has no lesion before it to
# remove and causes symptoms sooner, so screening does less against the fast
# pathway than against the slow one.

MODEL.STRATA <- c('30-44', '45-59', '60-74', '75-89')
START.AGE <- 30
END.AGE <- 89
AGES <- START.AGE:END.AGE
AGE.STRATUM <- MODEL.STRATA[findInterval(AGES, c(30, 45, 60, 75))]

STATES <- c('healthy', 'lgl', 'hgl', 'pre.slow', 'pre.fast',
            'clinical', 'cured', 'dead.cancer', 'dead.other')
ALIVE <- c('healthy', 'lgl', 'hgl', 'pre.slow', 'pre.fast', 'clinical', 'cured')
UNDIAGNOSED <- c('healthy', 'lgl', 'hgl', 'pre.slow', 'pre.fast')
LESION <- c('lgl', 'hgl')
PRECLINICAL <- c('pre.slow', 'pre.fast')

# The ages at which each strategy examines the cohort.
SCREENING.AGES <- list(
  no_intervention = integer(0),
  screening_5y = seq(45, 70, by = 5)
)

# The value of a parameter in a stratum, whether it arrives as one number or as
# a list (or vector) named by stratum.
stratum.value <- function(value, stratum) {
  if (length(value) == 1 && is.null(names(value))) return(as.numeric(unlist(value)))
  as.numeric(unlist(value)[[stratum]])
}

# Annual transitions among those who survive other causes, for one stratum.
disease.matrix <- function(pars, stratum) {
  v <- function(name) stratum.value(pars[[name]], stratum)
  P <- matrix(0, length(STATES), length(STATES), dimnames = list(STATES, STATES))

  P['healthy', 'lgl'] <- v('p.lgl.onset')
  P['healthy', 'pre.fast'] <- v('p.fast.onset')
  P['lgl', 'hgl'] <- v('p.lgl.progress')
  P['lgl', 'healthy'] <- v('p.lgl.regress')
  P['hgl', 'pre.slow'] <- v('p.hgl.progress')
  P['hgl', 'lgl'] <- v('p.hgl.regress')
  P['pre.slow', 'clinical'] <- v('p.symptoms.slow')
  P['pre.fast', 'clinical'] <- v('p.symptoms.fast')
  P['clinical', 'dead.cancer'] <- v('p.cancer.death')
  P['clinical', 'cured'] <- v('p.cure')

  if (any(!is.finite(P)) || any(P < 0)) stop('transition probabilities must be in [0, 1]')
  stay <- 1 - rowSums(P)
  if (any(stay < 0)) stop('the transitions out of a state add up to more than 1')
  diag(P) <- diag(P) + stay
  P['dead.cancer', ] <- 0
  P['dead.cancer', 'dead.cancer'] <- 1
  P['dead.other', ] <- 0
  P['dead.other', 'dead.other'] <- 1
  P
}

# Other-cause mortality by age, a Gompertz curve.
other.mortality <- function(pars, age) {
  pmin(1, pars$mortality.other.base * exp(pars$mortality.other.rate * (age - START.AGE)))
}

# Runs the cohort and returns its trace (a row per age, at the start of the
# year) and the yearly flows the outputs are counted from.
#
# At each age of `screening.ages` the year starts with an examination. Of those
# who attend, it finds a lesion with one sensitivity and a preclinical cancer
# with another. A lesion it finds is removed, which leaves the person healthy.
# A cancer it finds is diagnosed: treated while still without symptoms, it is
# cured at once with probability p.cure.screen.detected, and otherwise goes on
# as a diagnosed cancer.
simulate.cohort <- function(pars, screening.ages = integer(0)) {
  matrices <- lapply(MODEL.STRATA, function(s) disease.matrix(pars, s))
  names(matrices) <- MODEL.STRATA

  state <- setNames(numeric(length(STATES)), STATES)
  state['lgl'] <- pars$prevalence.lgl.start
  state['hgl'] <- pars$prevalence.hgl.start
  state['healthy'] <- 1 - state['lgl'] - state['hgl']

  n <- length(AGES)
  trace <- matrix(0, n, length(STATES), dimnames = list(AGES, STATES))
  diagnosed <- numeric(n)
  cancer.deaths <- numeric(n)
  screened <- removed <- screen.detected <- screen.cured <- numeric(n)
  for (i in seq_len(n)) {
    if (AGES[i] %in% screening.ages) {
      attending <- pars$screening.adherence
      found.lesions <- state[LESION] * attending * pars$screening.sens.lesion
      found.cancers <- state[PRECLINICAL] * attending * pars$screening.sens.cancer
      screened[i] <- attending * sum(state[UNDIAGNOSED])
      removed[i] <- sum(found.lesions)
      screen.detected[i] <- sum(found.cancers)
      screen.cured[i] <- screen.detected[i] * pars$p.cure.screen.detected
      state[LESION] <- state[LESION] - found.lesions
      state['healthy'] <- state[['healthy']] + removed[i]
      state[PRECLINICAL] <- state[PRECLINICAL] - found.cancers
      state['cured'] <- state[['cured']] + screen.cured[i]
      state['clinical'] <- state[['clinical']] + screen.detected[i] - screen.cured[i]
    }
    trace[i, ] <- state
    P <- matrices[[AGE.STRATUM[i]]]
    q <- other.mortality(pars, AGES[i])
    alive <- state[ALIVE]
    survivors <- state
    survivors[ALIVE] <- alive * (1 - q)
    diagnosed[i] <- survivors[['pre.slow']] * P['pre.slow', 'clinical'] +
      survivors[['pre.fast']] * P['pre.fast', 'clinical']
    cancer.deaths[i] <- survivors[['clinical']] * P['clinical', 'dead.cancer']
    state <- as.numeric(survivors %*% P)
    names(state) <- STATES
    state['dead.other'] <- state['dead.other'] + sum(alive) * q
  }
  list(trace = trace, diagnosed = diagnosed + screen.detected, cancer.deaths = cancer.deaths,
       screened = screened, removed = removed, screen.cured = screen.cured)
}

# The outputs the model is compared with registries and lesion surveys by, one
# value per age group:
# - clinical incidence, diagnoses per person-year among those never diagnosed,
#   at an examination as well as after symptoms;
# - cancer mortality, cancer deaths per person-year among the living;
# - lesion prevalence, the share of the never-diagnosed with a lesion of
#   either grade, averaged over the ages of the group.
cohort.outputs <- function(cohort) {
  trace <- cohort$trace
  undiagnosed <- rowSums(trace[, UNDIAGNOSED, drop = FALSE])
  alive <- rowSums(trace[, ALIVE, drop = FALSE])
  lesion <- rowSums(trace[, LESION, drop = FALSE])
  by.stratum <- function(x) tapply(x, factor(AGE.STRATUM, levels = MODEL.STRATA), sum)
  list(
    `Clinical incidence` = by.stratum(cohort$diagnosed) / by.stratum(undiagnosed),
    `Cancer mortality` = by.stratum(cohort$cancer.deaths) / by.stratum(alive),
    `Lesion prevalence` = by.stratum(lesion / undiagnosed) / as.numeric(table(AGE.STRATUM)[MODEL.STRATA])
  )
}

# Costs and QALYs of the cohort, discounted, for the summary the app expects:
# every year under treatment, and for a screening strategy the examinations,
# the lesions they remove and the cancers they find that are cured at once.
cohort.economics <- function(cohort, pars) {
  trace <- cohort$trace
  discount <- (1 + pars$discount)^-(AGES - START.AGE)
  in.treatment <- trace[, 'clinical']
  well <- rowSums(trace[, c('healthy', 'lgl', 'hgl', 'pre.slow', 'pre.fast', 'cured'), drop = FALSE])
  cost <- sum(discount * in.treatment * pars$cost.treatment) +
    sum(discount * (cohort$screened * pars$cost.screening +
                      cohort$removed * pars$cost.lesion.removal +
                      cohort$screen.cured * pars$cost.early.treatment))
  qalys <- sum(discount * (well + in.treatment * pars$utility.cancer))
  c(C = cost, E = qalys)
}

# The plots of the Base case tab for a cohort, each under the title of its tab:
# a numeric vector named by age, or a data frame of the ages and a line per
# column. Rates are among those never diagnosed, mortality among the living, as
# the outputs are.
#
# The new cancers of each pathway are the preclinical cancers it makes in a
# year, out of the trace and the transitions of each age: the fast pathway's
# share of them is what the constraints bound, and the data only pin down what
# the two pathways make together.
cohort.plots <- function(cohort, pars) {
  trace <- cohort$trace
  undiagnosed <- rowSums(trace[, UNDIAGNOSED, drop = FALSE])
  alive <- rowSums(trace[, ALIVE, drop = FALSE])
  survival <- 1 - other.mortality(pars, AGES)
  onset <- function(from, to) {
    p <- sapply(MODEL.STRATA, function(s) disease.matrix(pars, s)[from, to])
    trace[, from] * survival * p[AGE.STRATUM]
  }
  per.100k <- function(x) x / undiagnosed * 1e5
  list(
    `Clinical incidence per 100,000` = setNames(per.100k(cohort$diagnosed), AGES),
    `New cancers by pathway per 100,000` =
      data.frame(Age = AGES, `Slow pathway` = per.100k(onset('hgl', 'pre.slow')),
                 `Fast pathway` = per.100k(onset('healthy', 'pre.fast')), check.names = FALSE),
    `Lesion prevalence` =
      data.frame(Age = AGES, `Low-grade` = trace[, 'lgl'] / undiagnosed,
                 `High-grade` = trace[, 'hgl'] / undiagnosed, check.names = FALSE),
    `Cancer mortality per 100,000` = setNames(cohort$cancer.deaths / alive * 1e5, AGES)
  )
}

# The plots are drawn only when `make.plots` asks for them, so that every other
# run, the calibration's among them, returns what it did without them.
simulate <- function(strategies, pars, make.plots = FALSE) {
  outputs <- list()
  traces <- plots <- list()
  summary <- data.frame(strategy = character(0), C = numeric(0), E = numeric(0))
  for (strategy in strategies) {
    ages <- SCREENING.AGES[[strategy]]
    if (is.null(ages)) stop(sprintf("Unknown strategy '%s'", strategy))
    cohort <- simulate.cohort(pars, ages)
    outputs[[strategy]] <- cohort.outputs(cohort)
    traces[[strategy]] <- cohort$trace
    economics <- cohort.economics(cohort, pars)
    summary <- rbind(summary, data.frame(strategy = strategy, C = economics[['C']],
                                         E = economics[['E']]))
    if (make.plots) plots[[strategy]] <- cohort.plots(cohort, pars)
  }
  results <- list(summary = summary, outputs = outputs, cohort.info = traces)
  if (make.plots) results$plots <- plots
  results
}
