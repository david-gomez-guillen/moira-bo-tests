# Markov cohort model of a cancer with two natural histories. A birth cohort is
# followed in annual cycles from age 20 to 85.
#
# - Slow pathway: healthy -> low-grade lesion -> high-grade lesion -> preclinical
#   cancer. Both lesions can regress.
# - Fast pathway: healthy -> short-lived lesion -> preclinical cancer. Only people
#   carrying a susceptibility can take it, and clearing any lesion leaves them
#   immune to it.
#
# Both pathways end in the same preclinical cancer, which becomes clinical when
# it causes symptoms.
#
# Strategies: no intervention, which the calibration is run on, or an
# examination at fixed ages that removes the lesions it finds and diagnoses
# the preclinical cancers it finds.

START.AGE <- 20
END.AGE <- 85
AGES <- START.AGE:(END.AGE - 1)

MODEL.STRATA <- c('20-39', '40-54', '55-69', '70-84')
STRATUM.START <- c(20, 40, 55, 70)
STRATUM.OF.AGE <- findInterval(AGES, STRATUM.START)

# Clinical incidence is reported in two age groups per stratum, precursor
# prevalence in one, from the middle of the stratum. Cancer is too rare before
# 30 for the registry to give a rate.
INCIDENCE.GROUPS <- list(`30-34` = 30:34, `35-39` = 35:39, `40-47` = 40:47,
                         `48-54` = 48:54, `55-62` = 55:62, `63-69` = 63:69,
                         `70-77` = 70:77, `78-84` = 78:84)
PREVALENCE.GROUPS <- list(`32-36` = 32:36, `45-49` = 45:49, `60-64` = 60:64,
                          `75-79` = 75:79)

STATES <- c('healthy.susceptible', 'healthy', 'low.grade', 'high.grade', 'fast.lesion',
            'preclinical', 'clinical')
PRECURSORS <- c('low.grade', 'high.grade', 'fast.lesion')
UNDIAGNOSED <- c('healthy.susceptible', 'healthy', PRECURSORS, 'preclinical')

# The ages at which each strategy examines the cohort. A fast-pathway lesion
# lasts a little over a year, so an examination finds few of them unless
# examinations are frequent; the lesions of the slow pathway last for years.
SCREENING.AGES <- list(
  no_intervention = integer(0),
  screening_10y = seq(50, 70, by = 10),
  screening_5y = seq(50, 75, by = 5)
)

# The value of a parameter at every age of the simulation, whether it comes as
# one value or as one per stratum.
by.age <- function(value) {
  value <- unlist(value)
  if (length(value) == 1) return(rep(as.numeric(value), length(AGES)))
  as.numeric(value[MODEL.STRATA])[STRATUM.OF.AGE]
}

check.probability <- function(name, values) {
  if (anyNA(values) || any(values < 0 | values > 1))
    stop(sprintf('%s must be a probability at every age', name))
}

# Runs the cohort and returns, per cycle, the share of the cohort in each state
# at its start, the clinical diagnoses during it and the preclinical cancers
# each pathway produced during it, with the examinations of the strategy at the
# ages in `screening.ages`.
#
# An examination comes at the start of the year. Of those who attend, it finds
# a lesion of either pathway with one sensitivity and a preclinical cancer with
# another. A lesion it finds is removed, which leaves the person healthy and,
# like clearing it, no longer at risk of the fast pathway. A cancer it finds is
# diagnosed then, and counts as a diagnosis of that year.
run.cohort <- function(pars, screening.ages = integer(0)) {
  onset.slow <- by.age(pars$p.lgl.onset)
  onset.fast <- by.age(pars$p.fast.onset)
  hgl.progress <- by.age(pars$p.hgl.progress)
  check.probability('p.lgl.onset', onset.slow)
  check.probability('p.fast.onset', onset.fast)
  check.probability('p.hgl.progress', hgl.progress)
  if (any(onset.slow + onset.fast > 1))
    stop('The two onsets of a susceptible person cannot add up to more than 1')
  if (pars$p.lgl.progress + pars$p.lgl.regress > 1 ||
      any(hgl.progress + pars$p.hgl.regress > 1) ||
      pars$p.fast.progress + pars$p.fast.clear > 1)
    stop('The exits of a lesion cannot add up to more than 1')

  mortality <- pmin(1, pars$mortality.other.base *
                      exp(pars$mortality.other.rate * (AGES - START.AGE)))

  n <- length(AGES)
  trace <- matrix(0, n, length(STATES), dimnames = list(AGES, STATES))
  diagnoses <- slow.cancers <- fast.cancers <- numeric(n)
  screened <- removed <- screen.detected <- numeric(n)

  x <- setNames(numeric(length(STATES)), STATES)
  x['healthy.susceptible'] <- pars$prop.susceptible
  x['healthy'] <- 1 - pars$prop.susceptible

  for (i in seq_len(n)) {
    if (AGES[i] %in% screening.ages) {
      attending <- pars$screening.adherence
      found.lesions <- x[PRECURSORS] * attending * pars$screening.sens.lesion
      screened[i] <- attending * sum(x[UNDIAGNOSED])
      removed[i] <- sum(found.lesions)
      screen.detected[i] <- x[['preclinical']] * attending * pars$screening.sens.cancer
      x[PRECURSORS] <- x[PRECURSORS] - found.lesions
      x['healthy'] <- x[['healthy']] + removed[i]
      x['preclinical'] <- x[['preclinical']] - screen.detected[i]
      x['clinical'] <- x[['clinical']] + screen.detected[i]
    }
    trace[i, ] <- x
    to.preclinical.slow <- x[['high.grade']] * hgl.progress[i]
    to.preclinical.fast <- x[['fast.lesion']] * pars$p.fast.progress
    diagnoses[i] <- x[['preclinical']] * pars$p.symptoms
    slow.cancers[i] <- to.preclinical.slow
    fast.cancers[i] <- to.preclinical.fast

    y <- x
    y['healthy.susceptible'] <- x[['healthy.susceptible']] * (1 - onset.slow[i] - onset.fast[i])
    y['healthy'] <- x[['healthy']] * (1 - onset.slow[i]) +
      x[['low.grade']] * pars$p.lgl.regress +
      x[['fast.lesion']] * pars$p.fast.clear
    y['low.grade'] <- x[['low.grade']] * (1 - pars$p.lgl.progress - pars$p.lgl.regress) +
      (x[['healthy.susceptible']] + x[['healthy']]) * onset.slow[i] +
      x[['high.grade']] * pars$p.hgl.regress
    y['high.grade'] <- x[['high.grade']] * (1 - hgl.progress[i] - pars$p.hgl.regress) +
      x[['low.grade']] * pars$p.lgl.progress
    y['fast.lesion'] <- x[['fast.lesion']] * (1 - pars$p.fast.progress - pars$p.fast.clear) +
      x[['healthy.susceptible']] * onset.fast[i]
    y['preclinical'] <- x[['preclinical']] * (1 - pars$p.symptoms) +
      to.preclinical.slow + to.preclinical.fast
    y['clinical'] <- x[['clinical']] * (1 - pars$p.cancer.death) + diagnoses[i]
    x <- y * (1 - mortality[i])
  }
  list(trace = trace, diagnoses = diagnoses, slow.cancers = slow.cancers,
       fast.cancers = fast.cancers, screened = screened, removed = removed,
       screen.detected = screen.detected)
}

# Clinical incidence per 100,000 person-years of the undiagnosed, and the share
# of the undiagnosed with a precursor lesion, averaged over the ages of a group.
# Incidence counts every diagnosis, after symptoms or at an examination, as a
# registry would.
cohort.outputs <- function(cohort) {
  at.risk <- rowSums(cohort$trace[, UNDIAGNOSED, drop = FALSE])
  with.lesion <- rowSums(cohort$trace[, PRECURSORS, drop = FALSE]) / at.risk
  all.diagnoses <- cohort$diagnoses + cohort$screen.detected
  rows <- function(ages) match(ages, AGES)
  list(
    `Clinical incidence` = sapply(INCIDENCE.GROUPS, function(ages) {
      j <- rows(ages)
      sum(all.diagnoses[j]) / sum(at.risk[j]) * 1e5
    }),
    `Precursor prevalence` = sapply(PREVALENCE.GROUPS, function(ages) mean(with.lesion[rows(ages)]))
  )
}

# Costs and QALYs of a cohort, discounted to age 20: a one-off cost for every
# diagnosis, a lower utility for every year lived after it, and the cost of the
# examinations and of the lesions they remove.
cohort.summary <- function(cohort, pars) {
  discount <- (1 + pars$discount)^-(AGES - START.AGE)
  alive.undiagnosed <- rowSums(cohort$trace[, UNDIAGNOSED, drop = FALSE])
  alive.diagnosed <- cohort$trace[, 'clinical']
  cost <- sum((cohort$diagnoses + cohort$screen.detected) * pars$cost.diagnosis * discount) +
    sum((cohort$screened * pars$cost.screening + cohort$removed * pars$cost.lesion.removal) *
          discount)
  qalys <- sum((alive.undiagnosed + alive.diagnosed * pars$utility.after.diagnosis) * discount)
  c(C = cost, E = qalys)
}

# The share of the preclinical cancers of the cohort that each pathway made.
pathway.share <- function(cohort) {
  cancers <- sum(cohort$slow.cancers + cohort$fast.cancers)
  c(slow = sum(cohort$slow.cancers), fast = sum(cohort$fast.cancers)) /
    (if (cancers > 0) cancers else 1)
}

# The plots of the Base case tab for a cohort, each under the title of its tab:
# a numeric vector named by age, or a data frame of the ages and a line per
# column. Rates are among the undiagnosed, as the outputs are.
cohort.plots <- function(cohort) {
  at.risk <- rowSums(cohort$trace[, UNDIAGNOSED, drop = FALSE])
  per.100k <- function(x) x / at.risk * 1e5
  share <- function(state) cohort$trace[, state] / at.risk
  incidence <- data.frame(Age = AGES, `After symptoms` = per.100k(cohort$diagnoses),
                          check.names = FALSE)
  if (any(cohort$screen.detected > 0))
    incidence$`At an examination` <- per.100k(cohort$screen.detected)
  list(
    `Clinical incidence per 100,000` = incidence,
    # The preclinical cancers each pathway makes: where the two basins of the
    # calibration differ.
    `New cancers by pathway per 100,000` =
      data.frame(Age = AGES, `Slow pathway` = per.100k(cohort$slow.cancers),
                 `Fast pathway` = per.100k(cohort$fast.cancers), check.names = FALSE),
    `Precursor prevalence` =
      data.frame(Age = AGES, `Low-grade` = share('low.grade'), `High-grade` = share('high.grade'),
                 `Fast-pathway lesion` = share('fast.lesion'), check.names = FALSE),
    `Cumulative diagnoses per 100,000` =
      setNames(cumsum(cohort$diagnoses + cohort$screen.detected) * 1e5, AGES)
  )
}

# The plots are drawn only when `make.plots` asks for them, so that every other
# run, the calibration's among them, returns what it did without them.
simulate <- function(strategies, pars, make.plots = FALSE) {
  summary <- data.frame(strategy = character(0), C = numeric(0), E = numeric(0))
  outputs <- share <- traces <- plots <- list()
  for (strategy in strategies) {
    ages <- SCREENING.AGES[[strategy]]
    if (is.null(ages)) stop(sprintf("Unknown strategy '%s'", strategy))
    cohort <- run.cohort(pars, ages)
    ce <- cohort.summary(cohort, pars)
    summary <- rbind(summary, data.frame(strategy = strategy, C = ce[['C']], E = ce[['E']]))
    outputs[[strategy]] <- cohort.outputs(cohort)
    share[[strategy]] <- pathway.share(cohort)
    traces[[strategy]] <- cohort$trace
    if (make.plots) plots[[strategy]] <- cohort.plots(cohort)
  }
  results <- list(summary = summary, outputs = outputs, pathway.share = share,
                  cohort.info = traces)
  if (make.plots) results$plots <- plots
  results
}
