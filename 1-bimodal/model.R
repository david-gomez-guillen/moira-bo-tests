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
# each pathway produced during it.
run.cohort <- function(pars) {
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

  x <- setNames(numeric(length(STATES)), STATES)
  x['healthy.susceptible'] <- pars$prop.susceptible
  x['healthy'] <- 1 - pars$prop.susceptible

  for (i in seq_len(n)) {
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
       fast.cancers = fast.cancers)
}

# Clinical incidence per 100,000 person-years of the undiagnosed, and the share
# of the undiagnosed with a precursor lesion, averaged over the ages of a group.
cohort.outputs <- function(cohort) {
  at.risk <- rowSums(cohort$trace[, UNDIAGNOSED, drop = FALSE])
  with.lesion <- rowSums(cohort$trace[, PRECURSORS, drop = FALSE]) / at.risk
  rows <- function(ages) match(ages, AGES)
  list(
    `Clinical incidence` = sapply(INCIDENCE.GROUPS, function(ages) {
      j <- rows(ages)
      sum(cohort$diagnoses[j]) / sum(at.risk[j]) * 1e5
    }),
    `Precursor prevalence` = sapply(PREVALENCE.GROUPS, function(ages) mean(with.lesion[rows(ages)]))
  )
}

simulate <- function(strategies, pars) {
  cohort <- run.cohort(pars)
  outputs <- cohort.outputs(cohort)

  discount <- (1 + pars$discount)^-(AGES - START.AGE)
  alive.undiagnosed <- rowSums(cohort$trace[, UNDIAGNOSED, drop = FALSE])
  alive.diagnosed <- cohort$trace[, 'clinical']
  cost <- sum(cohort$diagnoses * pars$cost.diagnosis * discount)
  qalys <- sum((alive.undiagnosed + alive.diagnosed * pars$utility.after.diagnosis) * discount)

  cancers <- sum(cohort$slow.cancers + cohort$fast.cancers)
  share <- c(slow = sum(cohort$slow.cancers), fast = sum(cohort$fast.cancers)) /
    (if (cancers > 0) cancers else 1)

  list(
    summary = data.frame(strategy = strategies, C = cost, E = qalys),
    outputs = setNames(rep(list(outputs), length(strategies)), strategies),
    pathway.share = setNames(rep(list(share), length(strategies)), strategies),
    cohort.info = setNames(rep(list(cohort$trace), length(strategies)), strategies)
  )
}
