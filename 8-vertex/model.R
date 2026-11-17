# Cohort model of a precancerous high-grade lesion, annual cycles from age 20
# to 85. The cohort starts lesion-free; each year a lesion-free person may
# develop a lesion, a lesion may regress or become a preclinical cancer, and a
# preclinical cancer is diagnosed once it causes symptoms. Everyone alive also
# faces other-cause mortality.
#
# Strategies: no intervention, which the calibration is run on, or a screening
# test at a fixed interval, which removes the lesions it finds and diagnoses
# the preclinical cancers it finds.

MODEL.STRATA <- c('20-34', '35-49', '50-64', '65-85')
AGES <- 20:85
STRATUM.OF.AGE <- findInterval(AGES, c(20, 35, 50, 65))
YEARS.IN.STRATUM <- tabulate(STRATUM.OF.AGE, nbins = length(MODEL.STRATA))

# Age groups the cancer registry reports incidence in.
INCIDENCE.GROUPS <- list(`20-29` = 20:29, `30-34` = 30:34, `35-39` = 35:39,
                         `40-49` = 40:49, `50-59` = 50:59, `60-64` = 60:64,
                         `65-74` = 65:74, `75-85` = 75:85)

# Interval of the follow-up study of untreated lesions, in years.
FOLLOWUP.YEARS <- 2

# The ages at which each strategy screens the cohort.
SCREENING.AGES <- list(
  no_intervention = integer(0),
  screening_3y = seq(25, 64, by = 3),
  screening_5y = seq(30, 64, by = 5)
)

# One value per year of age, from a parameter given by stratum.
by.age <- function(value) {
  as.numeric(unlist(value)[MODEL.STRATA])[STRATUM.OF.AGE]
}

by.stratum <- function(value) {
  as.numeric(unlist(value)[MODEL.STRATA])
}

other.mortality <- function(pars) {
  pmin(1, pars$mortality.other.base * exp(pars$mortality.other.rate * (AGES - 20)))
}

# Runs the cohort and returns, for every year of age, the people at the start
# of the year by state and the cancers diagnosed during it.
#
# At each age of `screening.ages` the year starts with a screening test. Of
# those who attend, it finds a lesion with one sensitivity and a preclinical
# cancer with another. A lesion it finds is removed, which leaves the person
# lesion-free; a cancer it finds is diagnosed, and counts among the diagnoses
# of that year. Returns too, per year, the people screened and the lesions
# removed.
run.cohort <- function(pars, screening.ages = integer(0)) {
  onset <- by.age(pars$p.lesion.onset)
  regress <- by.age(pars$p.lesion.regress)
  progress <- by.age(pars$p.lesion.progress)
  symptoms <- pars$p.symptoms
  cancer.death <- pars$p.cancer.death
  survival <- 1 - other.mortality(pars)
  if (any(regress + progress > 1))
    stop('A lesion cannot regress and progress with a total probability above 1')

  n <- length(AGES)
  healthy <- lesion <- preclinical <- diagnosed <- diagnoses <- numeric(n)
  screened <- removed <- screen.detected <- numeric(n)
  h <- 1; l <- 0; p <- 0; d <- 0
  for (i in seq_len(n)) {
    if (AGES[i] %in% screening.ages) {
      attending <- pars$screening.adherence
      screened[i] <- attending * (h + l + p)
      removed[i] <- l * attending * pars$screening.sens.lesion
      screen.detected[i] <- p * attending * pars$screening.sens.cancer
      h <- h + removed[i]; l <- l - removed[i]
      p <- p - screen.detected[i]; d <- d + screen.detected[i]
    }
    healthy[i] <- h; lesion[i] <- l; preclinical[i] <- p; diagnosed[i] <- d
    new.diagnoses <- p * symptoms
    diagnoses[i] <- new.diagnoses
    h.next <- h * (1 - onset[i]) + l * regress[i]
    l.next <- h * onset[i] + l * (1 - regress[i] - progress[i])
    p.next <- l * progress[i] + p * (1 - symptoms)
    d.next <- d * (1 - cancer.death) + new.diagnoses
    h <- h.next * survival[i]; l <- l.next * survival[i]
    p <- p.next * survival[i]; d <- d.next * survival[i]
  }
  list(healthy = healthy, lesion = lesion, preclinical = preclinical,
       diagnosed = diagnosed, diagnoses = diagnoses + screen.detected,
       screened = screened, removed = removed, symptomatic = diagnoses,
       screen.detected = screen.detected)
}

# Registry incidence per 100,000 person-years among the undiagnosed, lesion
# prevalence as a population examination study would find it, and the share of
# lesions still present at the end of the follow-up of untreated lesions.
cohort.outputs <- function(trace, pars) {
  at.risk <- trace$healthy + trace$lesion + trace$preclinical
  incidence <- sapply(INCIDENCE.GROUPS, function(ages) {
    k <- ages - 19
    1e5 * sum(trace$diagnoses[k]) / sum(at.risk[k])
  })
  prevalence <- setNames(sapply(seq_along(MODEL.STRATA), function(j) {
    k <- STRATUM.OF.AGE == j
    sum(trace$lesion[k]) / sum(at.risk[k])
  }), MODEL.STRATA)
  # A lesion found at the first visit is still there at the last one if it
  # neither regressed nor progressed in any year between them.
  persistence <- setNames((1 - by.stratum(pars$p.lesion.regress) -
                             by.stratum(pars$p.lesion.progress))^FOLLOWUP.YEARS,
                          MODEL.STRATA)
  list(`Cancer incidence` = incidence,
       `Lesion prevalence` = prevalence,
       `Lesion persistence` = persistence)
}

# Discounted costs of diagnosed cancers, of the screening tests and of the
# lesions they remove, and discounted QALYs of the cohort.
cohort.summary <- function(trace, pars) {
  discount <- (1 + pars$discount)^-(AGES - 20)
  alive.undiagnosed <- trace$healthy + trace$lesion + trace$preclinical
  costs <- sum(trace$diagnoses * pars$cost.diagnosis * discount) +
    sum((trace$screened * pars$cost.screening + trace$removed * pars$cost.lesion.removal) *
          discount)
  qalys <- sum((alive.undiagnosed + trace$diagnosed * pars$utility.after.diagnosis) * discount)
  c(C = costs, E = qalys)
}

# The plots of the Base case tab for a cohort, each under the title of its tab:
# a numeric vector named by age, or a data frame of the ages and a line per
# column. Rates and shares are among the undiagnosed, as the outputs are. The
# cumulative diagnoses at 85 are the lifetime risk the constraints bound.
cohort.plots <- function(trace) {
  at.risk <- trace$healthy + trace$lesion + trace$preclinical
  incidence <- data.frame(Age = AGES, `After symptoms` = trace$symptomatic / at.risk * 1e5,
                          check.names = FALSE)
  if (any(trace$screen.detected > 0))
    incidence$`At screening` <- trace$screen.detected / at.risk * 1e5
  list(
    `Cancer incidence per 100,000` = incidence,
    `Lesion prevalence` = setNames(trace$lesion / at.risk, AGES),
    `Preclinical cancer prevalence per 100,000` =
      setNames(trace$preclinical / at.risk * 1e5, AGES),
    `Cumulative diagnoses per 100,000` = setNames(cumsum(trace$diagnoses) * 1e5, AGES)
  )
}

# The plots are drawn only when `make.plots` asks for them, so that every other
# run, the calibration's among them, returns what it did without them.
simulate <- function(strategies, pars, make.plots = FALSE) {
  outputs <- plots <- list()
  summary <- data.frame(strategy = character(0), C = numeric(0), E = numeric(0))
  for (strategy in strategies) {
    ages <- SCREENING.AGES[[strategy]]
    if (is.null(ages)) stop(sprintf("Unknown strategy '%s'", strategy))
    trace <- run.cohort(pars, ages)
    outputs[[strategy]] <- cohort.outputs(trace, pars)
    ce <- cohort.summary(trace, pars)
    summary <- rbind(summary, data.frame(strategy = strategy, C = ce[['C']], E = ce[['E']]))
    if (make.plots) plots[[strategy]] <- cohort.plots(trace)
  }
  results <- list(summary = summary, outputs = outputs)
  if (make.plots) results$plots <- plots
  results
}
