# Cohort model of a precancerous high-grade lesion, annual cycles from age 20
# to 85. The cohort starts lesion-free; each year a lesion-free person may
# develop a lesion, a lesion may regress or become a preclinical cancer, and a
# preclinical cancer is diagnosed once it causes symptoms. Everyone alive also
# faces other-cause mortality.

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
run.cohort <- function(pars) {
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
  h <- 1; l <- 0; p <- 0; d <- 0
  for (i in seq_len(n)) {
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
       diagnosed = diagnosed, diagnoses = diagnoses)
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

# Discounted costs of diagnosed cancers and discounted QALYs of the cohort.
cohort.summary <- function(trace, pars) {
  discount <- (1 + pars$discount)^-(AGES - 20)
  alive.undiagnosed <- trace$healthy + trace$lesion + trace$preclinical
  costs <- sum(trace$diagnoses * pars$cost.diagnosis * discount)
  qalys <- sum((alive.undiagnosed + trace$diagnosed * pars$utility.after.diagnosis) * discount)
  c(C = costs, E = qalys)
}

simulate <- function(strategies, pars) {
  outputs <- list()
  summary <- data.frame(strategy = character(0), C = numeric(0), E = numeric(0))
  for (strategy in strategies) {
    trace <- run.cohort(pars)
    outputs[[strategy]] <- cohort.outputs(trace, pars)
    ce <- cohort.summary(trace, pars)
    summary <- rbind(summary, data.frame(strategy = strategy, C = ce[['C']], E = ce[['E']]))
  }
  list(summary = summary, outputs = outputs)
}
