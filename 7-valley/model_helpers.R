# Helpers for model.R: bookkeeping that the model needs but that is not needed
# to understand how it works. They turn age-dependent parameters into the values
# of one age, and summarize the results of a run per age group.


# ---- Age groups -----------------------------------------------------------------

# Results by age are reported per ten-year age group, the strata of the model.
MODEL.STRATA <- c('20-29', '30-39', '40-49', '50-59', '60-69', '70-79', '80-89')
STRATUM.SIZE <- 10

# The position of the stratum an age falls in.
stratum.index <- function(age) {
  min(length(MODEL.STRATA), max(1, (age - MODEL.AGE.START) %/% STRATUM.SIZE + 1))
}

stratum.of.age <- function(age) {
  MODEL.STRATA[stratum.index(age)]
}


# ---- Parameters -----------------------------------------------------------------

# A parameter is either one number for every age or one value per stratum.
# This gives the value that applies at an age, whichever of the two it is.
par.at.age <- function(value, age) {
  if (!is.list(value) && length(value) == 1) return(as.numeric(value))

  stratum <- stratum.of.age(age)
  if (stratum %in% names(value)) return(as.numeric(value[[stratum]]))

  # Values without stratum names are taken in the order of the strata.
  as.numeric(value[[min(length(value), stratum.index(age))]])
}

# Every parameter, as the plain number it is worth at an age. The rest of the
# model works on these and never has to think about strata.
parameters.at.age <- function(pars, age) {
  lapply(pars, par.at.age, age = age)
}


# ---- States ---------------------------------------------------------------------

# One value per state, all zero, to be filled in by name.
per.state <- function() setNames(rep(0, length(MODEL.STATES)), MODEL.STATES)


# ---- Results --------------------------------------------------------------------

# The outputs a calibration compares with observed data, per stratum, from what
# simulate() measured year by year. Incidence and prevalence are yearly rates
# over the people at risk, averaged over the years of the stratum. The
# late-stage share is the part of the cancers diagnosed in the stratum that
# were already advanced.
calibration.outputs <- function(ages, new.early, new.late, at.risk, with.lesion) {
  stratum <- factor(sapply(ages, stratum.of.age), levels = MODEL.STRATA)
  anyone.at.risk <- at.risk > 1e-12
  yearly.incidence <- ifelse(anyone.at.risk, (new.early + new.late) / at.risk, 0)
  yearly.prevalence <- ifelse(anyone.at.risk, with.lesion / at.risk, 0)

  diagnosed <- tapply(new.early + new.late, stratum, sum)
  diagnosed.late <- tapply(new.late, stratum, sum)
  late.share <- ifelse(diagnosed > 1e-12, diagnosed.late / diagnosed, 0)

  list(
    `Cancer incidence` = setNames(as.numeric(tapply(yearly.incidence, stratum, mean)), MODEL.STRATA),
    `Lesion prevalence` = setNames(as.numeric(tapply(yearly.prevalence, stratum, mean)), MODEL.STRATA),
    `Late-stage share` = setNames(as.numeric(late.share), MODEL.STRATA)
  )
}
