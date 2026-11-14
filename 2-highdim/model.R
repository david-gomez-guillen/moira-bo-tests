# Natural history of an HPV-related cancer, followed as two closed cohorts, women
# and men, from age 20 to 85 in annual cycles.
#
# States: healthy, persistent low-grade lesion, high-grade lesion, preclinical
# cancer, diagnosed cancer, and death from other causes. Every transition but
# symptomatic diagnosis and other-cause death takes one value per stratum, a
# stratum being a sex and an age group.

SEXES <- c('F', 'M')
AGE.GROUPS <- c('20_29', '30_39', '40_49', '50_59', '60_69', '70_85')
AGE.GROUP.START <- c(20, 30, 40, 50, 60, 70)
MODEL.STRATA <- as.vector(t(outer(SEXES, AGE.GROUPS, paste, sep='.')))
FIRST.AGE <- 20
LAST.AGE <- 84

CALIBRATED.TRANSITIONS <- c('p.lesion.onset', 'p.lesion.clearance', 'p.lesion.progression',
                            'p.high.grade.regression', 'p.invasion')

age.group.of <- function(age) AGE.GROUPS[findInterval(age, AGE.GROUP.START)]

# Other-cause mortality, Gompertz in age.
other.mortality <- function(pars, sex, age) {
  base <- if (sex == 'F') pars$mortality.base.F else pars$mortality.base.M
  pmin(1, base * exp(pars$mortality.gompertz.rate * (age - FIRST.AGE)))
}

# Value of a transition for a stratum, whether the parameter arrives by stratum
# or as a single value.
stratum.value <- function(pars, name, stratum) {
  value <- pars[[name]]
  if (is.list(value)) value <- value[[stratum]]
  as.numeric(value)
}

# The probabilities of leaving a state in one cycle, and what is left to stay.
# Exits that add up to more than the whole state cannot be applied as they are:
# by default they are scaled down to share it, and with strict = TRUE the
# simulation stops instead.
state.exits <- function(exits, state, sex, age, strict) {
  total <- sum(exits)
  if (total > 1) {
    if (strict)
      stop(sprintf('Transition probabilities out of the %s state add up to %.4f at age %d (%s)',
                   state, total, age, sex), call.=FALSE)
    exits <- exits / total
    total <- 1
  }
  c(exits, stay=1 - total)
}

simulate.cohort <- function(pars, sex, strict) {
  n.ages <- LAST.AGE - FIRST.AGE + 1
  trace <- matrix(0, n.ages, 6,
                  dimnames=list(FIRST.AGE:LAST.AGE,
                                c('healthy', 'low.grade', 'high.grade', 'preclinical',
                                  'diagnosed', 'dead.other')))
  new.diagnoses <- numeric(n.ages)
  state <- c(healthy=1, low.grade=0, high.grade=0, preclinical=0, diagnosed=0, dead.other=0)

  for (i in seq_len(n.ages)) {
    age <- FIRST.AGE + i - 1
    stratum <- paste(sex, age.group.of(age), sep='.')
    trace[i, ] <- state
    m <- other.mortality(pars, sex, age)
    p <- function(name) stratum.value(pars, name, stratum)

    healthy <- state.exits(c(onset=p('p.lesion.onset'), death=m),
                           'healthy', sex, age, strict)
    low <- state.exits(c(clearance=p('p.lesion.clearance'),
                         progression=p('p.lesion.progression'), death=m),
                       'low-grade lesion', sex, age, strict)
    high <- state.exits(c(regression=p('p.high.grade.regression'),
                          invasion=p('p.invasion'), death=m),
                        'high-grade lesion', sex, age, strict)
    pre <- state.exits(c(diagnosis=pars$p.symptomatic.diagnosis, death=m),
                       'preclinical cancer', sex, age, strict)

    new.diagnoses[i] <- state[['preclinical']] * pre[['diagnosis']]
    state <- c(
      healthy=state[['healthy']] * healthy[['stay']] + state[['low.grade']] * low[['clearance']],
      low.grade=state[['healthy']] * healthy[['onset']] + state[['low.grade']] * low[['stay']] +
        state[['high.grade']] * high[['regression']],
      high.grade=state[['low.grade']] * low[['progression']] + state[['high.grade']] * high[['stay']],
      preclinical=state[['high.grade']] * high[['invasion']] + state[['preclinical']] * pre[['stay']],
      diagnosed=state[['diagnosed']] + new.diagnoses[i],
      dead.other=state[['dead.other']] + m * (state[['healthy']] + state[['low.grade']] +
                                                state[['high.grade']] + state[['preclinical']])
    )
  }
  list(trace=trace, new.diagnoses=new.diagnoses)
}

# Outputs by stratum: the incidence of diagnosed cancer per person-year among
# those alive and never diagnosed, and the prevalence of lesions of either grade
# among them, averaged over the years of the age group.
summarise.cohort <- function(cohort, sex) {
  trace <- cohort$trace
  at.risk <- rowSums(trace[, c('healthy', 'low.grade', 'high.grade', 'preclinical')])
  groups <- age.group.of(as.numeric(rownames(trace)))
  lesions <- (trace[, 'low.grade'] + trace[, 'high.grade']) / at.risk
  incidence <- tapply(cohort$new.diagnoses, groups, sum) / tapply(at.risk, groups, sum)
  prevalence <- tapply(lesions, groups, mean)
  names(incidence) <- names(prevalence) <- paste(sex, names(incidence), sep='.')
  list(incidence=incidence[paste(sex, AGE.GROUPS, sep='.')],
       prevalence=prevalence[paste(sex, AGE.GROUPS, sep='.')],
       life.years=sum(rowSums(trace[, c('healthy', 'low.grade', 'high.grade',
                                        'preclinical', 'diagnosed')])),
       cancers=sum(cohort$new.diagnoses))
}

simulate <- function(pars, strict=FALSE) {
  cohorts <- lapply(SEXES, function(sex) {
    cohort <- simulate.cohort(pars, sex, strict)
    c(summarise.cohort(cohort, sex), list(trace=cohort$trace))
  })
  names(cohorts) <- SEXES
  list(
    outputs=list(
      `Cancer incidence`=unlist(unname(lapply(cohorts, `[[`, 'incidence'))),
      `Lesion prevalence`=unlist(unname(lapply(cohorts, `[[`, 'prevalence')))
    ),
    life.years=mean(sapply(cohorts, `[[`, 'life.years')),
    cancers=mean(sapply(cohorts, `[[`, 'cancers')),
    traces=lapply(cohorts, `[[`, 'trace')
  )
}
