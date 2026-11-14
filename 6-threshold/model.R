# Age- and sex-structured SIRS model of a curable sexually transmitted
# infection with partial immunity that wanes. Heterosexual transmission only:
# women acquire the infection from men and men from women.
#
# State, by sex and age group, as fractions of the group: susceptible (s),
# infected (i) and recovered with antibodies (r). The population is stationary:
# people enter at 15 susceptible, move up an age group at the rate of
# AGE.EXIT and leave the sexually active population after 64.
#
# Strategies: current practice, the endemic equilibrium the calibration is run
# on, or a programme that tests young people for the infection and treats
# those found infected, started at that equilibrium and followed for
# PROGRAMME.YEARS.

AGE.GROUPS <- c('15-24', '25-34', '35-44', '45-64')

# Rate of leaving each age group, the inverse of its width in years.
AGE.EXIT <- c(1/10, 1/10, 1/10, 1/20)

# Relative size of each age group in the stationary population, the same for
# both sexes.
GROUP.SIZE <- c(1, 1, 1, 2)

# Partnerships formed per year between men (rows) and women (columns) of each
# age group, at full sexual activity of both, per unit of population. Most
# partnerships are within an age group, and men's partners tend to be younger
# than they are.
PARTNERSHIP.MATRIX <- matrix(
  c(14.0, 2.4, 0.5, 0.2,
     7.0, 9.0, 1.4, 0.5,
     2.0, 4.5, 6.0, 1.4,
     0.7, 1.4, 3.5, 7.0),
  nrow = 4, byrow = TRUE, dimnames = list(men = AGE.GROUPS, women = AGE.GROUPS))

# Sex acts per year within an ongoing partnership.
ACTS.PER.YEAR <- 60

# Susceptibility of men per act relative to that of women of the same age.
RELATIVE.SUSCEPTIBILITY.MEN <- 0.6

# Annual rate at which an infection is treated, on top of its natural
# clearance: men are more often symptomatic and seek treatment sooner.
TREATMENT.RATE <- c(women = 0.15, men = 0.35)

# Share of every group infected when the simulation starts.
SEED.PREVALENCE <- 0.01

# Forward-Euler integration: step and horizon in years. The horizon is long
# enough for the slowest group, the oldest, to settle on its equilibrium.
TIME.STEP <- 0.05
HORIZON <- 300

# Years a testing programme is followed for, from the endemic equilibrium of
# current practice, and the age groups each strategy tests, by sex.
PROGRAMME.YEARS <- 20
SCREENED.GROUPS <- list(
  current_practice = list(women = character(0), men = character(0)),
  screening_young_women = list(women = '15-24', men = character(0)),
  screening_young = list(women = '15-24', men = '15-24')
)

# The rates at which each age group of a sex acquires the infection, per
# infected partner group, from the calibrated values by age group:
#   susceptibility  per-act probability of acquiring the infection (women);
#   activity        sexual activity, as a fraction of the full activity the
#                   partnership matrix describes;
#   clearance       annual rate at which an infection clears naturally.
#
# A partnership of a susceptible person with an infected one passes the
# infection on only if it does so before the partnership ends or the infection
# clears or is treated. Partnerships end when either partner moves on to a new
# one, so the more active both are, the shorter they last.
#
# `screening` holds, by sex, the annual rate at which a testing programme finds
# and treats the infections of each age group, on top of current practice.
NO.SCREENING <- list(women = 0, men = 0)
transmission.matrices <- function(susceptibility, activity, clearance, screening = NO.SCREENING) {
  ends.women <- clearance + TREATMENT.RATE[['women']] + screening$women
  ends.men <- clearance + TREATMENT.RATE[['men']] + screening$men
  partnerships <- PARTNERSHIP.MATRIX * outer(activity, activity)
  rate.men <- rowSums(partnerships) / GROUP.SIZE
  rate.women <- colSums(partnerships) / GROUP.SIZE
  partnership.ends <- outer(rate.men, rate.women, `+`)

  # Men (rows) acquiring from women (columns), and women (rows) acquiring
  # from men (columns): the rate at which the infection crosses an ongoing
  # partnership, against the rates at which the partnership or the infection
  # ends.
  crossing.men <- RELATIVE.SUSCEPTIBILITY.MEN * susceptibility * ACTS.PER.YEAR
  crossing.women <- susceptibility * ACTS.PER.YEAR
  pass.men <- crossing.men / (crossing.men + sweep(partnership.ends, 2, ends.women, `+`))
  pass.women <- crossing.women / (crossing.women + sweep(t(partnership.ends), 2, ends.men, `+`))

  list(men = partnerships * pass.men / GROUP.SIZE,
       women = t(partnerships) * pass.women / GROUP.SIZE)
}

# Basic reproduction number: the spectral radius of the next-generation matrix
# at the infection-free state. Transmission only runs between the sexes, so
# that matrix pairs a women-from-men block with a men-from-women one, and its
# spectral radius is the square root of that of their 4 x 4 product.
basic.reproduction.number <- function(susceptibility, activity, clearance,
                                      screening = NO.SCREENING) {
  k <- transmission.matrices(susceptibility, activity, clearance, screening)
  # Expected time spent infected in each age group, from an infection in each
  # age group, counting the time carried into the older groups by ageing.
  time.infected <- function(ending) {
    v <- diag(ending + AGE.EXIT)
    v[cbind(2:4, 1:3)] <- -AGE.EXIT[2:4]
    solve(v)
  }
  from.men <- k$women %*% time.infected(clearance + TREATMENT.RATE[['men']] + screening$men)
  from.women <- k$men %*% time.infected(clearance + TREATMENT.RATE[['women']] + screening$women)
  sqrt(max(Mod(eigen(from.men %*% from.women, only.values = TRUE)$values)))
}

# One Euler step of the three states of a sex, given its force of infection
# and the rate at which its infections end. Cleared and treated infections
# both leave antibodies and a period of immunity.
euler.step <- function(state, force, ending, waning, dt) {
  s <- state$s; i <- state$i; r <- state$r
  entering <- function(x, first) c(first, x[-length(x)])
  infections <- force * s
  recoveries <- ending * i
  losses <- waning * r
  list(
    s = s + dt * (AGE.EXIT * (entering(s, 1) - s) - infections + losses),
    i = i + dt * (AGE.EXIT * (entering(i, 0) - i) + infections - recoveries),
    r = r + dt * (AGE.EXIT * (entering(r, 0) - r) + recoveries - losses)
  )
}

# The prevalence and seroprevalence of each age group and sex in a state.
state.outputs <- function(women, men) {
  named <- function(x) setNames(x, AGE.GROUPS)
  list(
    `Prevalence (women)` = named(women$i),
    `Prevalence (men)` = named(men$i),
    `Seroprevalence (women)` = named(women$i + women$r),
    `Seroprevalence (men)` = named(men$i + men$r)
  )
}

# Integrates the model to its equilibrium and returns the prevalence and
# seroprevalence of each age group and sex, with the trajectory of the
# prevalences when asked for it, and with the equilibrium itself, the states of
# both sexes, when asked for that.
simulate.equilibrium <- function(susceptibility, activity, clearance, waning,
                                 dt = TIME.STEP, horizon = HORIZON, trace = FALSE,
                                 equilibrium = FALSE) {
  k <- transmission.matrices(susceptibility, activity, clearance)
  start <- list(s = rep(1 - SEED.PREVALENCE, 4), i = rep(SEED.PREVALENCE, 4), r = rep(0, 4))
  women <- start
  men <- start
  n.steps <- round(horizon / dt)
  history <- if (trace) matrix(NA_real_, n.steps + 1, 8) else NULL
  if (trace) history[1, ] <- c(women$i, men$i)
  for (step in seq_len(n.steps)) {
    force.women <- as.vector(k$women %*% men$i)
    force.men <- as.vector(k$men %*% women$i)
    women <- euler.step(women, force.women, clearance + TREATMENT.RATE[['women']], waning, dt)
    men <- euler.step(men, force.men, clearance + TREATMENT.RATE[['men']], waning, dt)
    if (trace) history[step + 1, ] <- c(women$i, men$i)
  }
  out <- state.outputs(women, men)
  if (trace) attr(out, 'trace') <- history
  if (equilibrium) attr(out, 'equilibrium') <- list(women = women, men = men)
  out
}

# Follows a testing programme from a starting state, the equilibrium of
# current practice, for `years`. Returns the outputs at the end and, per year
# of the programme, what happened in the population of each sex, as a share of
# the whole population aged 15-64: the person-years infected, the new
# infections, the infections treated, by current practice or by the
# programme, and the tests the programme does. A test finds an infection with
# probability `sensitivity`, and those tested are tested at `rate` a year.
simulate.programme <- function(susceptibility, activity, clearance, waning, start, tested,
                               rate, sensitivity, years = PROGRAMME.YEARS, dt = TIME.STEP) {
  testing <- lapply(tested, function(groups) rate * (AGE.GROUPS %in% groups))
  screening <- lapply(testing, function(r) r * sensitivity)
  k <- transmission.matrices(susceptibility, activity, clearance, screening)
  ending <- list(women = clearance + TREATMENT.RATE[['women']] + screening$women,
                 men = clearance + TREATMENT.RATE[['men']] + screening$men)
  weight <- GROUP.SIZE / (2 * sum(GROUP.SIZE))

  women <- start$women
  men <- start$men
  steps.per.year <- round(1 / dt)
  yearly <- matrix(0, years, 8, dimnames = list(NULL, c(
    'infected.women', 'infected.men', 'infections.women', 'infections.men',
    'treated.women', 'treated.men', 'tests.women', 'tests.men')))
  for (step in seq_len(years * steps.per.year)) {
    year <- (step - 1) %/% steps.per.year + 1
    force <- list(women = as.vector(k$women %*% men$i), men = as.vector(k$men %*% women$i))
    for (sex in c('women', 'men')) {
      x <- if (sex == 'women') women else men
      treated <- (TREATMENT.RATE[[sex]] + screening[[sex]]) * x$i
      yearly[year, paste0(c('infected.', 'infections.', 'treated.', 'tests.'), sex)] <-
        yearly[year, paste0(c('infected.', 'infections.', 'treated.', 'tests.'), sex)] +
        dt * c(sum(weight * x$i), sum(weight * force[[sex]] * x$s),
               sum(weight * treated), sum(weight * testing[[sex]] * (x$s + x$i + x$r)))
    }
    women.next <- euler.step(women, force$women, ending$women, waning, dt)
    men <- euler.step(men, force$men, ending$men, waning, dt)
    women <- women.next
  }
  list(outputs = state.outputs(women, men), yearly = yearly,
       reproduction.number = basic.reproduction.number(susceptibility, activity, clearance,
                                                       screening))
}

# Costs and QALYs per person aged 15-64 over the years of a programme,
# discounted: the tests, the treatments and the sequelae of infection in women,
# which arise at a rate per year of infection, each with a cost, the sequelae
# with a loss of QALYs too. Every person-year is a QALY but for that loss.
programme.summary <- function(yearly, pars) {
  discount <- (1 + pars$discount)^-(seq_len(nrow(yearly)) - 0.5)
  sequelae <- yearly[, 'infected.women'] * pars$r.sequelae
  cost <- (yearly[, 'tests.women'] + yearly[, 'tests.men']) * pars$cost.test +
    (yearly[, 'treated.women'] + yearly[, 'treated.men']) * pars$cost.treatment +
    sequelae * pars$cost.sequelae
  c(C = sum(discount * cost), E = sum(discount * (1 - sequelae * pars$qaly.loss.sequelae)))
}

# The model run on a parameter list, whose four calibrated parameters hold a
# value per age group. The outputs of current practice are those of its
# endemic equilibrium, the outputs of a programme those at its end, and the
# costs and effects of every strategy are those of the years of a programme
# from that equilibrium.
simulate <- function(strategies, pars) {
  by.age <- function(name) as.numeric(unlist(pars[[name]])[AGE.GROUPS])
  susceptibility <- by.age('p.susceptibility')
  activity <- by.age('activity.level')
  clearance <- by.age('r.clearance')
  waning <- by.age('r.waning')
  endemic <- simulate.equilibrium(susceptibility, activity, clearance, waning,
                                  equilibrium = TRUE)
  start <- attr(endemic, 'equilibrium')
  attr(endemic, 'equilibrium') <- NULL

  summary <- data.frame(strategy = character(0), C = numeric(0), E = numeric(0))
  outputs <- yearly <- reproduction <- list()
  for (strategy in strategies) {
    tested <- SCREENED.GROUPS[[strategy]]
    if (is.null(tested)) stop(sprintf("Unknown strategy '%s'", strategy))
    programme <- simulate.programme(susceptibility, activity, clearance, waning, start, tested,
                                    rate = pars$screening.rate,
                                    sensitivity = pars$test.sensitivity)
    ce <- programme.summary(programme$yearly, pars)
    summary <- rbind(summary, data.frame(strategy = strategy, C = ce[['C']], E = ce[['E']]))
    outputs[[strategy]] <- if (strategy == 'current_practice') endemic else programme$outputs
    yearly[[strategy]] <- data.frame(year = seq_len(nrow(programme$yearly)), programme$yearly)
    reproduction[[strategy]] <- programme$reproduction.number
  }
  list(summary = summary, outputs = outputs, programme = yearly,
       reproduction.number = reproduction)
}
