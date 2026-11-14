# Age- and sex-structured SIRS model of a curable sexually transmitted
# infection with partial immunity that wanes. Heterosexual transmission only:
# women acquire the infection from men and men from women.
#
# State, by sex and age group, as fractions of the group: susceptible (s),
# infected (i) and recovered with antibodies (r). The population is stationary:
# people enter at 15 susceptible, move up an age group at the rate of
# AGE.EXIT and leave the sexually active population after 64.

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
transmission.matrices <- function(susceptibility, activity, clearance) {
  ends.women <- clearance + TREATMENT.RATE[['women']]
  ends.men <- clearance + TREATMENT.RATE[['men']]
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
basic.reproduction.number <- function(susceptibility, activity, clearance) {
  k <- transmission.matrices(susceptibility, activity, clearance)
  # Expected time spent infected in each age group, from an infection in each
  # age group, counting the time carried into the older groups by ageing.
  time.infected <- function(ending) {
    v <- diag(ending + AGE.EXIT)
    v[cbind(2:4, 1:3)] <- -AGE.EXIT[2:4]
    solve(v)
  }
  from.men <- k$women %*% time.infected(clearance + TREATMENT.RATE[['men']])
  from.women <- k$men %*% time.infected(clearance + TREATMENT.RATE[['women']])
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

# Integrates the model to its equilibrium and returns the prevalence and
# seroprevalence of each age group and sex, with the trajectory of the
# prevalences when asked for it.
simulate.equilibrium <- function(susceptibility, activity, clearance, waning,
                                 dt = TIME.STEP, horizon = HORIZON, trace = FALSE) {
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
  named <- function(x) setNames(x, AGE.GROUPS)
  out <- list(
    `Prevalence (women)` = named(women$i),
    `Prevalence (men)` = named(men$i),
    `Seroprevalence (women)` = named(women$i + women$r),
    `Seroprevalence (men)` = named(men$i + men$r)
  )
  if (trace) attr(out, 'trace') <- history
  out
}

# The model run on a parameter list, whose four calibrated parameters hold a
# value per age group.
simulate <- function(strategies, pars) {
  by.age <- function(name) as.numeric(unlist(pars[[name]])[AGE.GROUPS])
  outputs <- simulate.equilibrium(by.age('p.susceptibility'), by.age('activity.level'),
                                  by.age('r.clearance'), by.age('r.waning'))
  list(outputs = setNames(rep(list(outputs), length(strategies)), strategies))
}
