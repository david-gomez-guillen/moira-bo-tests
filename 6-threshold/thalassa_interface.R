source('model.R')

MODEL.STRATA <- AGE.GROUPS

get.overview <- function() {
  paste(readLines('overview.md'), collapse = '\n')
}

get.strata <- function() {
  MODEL.STRATA
}

get.strategies <- function() {
  list(
    list(name = 'current_practice', display.name = 'Current practice',
         description = 'Transmission at its endemic equilibrium under current testing and treatment.')
  )
}

get.model.settings <- function() {
  list(cost.unit = '€', effectiveness.unit = 'QALY')
}

# Base values: a prior from the behavioural and clinical literature, one value
# per age group.
BASE.VALUES <- list(
  p.susceptibility = c(0.10, 0.08, 0.07, 0.06),
  activity.level = c(0.48, 0.40, 0.32, 0.25),
  r.clearance = c(0.50, 0.45, 0.40, 0.40),
  r.waning = c(0.15, 0.15, 0.12, 0.10)
)

stratified.parameter <- function(name, display.name, values, class) {
  lapply(seq_along(MODEL.STRATA), function(i) {
    list(name = name,
         display.name = display.name,
         base.value = values[i],
         stratum = MODEL.STRATA[i],
         class = class,
         distribution = 'beta',
         min.value = 0,
         max.value = 1)
  })
}

MODEL.PARAMETERS <- c(
  stratified.parameter('p.susceptibility',
                       'Per-act probability of acquiring the infection from an infected partner (women; men are 0.6 times as susceptible)',
                       BASE.VALUES$p.susceptibility, 'Transmission'),
  stratified.parameter('activity.level',
                       'Sexual activity, as a fraction of the full activity of the partnership matrix',
                       BASE.VALUES$activity.level, 'Sexual behaviour'),
  stratified.parameter('r.clearance',
                       'Annual rate of natural clearance of an infection, besides treatment',
                       BASE.VALUES$r.clearance, 'Natural history'),
  stratified.parameter('r.waning',
                       'Annual rate of losing antibodies and immunity after an infection',
                       BASE.VALUES$r.waning, 'Natural history')
)

get.parameters <- function() {
  MODEL.PARAMETERS
}

default.parameters <- function() {
  pars <- list()
  for (descriptor in MODEL.PARAMETERS) {
    if (is.null(pars[[descriptor$name]])) pars[[descriptor$name]] <- list()
    pars[[descriptor$name]][[descriptor$stratum]] <- descriptor$base.value
  }
  pars
}

run.simulation <- function(strategies, pars) {
  simulate(strategies, pars)
}

# Prevalence and seroprevalence surveys by age group and sex, one row per
# value: the series, the age group and the value.
read.targets <- function(file = 'data/targets.csv') {
  rows <- read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
  target <- list()
  for (series in unique(rows$series)) {
    block <- rows[rows$series == series, ]
    target[[series]] <- as.list(setNames(block$value, block$label))
  }
  target
}

CALIB.TARGET <- read.targets()

CALIB.PARAMS <- c('p.susceptibility', 'activity.level', 'r.clearance', 'r.waning')

# Relative squared error, averaged within each series and summed over series,
# so that the prevalences and the seroprevalences weigh alike. A simulation in
# which the infection dies out has every output at 0 and an error of one per
# series.
calibration.output.error <- function(output, target) {
  error <- 0
  for (series in names(target)) {
    target.values <- unlist(target[[series]])
    simulated <- unlist(output[[series]])[names(target.values)]
    error <- error + mean(((simulated - target.values) / target.values)^2)
  }
  error
}

calibration.error <- function(pars, target) {
  strategy <- 'current_practice'
  tryCatch({
    simulated <- run.simulation(strategy, pars)$outputs[[strategy]]
    output <- list()
    for (series in names(target))
      output[[series]] <- simulated[[series]][names(unlist(target[[series]]))]
    if (!all(is.finite(unlist(output)))) stop('non-finite output')
    list(error = calibration.output.error(output, target), output = output)
  }, error = function(e) {
    output <- lapply(target, function(series) {
      setNames(rep(NA_real_, length(series)), names(series))
    })
    list(error = Inf, output = output)
  })
}

# The initial guess scaled by a uniform factor per value, within [0, 1]. The
# calibration vector is laid out parameter by parameter, and within a
# parameter stratum by stratum.
generate.training.dataset <- function(initial_guess, n, ...) {
  variation <- list(...)$variation
  if (is.null(variation)) variation <- 1
  k <- length(initial_guess)
  factors <- matrix(runif(n * k, min = max(0, 1 - variation), max = 1 + variation), nrow = n)
  pmin(pmax(sweep(factors, 2, initial_guess, `*`), 0), 1)
}

# The constrained scheme: the infection is known to be endemic, so the
# parameters must give it a basic reproduction number of at least MIN.R0.
MIN.R0 <- 1.05

# The calibration vector holds the four parameters as blocks of one value per
# age group, in the order of CALIB.PARAMS. R0 needs every age group, so a
# vector covering only the age groups calibrated so far (`partial`) answers as
# satisfied. In the c(x) <= 0 convention.
calibration.constraints <- function(x, partial = FALSE) {
  x <- as.numeric(x)
  n <- length(MODEL.STRATA)
  if (length(x) < length(CALIB.PARAMS) * n) return(-1)
  block <- function(j) x[((j - 1) * n + 1):(j * n)]
  r0 <- tryCatch(basic.reproduction.number(block(1), block(2), block(3)),
                 error = function(e) NA_real_)
  if (!is.finite(r0)) return(1)
  MIN.R0 - r0
}

calibration.constraints.description <- function() {
  n <- length(MODEL.STRATA)
  paste0(
    'The vector holds ', paste(CALIB.PARAMS, collapse = ', '), ', in that order, each ',
    'as a block of ', n, ' values over the age groups ', paste(MODEL.STRATA, collapse = ', '), '.\n',
    'The infection is endemic, so the basic reproduction number R0 must be at least ', MIN.R0,
    '. R0 is the spectral radius of the next-generation matrix at the infection-free state, ',
    'computed by basic.reproduction.number() in model.R from p.susceptibility, activity.level ',
    'and r.clearance of every age group; r.waning does not enter it. Below R0 = 1 the infection ',
    'dies out and every prevalence is 0.'
  )
}

# Same as generate.training.dataset, keeping only the rows that satisfy the
# constraint.
generate.constrained.training.dataset <- function(initial_guess, n, ...) {
  rows <- matrix(NA_real_, nrow = 0, ncol = length(initial_guess))
  attempts <- 0
  while (nrow(rows) < n && attempts < 100) {
    attempts <- attempts + 1
    candidates <- generate.training.dataset(initial_guess, max(100, ceiling(1.5 * (n - nrow(rows)))), ...)
    feasible <- apply(candidates, 1, function(x) calibration.constraints(x) <= 0)
    rows <- rbind(rows, candidates[feasible, , drop = FALSE])
  }
  rows[seq_len(min(n, nrow(rows))), , drop = FALSE]
}

# Simulated against surveyed prevalence, by age group and sex. The plot is
# given the values of the evaluation, the simulated outputs named
# 'sim_<series>.<age group>'.
plot.prevalence.fit <- function(data) {
  if (!requireNamespace('ggplot2', quietly = TRUE)) return(NULL)
  rows <- list()
  for (sex in c('women', 'men')) {
    series <- sprintf('Prevalence (%s)', sex)
    prefix <- paste0('sim_', series, '.')
    simulated <- data[startsWith(names(data), prefix)]
    if (length(simulated) == 0) next
    groups <- sub(prefix, '', names(simulated), fixed = TRUE)
    observed <- unlist(CALIB.TARGET[[series]])[groups]
    rows[[sex]] <- data.frame(group = rep(groups, 2), sex = sex,
                              value = c(as.numeric(simulated), as.numeric(observed)),
                              source = rep(c('Simulated', 'Survey'), each = length(groups)))
  }
  if (length(rows) == 0) return(NULL)
  df <- do.call(rbind, rows)
  df$group <- factor(df$group, levels = MODEL.STRATA)
  df$line <- paste(df$sex, df$source)
  df$text <- sprintf('%s, %s, %s: %.2f%%', df$group, df$sex, df$source, 100 * df$value)
  ggplot2::ggplot(df, ggplot2::aes(x = group, y = 100 * value, colour = sex,
                                   linetype = source, group = line, text = text)) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::labs(x = 'Age group', y = 'Prevalence (%)', colour = '', linetype = '') +
    ggplot2::theme_minimal()
}

initial.guess <- function() {
  base <- default.parameters()
  unname(unlist(lapply(CALIB.PARAMS, function(param) unlist(base[[param]])[MODEL.STRATA])))
}

get.calibration.schemes <- function() {
  list(
    standard = list(
      description = 'Transmission and natural history by age group',
      parameters = CALIB.PARAMS,
      target = CALIB.TARGET,
      strata = MODEL.STRATA,
      initial_guess = initial.guess(),
      error_function = calibration.error,
      output_error = calibration.output.error,
      latent_space_training_set = generate.training.dataset,
      evaluation.plots = list(`Prevalence by age and sex` = plot.prevalence.fit)
    ),
    constrained = list(
      description = 'Transmission and natural history by age group, endemic',
      parameters = CALIB.PARAMS,
      target = CALIB.TARGET,
      strata = MODEL.STRATA,
      initial_guess = initial.guess(),
      error_function = calibration.error,
      output_error = calibration.output.error,
      constraints = calibration.constraints,
      constraints_description = calibration.constraints.description,
      latent_space_training_set = generate.constrained.training.dataset,
      evaluation.plots = list(`Prevalence by age and sex` = plot.prevalence.fit)
    )
  )
}
