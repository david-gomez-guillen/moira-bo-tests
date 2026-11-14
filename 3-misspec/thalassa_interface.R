source('model.R')

get.overview <- function() {
  paste(readLines('overview.md'), collapse = '\n')
}

get.strategies <- function() {
  list(
    list(name = 'no_intervention', display.name = 'No intervention',
         description = 'Natural history only: cancers are diagnosed when they give symptoms.')
  )
}

get.strata <- function() {
  MODEL.STRATA
}

get.model.settings <- function() {
  list(cost.unit = '€', effectiveness.unit = 'QALY')
}

# A parameter with one value per age stratum, as one entry per stratum.
stratified.parameter <- function(name, display.name, values, class) {
  lapply(seq_along(MODEL.STRATA), function(i) {
    list(name = name, display.name = display.name, base.value = values[i],
         stratum = MODEL.STRATA[i], class = class, distribution = 'beta',
         min.value = 0, max.value = 1)
  })
}

probability <- function(name, display.name, base.value, class) {
  list(name = name, display.name = display.name, base.value = base.value,
       class = class, distribution = 'beta', min.value = 0, max.value = 1)
}

# Base values of the three age-dependent parameters, one per stratum, from the
# literature on lesion onset and on the natural history of comparable cancers.
BASE.LESION.ONSET <- c(0.008, 0.017, 0.019, 0.023)
BASE.LESION.PROGRESS <- c(0.018, 0.018, 0.015, 0.012)
BASE.STAGE.PROGRESS <- c(0.55, 0.65, 0.65, 0.55)

MODEL.PARAMETERS <- c(
  stratified.parameter('p.lesion.onset',
                       'Annual probability that a healthy person develops a precancerous lesion',
                       BASE.LESION.ONSET, 'Natural history'),
  stratified.parameter('p.lesion.progress',
                       'Annual probability that a lesion becomes a preclinical localized cancer',
                       BASE.LESION.PROGRESS, 'Natural history'),
  stratified.parameter('p.stage.progress',
                       'Annual probability that a preclinical localized cancer becomes advanced',
                       BASE.STAGE.PROGRESS, 'Cancer progression'),
  list(
    probability('p.lesion.regress', 'Annual probability that a lesion regresses',
                0.04, 'Natural history'),
    probability('p.symptoms.early',
                'Annual probability that a preclinical localized cancer gives symptoms and is diagnosed',
                0.15, 'Cancer progression'),
    probability('p.symptoms.late',
                'Annual probability that a preclinical advanced cancer gives symptoms and is diagnosed',
                0.55, 'Cancer progression'),
    probability('mortality.other.base', 'Annual probability of other-cause death at age 20',
                0.0005, 'General'),
    list(name = 'mortality.other.rate', display.name = 'Gompertz rate of other-cause mortality',
         base.value = 0.085, class = 'General', distribution = 'lognormal', min.value = 0),
    list(name = 'cost.diagnosis', display.name = 'Cost of diagnosing and treating a cancer',
         base.value = 25000, class = 'Costs', distribution = 'gamma',
         min.value = 0, max.value = Inf),
    probability('utility.after.diagnosis', 'Utility of a year lived after a cancer diagnosis',
                0.75, 'Utilities'),
    list(name = 'discount', display.name = 'Discount rate for costs and effects',
         base.value = 0.03, class = 'General', min.value = 0, max.value = 1)
  )
)

get.parameters <- function() {
  MODEL.PARAMETERS
}

# The base values in the shape simulate() expects: a named list, a stratified
# parameter being itself a list named by stratum.
default.parameters <- function() {
  pars <- list()
  for (p in MODEL.PARAMETERS) {
    if (is.null(p$stratum)) {
      pars[[p$name]] <- p$base.value
    } else {
      if (is.null(pars[[p$name]])) pars[[p$name]] <- list()
      pars[[p$name]][[p$stratum]] <- p$base.value
    }
  }
  pars
}

get.model.states <- function() {
  list(
    nodes = data.frame(
      id = c('healthy', 'lesion', 'pre.early', 'pre.late', 'diagnosed', 'dead.other'),
      label = c('Healthy', 'Lesion', 'Preclinical\nlocalized', 'Preclinical\nadvanced',
                'Diagnosed', 'Other-cause\ndeath'),
      description = c('No lesion. The whole cohort starts here at age 20.',
                      'Precancerous lesion, which can regress or progress.',
                      'Undiagnosed localized cancer.', 'Undiagnosed advanced cancer.',
                      'Cancer diagnosed after symptoms.', 'Death from other causes.'),
      level = c(1, 2, 3, 4, 5, 6)
    ),
    edges = data.frame(
      from = c('healthy', 'lesion', 'lesion', 'pre.early', 'pre.early', 'pre.late'),
      to = c('lesion', 'healthy', 'pre.early', 'pre.late', 'diagnosed', 'diagnosed'),
      label = c('Onset', 'Regression', 'Progression', 'Stage progression', 'Symptoms', 'Symptoms'),
      description = c('p.lesion.onset', 'p.lesion.regress', 'p.lesion.progress',
                      'p.stage.progress', 'p.symptoms.early', 'p.symptoms.late')
    )
  )
}

get.code.sample <- function() {
  list(`model.R` = paste(readLines('model.R'), collapse = '\n'))
}

run.simulation <- function(strategies, pars) {
  simulate(strategies, pars)
}

# Cancer incidence comes from the cancer registry, lesion prevalence from a
# population study that examined a sample of every age group, and the stage at
# diagnosis from the registry's staging records. Each value is an estimate with
# its own sampling error.
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
CALIB.PARAMS <- c('p.lesion.onset', 'p.lesion.progress', 'p.stage.progress')

# Relative squared error, averaged within each series and summed over series,
# so that incidence per 100,000 and the proportions weigh alike.
calibration.output.error <- function(output, target) {
  error <- 0
  for (series in names(target)) {
    target.values <- unlist(target[[series]])
    error <- error + mean(((output[[series]][names(target.values)] - target.values) / target.values)^2)
  }
  error
}

calibration.error <- function(pars, target) {
  strategy <- 'no_intervention'
  tryCatch({
    simulated <- run.simulation(strategy, pars)$outputs[[strategy]]
    output <- list()
    for (series in names(target))
      output[[series]] <- simulated[[series]][names(unlist(target[[series]]))]
    error <- calibration.output.error(output, target)
    if (!is.finite(error)) stop('the simulation gave no usable output')
    list(error = error, output = output)
  }, error = function(e) {
    output <- lapply(target, function(series) setNames(rep(NA_real_, length(series)), names(series)))
    list(error = Inf, output = output)
  })
}

# The calibration vector is laid out parameter by parameter, and within a
# parameter stratum by stratum.
initial.guess <- function() {
  base <- default.parameters()
  unname(unlist(lapply(CALIB.PARAMS, function(param) unlist(base[[param]])[MODEL.STRATA])))
}

# The initial guess scaled by a uniform factor per value, within [0, 1].
generate.training.dataset <- function(initial_guess, n, ...) {
  variation <- list(...)$variation
  if (is.null(variation)) variation <- 1
  k <- length(initial_guess)
  factors <- matrix(runif(n * k, min = max(0, 1 - variation), max = 1 + variation), nrow = n)
  pmin(pmax(sweep(factors, 2, initial_guess, `*`), 0), 1)
}

# The mean time a lesion lasts before it regresses or progresses, 1 / (regression
# + progression), lies between 15.5 and 24 years in every stratum, the range
# reported for this lesion by follow-up studies. One value per bound and
# stratum, in years, in the c(x) <= 0 convention.
#
# `partial` says the vector covers only the strata a stepwise calibration has
# got through, in the same layout over those strata. Each constraint is read
# within one stratum, so the vector answers for the strata it holds.
DWELL.MIN <- 15.5
DWELL.MAX <- 24
LESION.REGRESS <- default.parameters()$p.lesion.regress
calibration.constraints <- function(x, partial = FALSE) {
  x <- as.numeric(x)
  n.strata <- length(x) %/% length(CALIB.PARAMS)
  j <- match('p.lesion.progress', CALIB.PARAMS)
  progress <- x[(j - 1) * n.strata + seq_len(n.strata)]
  dwell <- 1 / (LESION.REGRESS + progress)
  as.vector(rbind(DWELL.MIN - dwell, dwell - DWELL.MAX))
}

calibration.constraints.llm <- function() {
  paste0(
    'The vector holds ', paste(CALIB.PARAMS, collapse = ', '), ', in that order, each ',
    'as a block of ', length(MODEL.STRATA), ' values running over the age groups ',
    paste(MODEL.STRATA, collapse = ', '), '.\n',
    'In every age group, the mean lesion dwell time, 1 / (p.lesion.regress + ',
    'p.lesion.progress), lies between ', DWELL.MIN, ' and ', DWELL.MAX, ' years. ',
    'p.lesion.regress is not calibrated and stays at ', LESION.REGRESS, '.'
  )
}

# The training set of the constrained scheme: draws as above, of which only
# those meeting the constraints are kept.
generate.feasible.training.dataset <- function(initial_guess, n, ...) {
  rows <- list()
  kept <- 0
  attempts <- 0
  while (kept < n && attempts < 100) {
    attempts <- attempts + 1
    draws <- generate.training.dataset(initial_guess, 4 * n, ...)
    feasible <- apply(draws, 1, function(x) all(calibration.constraints(x) <= 0))
    rows[[attempts]] <- draws[feasible, , drop = FALSE]
    kept <- kept + sum(feasible)
  }
  dataset <- do.call(rbind, rows)
  dataset[seq_len(min(n, nrow(dataset))), , drop = FALSE]
}

# Simulated against observed lesion prevalence, by age group. The plot is given
# the values of the evaluation that belong to a stratum, named after it.
plot.prevalence.fit <- function(data) {
  if (!requireNamespace('ggplot2', quietly = TRUE)) return(NULL)
  prefix <- 'sim_Lesion prevalence.'
  simulated <- data[startsWith(names(data), prefix)]
  if (length(simulated) == 0) return(NULL)
  groups <- sub(prefix, '', names(simulated), fixed = TRUE)
  observed <- unlist(CALIB.TARGET$`Lesion prevalence`)[groups]
  df <- data.frame(group = factor(rep(groups, 2), levels = groups),
                   value = c(as.numeric(simulated), as.numeric(observed)),
                   source = rep(c('Simulated', 'Target'), each = length(groups)))
  df$text <- sprintf('%s, %s: %.1f%%', df$group, df$source, 100 * df$value)
  ggplot2::ggplot(df, ggplot2::aes(x = group, y = value, colour = source,
                                   group = source, text = text)) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::labs(x = 'Age group', y = 'Lesion prevalence', colour = '') +
    ggplot2::theme_minimal()
}

get.calibration.schemes <- function() {
  list(
    standard = list(
      description = 'Lesion onset, progression and stage progression by age',
      parameters = CALIB.PARAMS,
      target = CALIB.TARGET,
      strata = MODEL.STRATA,
      initial_guess = initial.guess(),
      error_function = calibration.error,
      output_error = calibration.output.error,
      latent_space_training_set = generate.training.dataset,
      evaluation.plots = list(`Lesion prevalence` = plot.prevalence.fit)
    ),
    constrained = list(
      description = 'Lesion onset, progression and stage progression by age, within the lesion dwell time',
      parameters = CALIB.PARAMS,
      target = CALIB.TARGET,
      strata = MODEL.STRATA,
      initial_guess = initial.guess(),
      error_function = calibration.error,
      output_error = calibration.output.error,
      constraints = calibration.constraints,
      constraints_description = calibration.constraints.llm,
      latent_space_training_set = generate.feasible.training.dataset,
      evaluation.plots = list(`Lesion prevalence` = plot.prevalence.fit)
    )
  )
}
