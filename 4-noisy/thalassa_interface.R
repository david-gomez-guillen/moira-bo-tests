source('model.R')

get.overview <- function() {
  paste(readLines('overview.md'), collapse = '\n')
}

get.strategies <- function() {
  list(
    list(name = 'no_intervention', display.name = 'No intervention',
         description = 'Natural history only: cancers are diagnosed when they cause symptoms.')
  )
}

probability <- function(name, display.name, base.value, class) {
  list(name = name, display.name = display.name, base.value = base.value,
       class = class, distribution = 'beta', min.value = 0, max.value = 1)
}

# One entry per stratum, as the app expects from get.parameters().
stratified.parameter <- function(name, display.name, values, class) {
  lapply(seq_along(MODEL.STRATA), function(i) {
    list(name = name, display.name = display.name, base.value = values[i],
         stratum = MODEL.STRATA[i], class = class, distribution = 'beta',
         min.value = 0, max.value = 1)
  })
}

MODEL.PARAMETERS <- c(
  stratified.parameter('p.lesion.onset', 'Annual probability of lesion onset',
                       c(0.016, 0.025, 0.030, 0.032), 'Natural history'),
  stratified.parameter('p.lesion.regress', 'Annual probability that a lesion regresses',
                       c(0.120, 0.088, 0.064, 0.052), 'Natural history'),
  stratified.parameter('p.lesion.progress',
                       'Annual probability that a lesion becomes a preclinical cancer',
                       c(0.0046, 0.0050, 0.0060, 0.0072), 'Natural history'),
  stratified.parameter('p.symptoms',
                       'Annual probability that a preclinical cancer causes symptoms and is diagnosed',
                       c(0.15, 0.17, 0.18, 0.21), 'Natural history'),
  list(
    list(name = 'cost.diagnosis', display.name = 'Cost of diagnosing and treating a clinical cancer',
         base.value = 25000, class = 'Costs', distribution = 'gamma',
         min.value = 0, max.value = Inf),
    probability('utility.after.diagnosis', 'Utility of a year lived after a cancer diagnosis',
                0.75, 'Utilities'),
    list(name = 'discount', display.name = 'Discount rate for costs and effects',
         base.value = 0.03, class = 'General', min.value = 0, max.value = 1)
  )
)

get.parameters <- function() MODEL.PARAMETERS

get.strata <- function() MODEL.STRATA

get.model.settings <- function() list(cost.unit = '€', effectiveness.unit = 'QALY')

get.model.states <- function() {
  list(
    nodes = data.frame(
      id = c('healthy', 'lesion', 'preclinical', 'clinical', 'dead'),
      label = c('Healthy', 'Lesion', 'Preclinical\ncancer', 'Clinical\ncancer', 'Other-cause\ndeath'),
      description = c('No lesion.', 'Precancerous lesion, which can regress or progress.',
                      'Undiagnosed cancer, without symptoms.',
                      'Cancer diagnosed after symptoms.', 'Death from other causes.'),
      level = c(1, 2, 3, 4, 5)
    ),
    edges = data.frame(
      from = c('healthy', 'lesion', 'lesion', 'preclinical'),
      to = c('lesion', 'healthy', 'preclinical', 'clinical'),
      label = c('Onset', 'Regression', 'Progression', 'Symptoms'),
      description = c('p.lesion.onset', 'p.lesion.regress', 'p.lesion.progress', 'p.symptoms')
    )
  )
}

# Base values in the shape run.simulation() receives: stratified parameters
# as lists named by stratum.
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

run.simulation <- function(strategies, pars) {
  for (name in CALIB.PARAMS) {
    value <- as.numeric(unlist(pars[[name]])[MODEL.STRATA])
    if (anyNA(value) || any(value < 0 | value > 1))
      stop(sprintf('%s must be a probability in every age group', name))
  }
  simulate(strategies, pars)
}

# Clinical incidence comes from the cancer registry, both prevalences from a
# population study that examined people of every age group for lesions and
# undiagnosed cancers, and lesion persistence from a follow-up of the people
# that study found with a lesion.
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

CALIB.PARAMS <- c('p.lesion.onset', 'p.lesion.regress', 'p.lesion.progress', 'p.symptoms')

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

# Simulated against registry incidence, by age group. The plot is given the
# values of the evaluation that belong to a stratum, named after it.
plot.incidence.fit <- function(data) {
  if (!requireNamespace('ggplot2', quietly = TRUE)) return(NULL)
  prefix <- 'sim_Clinical incidence.'
  simulated <- data[startsWith(names(data), prefix)]
  if (length(simulated) == 0) return(NULL)
  groups <- sub(prefix, '', names(simulated), fixed = TRUE)
  observed <- unlist(CALIB.TARGET$`Clinical incidence`)[groups]
  df <- data.frame(group = factor(rep(groups, 2), levels = groups),
                   value = c(as.numeric(simulated), as.numeric(observed)),
                   source = rep(c('Simulated', 'Target'), each = length(groups)))
  df$text <- sprintf('%s, %s: %.1f per 100,000', df$group, df$source, df$value)
  ggplot2::ggplot(df, ggplot2::aes(x = group, y = value, colour = source,
                                   group = source, text = text)) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::labs(x = 'Age group', y = 'Clinical incidence per 100,000', colour = '') +
    ggplot2::theme_minimal()
}

initial.guess <- function() {
  base <- default.parameters()
  unname(unlist(lapply(CALIB.PARAMS, function(param) unlist(base[[param]])[MODEL.STRATA])))
}

# A lesion is left, by regression or progression, with an annual probability of
# at most LESION.EXIT.MAX in every age group. The vector holds the strata
# calibrated so far when `partial` is TRUE, in the same layout, and each
# constraint lies within one stratum, so it is answered for those strata.
LESION.EXIT.MAX <- 0.18

calibration.constraints <- function(x, partial = FALSE) {
  x <- as.numeric(x)
  n.strata <- length(x) %/% length(CALIB.PARAMS)
  block <- function(param) {
    j <- match(param, CALIB.PARAMS)
    x[(j - 1) * n.strata + seq_len(n.strata)]
  }
  block('p.lesion.regress') + block('p.lesion.progress') - LESION.EXIT.MAX
}

calibration.constraints.llm <- function() {
  n <- length(MODEL.STRATA)
  paste0(
    'The vector holds ', paste(CALIB.PARAMS, collapse = ', '), ', in that order, ',
    'each as a block of ', n, ' values running over the age groups ',
    paste(MODEL.STRATA, collapse = ', '), '.\n',
    'In every age group, p.lesion.regress + p.lesion.progress is at most ',
    LESION.EXIT.MAX, ': the annual probability of leaving the lesion state.'
  )
}

# Same, but only points that meet the constraints, drawn by rejection.
generate.constrained.training.dataset <- function(initial_guess, n, ...) {
  rows <- matrix(numeric(0), nrow = 0, ncol = length(initial_guess))
  while (nrow(rows) < n) {
    draws <- generate.training.dataset(initial_guess, 2 * n, ...)
    feasible <- apply(draws, 1, function(x) all(calibration.constraints(x) <= 0))
    rows <- rbind(rows, draws[feasible, , drop = FALSE])
  }
  rows[seq_len(n), , drop = FALSE]
}

get.calibration.schemes <- function() {
  standard <- list(
    description = 'Natural history by age group',
    parameters = CALIB.PARAMS,
    target = CALIB.TARGET,
    strata = MODEL.STRATA,
    initial_guess = initial.guess(),
    error_function = calibration.error,
    output_error = calibration.output.error,
    latent_space_training_set = generate.training.dataset,
    evaluation.plots = list(`Clinical incidence` = plot.incidence.fit)
  )
  constrained <- standard
  constrained$description <- 'Natural history by age group, bounded lesion exit'
  constrained$constraints <- calibration.constraints
  constrained$constraints_description <- calibration.constraints.llm
  constrained$latent_space_training_set <- generate.constrained.training.dataset
  list(standard = standard, constrained = constrained)
}
