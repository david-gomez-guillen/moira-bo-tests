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
  stratified.parameter('p.lesion.onset',
                       'Annual probability that a lesion-free person develops a high-grade lesion',
                       c(0.0021, 0.0040, 0.0050, 0.0060), 'Natural history'),
  stratified.parameter('p.lesion.regress',
                       'Annual probability that a high-grade lesion regresses',
                       c(0.18, 0.16, 0.14, 0.14), 'Natural history'),
  stratified.parameter('p.lesion.progress',
                       'Annual probability that a high-grade lesion becomes a preclinical cancer',
                       c(0.032, 0.040, 0.042, 0.042), 'Natural history'),
  list(
    probability('p.symptoms',
                'Annual probability that a preclinical cancer causes symptoms and is diagnosed',
                0.25, 'Natural history'),
    probability('p.cancer.death', 'Annual probability of death from a diagnosed cancer',
                0.08, 'Cancer survival'),
    probability('mortality.other.base', 'Other-cause mortality at age 20', 0.0005, 'General'),
    list(name = 'mortality.other.rate', display.name = 'Gompertz rate of other-cause mortality',
         base.value = 0.085, class = 'General', distribution = 'lognormal', min.value = 0),
    list(name = 'cost.diagnosis', display.name = 'Cost of diagnosing and treating a cancer',
         base.value = 20000, class = 'Costs', distribution = 'gamma',
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
      id = c('healthy', 'lesion', 'preclinical', 'diagnosed'),
      label = c('Lesion-free', 'High-grade\nlesion', 'Preclinical\ncancer', 'Diagnosed\ncancer'),
      description = c('No lesion. The whole cohort starts here at age 20.',
                      'Precancerous high-grade lesion, which can regress or progress.',
                      'Undiagnosed cancer, without symptoms.',
                      'Cancer diagnosed after symptoms.'),
      level = c(1, 2, 3, 4)
    ),
    edges = data.frame(
      from = c('healthy', 'lesion', 'lesion', 'preclinical'),
      to = c('lesion', 'healthy', 'preclinical', 'diagnosed'),
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

# Incidence comes from the cancer registry, prevalence from a population study
# that examined people of every age group, and persistence from a cohort of
# untreated lesions examined again two years after they were found.
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

CALIB.PARAMS <- c('p.lesion.onset', 'p.lesion.regress', 'p.lesion.progress')

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
    if (!is.finite(error)) stop('the error is not finite')
    list(error = error, output = output)
  }, error = function(e) {
    output <- lapply(target, function(series) {
      setNames(rep(NA_real_, length(series)), names(series))
    })
    list(error = Inf, output = output)
  })
}

# What the clinical literature says about the natural history, whatever the
# targets: a high-grade lesion lasts between 3 and 6 years on average before it
# regresses or progresses, at most a quarter of lesions progress rather than
# regress, and at most 32% of people develop a lesion by age 85.
SOJOURN.MIN <- 3
SOJOURN.MAX <- 6
PROGRESSION.SHARE.MAX <- 0.25
LIFETIME.RISK.MAX <- 0.32

# The constraints of every row of X, a matrix of calibration vectors laid out
# parameter by parameter and, within a parameter, stratum by stratum, over the
# first ncol(X) / 3 strata. One column per constraint, c <= 0 when it is met:
# the mean sojourn bounds and the progression share for each of those strata,
# then, when every stratum is there, the lifetime risk.
constraint.matrix <- function(X, partial = FALSE) {
  if (!is.matrix(X)) X <- matrix(as.numeric(X), nrow = 1)
  k <- ncol(X) %/% length(CALIB.PARAMS)
  block <- function(j) X[, (j - 1) * k + seq_len(k), drop = FALSE]
  onset <- block(1); regress <- block(2); progress <- block(3)
  # A floor on the exit probability keeps the sojourn finite at the corner of
  # the box where both are zero.
  exit <- pmax(regress + progress, 1e-3)
  sojourn <- 1 / exit
  out <- cbind(sojourn - SOJOURN.MAX, SOJOURN.MIN - sojourn,
               progress / exit - PROGRESSION.SHARE.MAX)
  if (!partial && k == length(MODEL.STRATA)) {
    never <- exp(log1p(-pmin(onset, 1 - 1e-12)) %*% YEARS.IN.STRATUM)
    out <- cbind(out, (1 - never) - LIFETIME.RISK.MAX)
  }
  out
}

# For the scheme, on one calibration vector. With partial = TRUE the vector
# covers the strata a stepwise run has reached so far: those decide their own
# sojourn and progression share, and the lifetime risk waits for every stratum.
calibration.constraints <- function(x, partial = FALSE) {
  as.numeric(constraint.matrix(matrix(as.numeric(x), nrow = 1), partial))
}

calibration.constraints.description <- function() {
  paste0(
    'The vector holds ', paste(CALIB.PARAMS, collapse = ', '), ', in that order, ',
    'each as a block of ', length(MODEL.STRATA), ' values over the age groups ',
    paste(MODEL.STRATA, collapse = ', '), '.\n',
    'In every age group, the mean sojourn of a high-grade lesion, ',
    '1 / (p.lesion.regress + p.lesion.progress), is between ', SOJOURN.MIN, ' and ',
    SOJOURN.MAX, ' years.\n',
    'In every age group, the share of lesions that progress rather than regress, ',
    'p.lesion.progress / (p.lesion.regress + p.lesion.progress), is at most ',
    PROGRESSION.SHARE.MAX, '.\n',
    'The lifetime risk of a first lesion between ages 20 and 85 in the absence of death, ',
    '1 - prod((1 - p.lesion.onset)^years), with 15, 15, 15 and 21 years in the four age groups, ',
    'is at most ', LIFETIME.RISK.MAX, '.'
  )
}

# Feasible samples for the latent space: draws from the search box around the
# initial guess, the scope given by `variation`, kept when they meet every
# constraint. Drawn in batches until there are n.
generate.constrained.training.dataset <- function(initial_guess, n, ...) {
  variation <- list(...)$variation
  if (is.null(variation)) variation <- 1
  lower <- pmax(0, initial_guess * (1 - variation))
  upper <- pmin(1, initial_guess * (1 + variation))
  d <- length(initial_guess)
  kept <- matrix(NA_real_, nrow = 0, ncol = d)
  batch <- max(5000, 50 * n)
  for (attempt in 1:200) {
    draws <- matrix(runif(batch * d, rep(lower, each = batch), rep(upper, each = batch)),
                    ncol = d)
    feasible <- rowSums(constraint.matrix(draws) > 0) == 0
    kept <- rbind(kept, draws[feasible, , drop = FALSE])
    if (nrow(kept) >= n) break
  }
  if (nrow(kept) < n) stop('Too few feasible draws for the training set')
  kept[seq_len(n), , drop = FALSE]
}

# Simulated against observed lesion prevalence, by age group. The plot is
# given the values of the evaluation that belong to a stratum, named after it.
plot.prevalence.fit <- function(data) {
  if (!requireNamespace('ggplot2', quietly = TRUE)) return(NULL)
  prefix <- 'sim_Lesion prevalence.'
  simulated <- data[startsWith(names(data), prefix)]
  if (length(simulated) == 0) return(NULL)
  groups <- sub(prefix, '', names(simulated), fixed = TRUE)
  observed <- unlist(CALIB.TARGET$`Lesion prevalence`)[groups]
  df <- data.frame(group = factor(rep(groups, 2), levels = groups),
                   value = 100 * c(as.numeric(simulated), as.numeric(observed)),
                   source = rep(c('Simulated', 'Target'), each = length(groups)))
  df$text <- sprintf('%s, %s: %.2f%%', df$group, df$source, df$value)
  ggplot2::ggplot(df, ggplot2::aes(x = group, y = value, colour = source,
                                   group = source, text = text)) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::labs(x = 'Age group', y = 'Lesion prevalence (%)', colour = '') +
    ggplot2::theme_minimal()
}

initial.guess <- function() {
  base <- default.parameters()
  unname(unlist(lapply(CALIB.PARAMS, function(param) unlist(base[[param]])[MODEL.STRATA])))
}

get.calibration.schemes <- function() {
  list(
    constrained = list(
      description = 'Lesion natural history within clinical bounds',
      parameters = CALIB.PARAMS,
      target = CALIB.TARGET,
      strata = MODEL.STRATA,
      initial_guess = initial.guess(),
      error_function = calibration.error,
      output_error = calibration.output.error,
      constraints = calibration.constraints,
      constraints_description = calibration.constraints.description,
      latent_space_training_set = generate.constrained.training.dataset,
      evaluation.plots = list(`Lesion prevalence` = plot.prevalence.fit)
    )
  )
}
