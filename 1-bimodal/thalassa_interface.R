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

# The base values of the calibrated parameters are the prior of the modeller:
# a fast pathway that carries a large part of the cancers.
BASE.LGL.ONSET <- c(0.0039, 0.0063, 0.0080, 0.0098)
BASE.FAST.ONSET <- c(0.0060, 0.0160, 0.0210, 0.0204)
BASE.HGL.PROGRESS <- c(0.065, 0.080, 0.098, 0.107)

MODEL.PARAMETERS <- c(
  stratified.parameter('p.lgl.onset',
                       'Annual probability that a healthy person develops a low-grade lesion',
                       BASE.LGL.ONSET, 'Natural history (slow pathway)'),
  stratified.parameter('p.fast.onset',
                       'Annual probability that a susceptible healthy person develops a fast-pathway lesion',
                       BASE.FAST.ONSET, 'Natural history (fast pathway)'),
  stratified.parameter('p.hgl.progress',
                       'Annual probability that a high-grade lesion becomes a preclinical cancer',
                       BASE.HGL.PROGRESS, 'Natural history (slow pathway)'),
  list(
    probability('p.lgl.progress', 'Annual probability that a low-grade lesion becomes high-grade',
                0.06, 'Natural history (slow pathway)'),
    probability('p.lgl.regress', 'Annual probability that a low-grade lesion regresses',
                0.06, 'Natural history (slow pathway)'),
    probability('p.hgl.regress', 'Annual probability that a high-grade lesion regresses to low-grade',
                0.063, 'Natural history (slow pathway)'),
    probability('prop.susceptible', 'Proportion of the cohort susceptible to the fast pathway',
                0.39, 'Natural history (fast pathway)'),
    probability('p.fast.progress', 'Annual probability that a fast-pathway lesion becomes a preclinical cancer',
                0.39, 'Natural history (fast pathway)'),
    probability('p.fast.clear', 'Annual probability that a fast-pathway lesion clears',
                0.41, 'Natural history (fast pathway)'),
    probability('p.symptoms', 'Annual probability that a preclinical cancer causes symptoms and is diagnosed',
                0.31, 'Cancer'),
    probability('p.cancer.death', 'Annual probability of cancer death after diagnosis',
                0.08, 'Cancer'),
    probability('mortality.other.base', 'Other-cause mortality at age 20', 0.0005, 'General'),
    list(name = 'mortality.other.rate', display.name = 'Gompertz rate of other-cause mortality',
         base.value = 0.085, class = 'General', distribution = 'lognormal', min.value = 0),
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
      id = c('healthy.susceptible', 'healthy', 'low.grade', 'high.grade', 'fast.lesion',
             'preclinical', 'clinical'),
      label = c('Healthy\n(susceptible)', 'Healthy', 'Low-grade\nlesion', 'High-grade\nlesion',
                'Fast-pathway\nlesion', 'Preclinical\ncancer', 'Clinical\ncancer'),
      description = c('No lesion, and at risk of both pathways.',
                      'No lesion, and at risk of the slow pathway only.',
                      'Low-grade lesion, the first step of the slow pathway.',
                      'High-grade lesion, the precursor of slow-pathway cancer.',
                      'Short-lived lesion of the fast pathway.',
                      'Undiagnosed cancer, without symptoms.',
                      'Cancer diagnosed after symptoms.'),
      level = c(1, 1, 2, 3, 3, 4, 5)
    ),
    edges = data.frame(
      from = c('healthy.susceptible', 'healthy.susceptible', 'healthy', 'low.grade', 'low.grade',
               'high.grade', 'high.grade', 'fast.lesion', 'fast.lesion', 'preclinical'),
      to = c('low.grade', 'fast.lesion', 'low.grade', 'high.grade', 'healthy',
             'preclinical', 'low.grade', 'preclinical', 'healthy', 'clinical'),
      label = c('Onset', 'Onset', 'Onset', 'Progression', 'Regression',
                'Progression', 'Regression', 'Progression', 'Clearance', 'Symptoms'),
      description = c('p.lgl.onset', 'p.fast.onset', 'p.lgl.onset', 'p.lgl.progress', 'p.lgl.regress',
                      'p.hgl.progress', 'p.hgl.regress', 'p.fast.progress', 'p.fast.clear', 'p.symptoms')
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
  simulate(strategies, pars)
}

# Clinical incidence comes from the cancer registry. Precursor prevalence comes
# from a screening study that examined a sample of each age group for lesions
# of either kind.
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

CALIB.PARAMS <- c('p.lgl.onset', 'p.fast.onset', 'p.hgl.progress')

# Relative squared error, averaged within each series and summed over series,
# so that incidence per 100,000 and a proportion weigh alike.
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

# Simulated against target for both series. The plot is given the values of
# the evaluation that belong to a stratum, which are the calibrated parameters
# (the age groups of the outputs are not strata), so it runs the model on them
# with every other parameter at its base value.
plot.fit <- function(data) {
  if (!requireNamespace('ggplot2', quietly = TRUE)) return(NULL)
  pars <- default.parameters()
  for (param in CALIB.PARAMS) for (stratum in MODEL.STRATA) {
    value <- data[paste0(param, '.', stratum)]
    if (length(value) == 1 && !is.na(value)) pars[[param]][[stratum]] <- as.numeric(value)
  }
  simulated <- tryCatch(run.simulation('no_intervention', pars)$outputs$no_intervention,
                        error = function(e) NULL)
  if (is.null(simulated)) return(NULL)
  df <- do.call(rbind, lapply(names(CALIB.TARGET), function(series) {
    observed <- unlist(CALIB.TARGET[[series]])
    groups <- names(observed)
    data.frame(series = series, group = rep(groups, 2),
               value = c(as.numeric(simulated[[series]][groups]), as.numeric(observed)),
               source = rep(c('Simulated', 'Target'), each = length(groups)))
  }))
  df$group <- factor(df$group, levels = unique(df$group))
  df$text <- sprintf('%s, %s: %s', df$group, df$source, signif(df$value, 4))
  ggplot2::ggplot(df, ggplot2::aes(x = group, y = value, colour = source,
                                   group = source, text = text)) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::facet_wrap(~series, scales = 'free') +
    ggplot2::labs(x = 'Age group', y = NULL, colour = '') +
    ggplot2::theme_minimal()
}

# The share of the cancers of a stratum that come from the fast pathway, from
# its three calibrated values alone: each pathway's onset times the chance that
# its lesion goes on to become a preclinical cancer. Pure arithmetic, so the
# constraint below never runs the model.
FAST.SHARE.FIXED <- default.parameters()
stratum.fast.share <- function(lgl.onset, fast.onset, hgl.progress, fixed = FAST.SHARE.FIXED) {
  to.high <- fixed$p.lgl.progress / (fixed$p.lgl.progress + fixed$p.lgl.regress)
  to.cancer <- hgl.progress / (hgl.progress + fixed$p.hgl.regress)
  # A high-grade lesion that regresses becomes low-grade and can progress again.
  slow <- lgl.onset * to.high * to.cancer / (1 - to.high * (1 - to.cancer))
  fast <- fixed$prop.susceptible * fast.onset *
    fixed$p.fast.progress / (fixed$p.fast.progress + fixed$p.fast.clear)
  fast / (slow + fast)
}

# The constrained scheme: a mix of the two pathways near half and half at 70-84
# is only credible if the fast pathway already makes a good part of the cancers
# at 55-69. One value, in the c(x) <= 0 convention: positive only when the fast
# share at 70-84 is between 35% and 65% and the one at 55-69 is below 40%.
#
# The vector is laid out parameter by parameter, each over the strata. With
# `partial`, it covers only the strata calibrated so far, and the constraint,
# which needs the last two, answers as satisfied until they are there.
calibration.constraints <- function(x, partial = FALSE) {
  x <- as.numeric(x)
  n <- length(x) %/% length(CALIB.PARAMS)
  if (n < length(MODEL.STRATA)) return(-1)
  share <- function(k) stratum.fast.share(x[k], x[n + k], x[2 * n + k])
  late <- share(4)
  previous <- share(3)
  min((late - 0.35) * (0.65 - late) / 0.15^2, (0.40 - previous) / 0.40)
}

calibration.constraints.description <- function() {
  paste0(
    'The vector holds ', paste(CALIB.PARAMS, collapse = ', '), ', in that order, each as a ',
    'block of ', length(MODEL.STRATA), ' values over the age groups ',
    paste(MODEL.STRATA, collapse = ', '), '.\n',
    'The fast-pathway share of the cancers of an age group is fast / (slow + fast), where ',
    'slow = p.lgl.onset * a * b / (1 - a * (1 - b)), a = p.lgl.progress / (p.lgl.progress + ',
    'p.lgl.regress), b = p.hgl.progress / (p.hgl.progress + p.hgl.regress), and fast = ',
    'prop.susceptible * p.fast.onset * p.fast.progress / (p.fast.progress + p.fast.clear).\n',
    'If the fast-pathway share at 70-84 is between 35% and 65%, the share at 55-69 must be ',
    'at least 40%.'
  )
}

# Feasible draws only: the uniform draws of generate.training.dataset(), with
# those that break the constraint drawn again.
generate.constrained.training.dataset <- function(initial_guess, n, ...) {
  rows <- matrix(numeric(0), 0, length(initial_guess))
  attempts <- 0
  while (nrow(rows) < n && attempts < 100) {
    attempts <- attempts + 1
    draw <- generate.training.dataset(initial_guess, n, ...)
    rows <- rbind(rows, draw[apply(draw, 1, calibration.constraints) <= 0, , drop = FALSE])
  }
  rows[seq_len(min(n, nrow(rows))), , drop = FALSE]
}

initial.guess <- function() {
  base <- default.parameters()
  unname(unlist(lapply(CALIB.PARAMS, function(param) unlist(base[[param]])[MODEL.STRATA])))
}

get.calibration.schemes <- function() {
  list(
    standard = list(
      description = 'Two-pathway natural history by age group',
      parameters = CALIB.PARAMS,
      target = CALIB.TARGET,
      strata = MODEL.STRATA,
      initial_guess = initial.guess(),
      error_function = calibration.error,
      output_error = calibration.output.error,
      latent_space_training_set = generate.training.dataset,
      evaluation.plots = list(`Fit to the targets` = plot.fit)
    ),
    constrained = list(
      description = 'Two-pathway natural history, consistent pathway shares at 55-84',
      parameters = CALIB.PARAMS,
      target = CALIB.TARGET,
      strata = MODEL.STRATA,
      initial_guess = initial.guess(),
      error_function = calibration.error,
      output_error = calibration.output.error,
      constraints = calibration.constraints,
      constraints_description = calibration.constraints.description,
      latent_space_training_set = generate.constrained.training.dataset,
      evaluation.plots = list(`Fit to the targets` = plot.fit)
    )
  )
}
