source('model.R')

get.overview <- function() {
  paste(readLines('overview.md'), collapse = '\n')
}

get.strategies <- function() {
  list(list(name = 'no_intervention', display.name = 'No intervention',
            description = 'Natural history only: cancers are diagnosed when they cause symptoms.'))
}

get.strata <- function() {
  MODEL.STRATA
}

# One entry per parameter, the stratified ones with a value per age group.
parameter <- function(name, display.name, value, class, distribution = 'probability',
                      max.value = 1) {
  if (length(value) > 1) value <- as.list(setNames(value, MODEL.STRATA))
  list(name = name, display.name = display.name, base.value = value, class = class,
       distribution = distribution, min.value = 0, max.value = max.value)
}

MODEL.PARAMETERS <- list(
  parameter('p.lgl.onset', 'Annual probability of developing a low-grade lesion',
            c(0.018, 0.025, 0.030, 0.035), 'Natural history (slow pathway)'),
  parameter('p.lgl.progress', 'Annual probability that a low-grade lesion becomes high-grade',
            c(0.045, 0.060, 0.075, 0.090), 'Natural history (slow pathway)'),
  parameter('p.hgl.progress', 'Annual probability that a high-grade lesion becomes a preclinical cancer',
            c(0.006, 0.009, 0.012, 0.015), 'Natural history (slow pathway)'),
  parameter('p.fast.onset', 'Annual probability of developing a preclinical cancer directly, without a detectable lesion',
            c(0.00005, 0.0005, 0.0010, 0.0018), 'Natural history (fast pathway)'),
  parameter('p.lgl.regress', 'Annual probability that a low-grade lesion regresses',
            0.06, 'Natural history (slow pathway)'),
  parameter('p.hgl.regress', 'Annual probability that a high-grade lesion regresses to low grade',
            0.03, 'Natural history (slow pathway)'),
  parameter('p.symptoms.slow', 'Annual probability that a preclinical slow-pathway cancer causes symptoms and is diagnosed',
            0.30, 'Natural history (slow pathway)'),
  parameter('p.symptoms.fast', 'Annual probability that a preclinical fast-pathway cancer causes symptoms and is diagnosed',
            0.55, 'Natural history (fast pathway)'),
  parameter('p.cancer.death', 'Annual probability of cancer death after diagnosis',
            0.06, 'Cancer survival'),
  parameter('p.cure', 'Annual probability of cure after diagnosis',
            0.11, 'Cancer survival'),
  parameter('prevalence.lgl.start', 'Prevalence of low-grade lesions at age 30',
            0.05, 'Cohort at entry'),
  parameter('prevalence.hgl.start', 'Prevalence of high-grade lesions at age 30',
            0.004, 'Cohort at entry'),
  parameter('mortality.other.base', 'Annual other-cause mortality at age 30',
            0.0008, 'General'),
  parameter('mortality.other.rate', 'Gompertz rate of other-cause mortality',
            0.088, 'General', distribution = 'lognormal', max.value = Inf),
  parameter('cost.treatment', 'Annual cost of treating a diagnosed cancer',
            25000, 'Costs', distribution = 'cost', max.value = Inf),
  parameter('utility.cancer', 'Utility of a year with a diagnosed cancer under treatment',
            0.70, 'Utilities', distribution = 'utility'),
  list(name = 'discount', display.name = 'Discount rate for costs and effects',
       base.value = 0.03, class = 'General', min.value = 0, max.value = 1)
)

get.parameters <- function() {
  MODEL.PARAMETERS
}

# Base values as simulate() takes them: stratified ones as lists by age group.
default.parameters <- function() {
  pars <- list()
  for (p in MODEL.PARAMETERS) pars[[p$name]] <- p$base.value
  pars
}

get.model.settings <- function() {
  list(cost.unit = '€', effectiveness.unit = 'QALY')
}

run.simulation <- function(strategies, pars) {
  simulate(strategies, pars)
}

# Registry incidence and cancer mortality, and the lesion prevalence of a
# survey examining people of every age group, all by age group.
read.targets <- function(file = 'data/targets.csv') {
  rows <- read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
  target <- list()
  for (series in unique(rows$series)) {
    part <- rows[rows$series == series, ]
    target[[series]] <- as.list(setNames(part$value, part$label))
  }
  target
}

CALIB.TARGET <- read.targets()
CALIB.PARAMS <- c('p.lgl.onset', 'p.lgl.progress', 'p.hgl.progress', 'p.fast.onset')

calibration.error <- function(pars, target) {
  strategy <- 'no_intervention'
  tryCatch({
    simulated <- run.simulation(strategy, pars)$outputs[[strategy]]
    output <- list()
    for (series in names(target))
      output[[series]] <- simulated[[series]][names(unlist(target[[series]]))]
    error <- calibration.output.error(output, target)
    if (!is.finite(error)) stop('the simulation gave no finite error')
    list(error = error, output = output)
  }, error = function(e) {
    output <- lapply(target, function(series) {
      setNames(rep(NA_real_, length(series)), names(series))
    })
    list(error = Inf, output = output)
  })
}

# Relative squared error, averaged within each series and summed over series.
# Mortality is two orders of magnitude below incidence and four below
# prevalence, so only a scale-free error lets every series count.
calibration.output.error <- function(output, target) {
  error <- 0
  for (series in names(target)) {
    observed <- unlist(target[[series]])
    error <- error + mean(((output[[series]][names(observed)] - observed) / observed)^2)
  }
  error
}

# The initial guess scaled by a uniform random factor per value, within [0, 1].
generate.training.dataset <- function(initial_guess, n, ...) {
  variation <- list(...)$variation
  if (is.null(variation)) variation <- 1
  k <- length(initial_guess)
  factors <- matrix(runif(n * k, min = max(0, 1 - variation), max = 1 + variation), nrow = n)
  pmin(pmax(sweep(factors, 2, initial_guess, `*`), 0), 1)
}

# What is known of the fast pathway apart from the targets: in every age group
# it gives rise to between 10% and 70% of new cancers.
FAST.SHARE.MIN <- 0.10
FAST.SHARE.MAX <- 0.70

# The share of new preclinical cancers that come from the fast pathway, in each
# age group, for each row of X, a calibration vector over the first n.strata
# age groups. It follows the healthy and lesion states of the cohort from age
# 30 with the calibrated probabilities and the fixed regressions. Other-cause
# death strikes every state alike, so it leaves the share unchanged and is left
# out. Each age group depends on those before it only.
fast.cancer.share <- function(X, n.strata) {
  X <- matrix(X, ncol = length(CALIB.PARAMS) * n.strata)
  base <- default.parameters()
  block <- function(j) X[, (j - 1) * n.strata + seq_len(n.strata), drop = FALSE]
  onset <- block(1); progress <- block(2); malignant <- block(3); fast <- block(4)
  healthy <- rep(1 - base$prevalence.lgl.start - base$prevalence.hgl.start, nrow(X))
  low <- rep(base$prevalence.lgl.start, nrow(X))
  high <- rep(base$prevalence.hgl.start, nrow(X))
  from.fast <- from.slow <- matrix(0, nrow(X), n.strata)
  groups <- match(AGE.STRATUM, MODEL.STRATA)
  for (i in which(groups <= n.strata)) {
    g <- groups[i]
    from.fast[, g] <- from.fast[, g] + healthy * fast[, g]
    from.slow[, g] <- from.slow[, g] + high * malignant[, g]
    next.healthy <- healthy * (1 - onset[, g] - fast[, g]) + low * base$p.lgl.regress
    next.low <- low * (1 - progress[, g] - base$p.lgl.regress) + healthy * onset[, g] +
      high * base$p.hgl.regress
    high <- high * (1 - malignant[, g] - base$p.hgl.regress) + low * progress[, g]
    healthy <- next.healthy
    low <- next.low
  }
  from.fast / (from.fast + from.slow)
}

# In the c(x) <= 0 convention, two per age group, in the order of the age
# groups: the share above its minimum, then below its maximum. With partial,
# x covers the age groups calibrated so far, which decide their own shares.
calibration.constraints <- function(x, partial = FALSE) {
  x <- as.numeric(x)
  n.strata <- if (partial) length(x) %/% length(CALIB.PARAMS) else length(MODEL.STRATA)
  share <- fast.cancer.share(x, n.strata)[1, ]
  share[!is.finite(share)] <- 0
  as.numeric(rbind(FAST.SHARE.MIN - share, share - FAST.SHARE.MAX))
}

calibration.constraints.llm <- function() {
  paste0(
    'The vector holds ', paste(CALIB.PARAMS, collapse = ', '), ', in that order, each as ',
    'a block of ', length(MODEL.STRATA), ' values running over the age groups ',
    paste(MODEL.STRATA, collapse = ', '), '.\n',
    'In every age group, the fast pathway gives rise to between ',
    100 * FAST.SHARE.MIN, '% and ', 100 * FAST.SHARE.MAX, '% of new preclinical cancers. ',
    'The share is that of the cohort followed from age 30 through the healthy and lesion ',
    'states: new fast-pathway cancers are those healthy times p.fast.onset, new ',
    'slow-pathway ones those with a high-grade lesion times p.hgl.progress, both summed ',
    'over the ages of the group. It depends on the age group and those before it.'
  )
}

# The training sampler of the constrained scheme: draws of the unconstrained
# one, keeping the feasible ones.
generate.constrained.training.dataset <- function(initial_guess, n, ...) {
  kept <- matrix(numeric(0), 0, length(initial_guess))
  for (attempt in 1:50) {
    draws <- generate.training.dataset(initial_guess, 4 * n, ...)
    share <- fast.cancer.share(draws, length(MODEL.STRATA))
    feasible <- rowSums(!is.finite(share) | share < FAST.SHARE.MIN | share > FAST.SHARE.MAX) == 0
    kept <- rbind(kept, draws[feasible, , drop = FALSE])
    if (nrow(kept) >= n) break
  }
  kept[seq_len(min(n, nrow(kept))), , drop = FALSE]
}

# Simulated against registry incidence, by age group.
plot.incidence.fit <- function(data) {
  if (!requireNamespace('ggplot2', quietly = TRUE)) return(NULL)
  prefix <- 'sim_Clinical incidence.'
  simulated <- data[startsWith(names(data), prefix)]
  if (length(simulated) == 0) return(NULL)
  groups <- sub(prefix, '', names(simulated), fixed = TRUE)
  observed <- unlist(CALIB.TARGET$`Clinical incidence`)[groups]
  df <- data.frame(group = factor(rep(groups, 2), levels = groups),
                   value = 1e5 * c(as.numeric(simulated), as.numeric(observed)),
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

get.calibration.schemes <- function() {
  list(
    standard = list(
      description = 'Natural history of both pathways by age group',
      parameters = CALIB.PARAMS,
      target = CALIB.TARGET,
      strata = MODEL.STRATA,
      initial_guess = initial.guess(),
      error_function = calibration.error,
      output_error = calibration.output.error,
      latent_space_training_set = generate.training.dataset,
      evaluation.plots = list(`Clinical incidence` = plot.incidence.fit)
    ),
    constrained = list(
      description = 'Natural history with a bounded fast-pathway share',
      parameters = CALIB.PARAMS,
      target = CALIB.TARGET,
      strata = MODEL.STRATA,
      initial_guess = initial.guess(),
      error_function = calibration.error,
      output_error = calibration.output.error,
      constraints = calibration.constraints,
      constraints_description = calibration.constraints.llm,
      latent_space_training_set = generate.constrained.training.dataset,
      evaluation.plots = list(`Clinical incidence` = plot.incidence.fit)
    )
  )
}
