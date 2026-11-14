source('model.R')

get.overview <- function() {
  paste(readLines('overview.md'), collapse='\n')
}

get.strata <- function() {
  MODEL.STRATA
}

get.strategies <- function() {
  list(
    list(name='natural_history', display.name='Natural history',
         description='No screening and no vaccination: cancers are diagnosed when they cause symptoms.')
  )
}

get.model.settings <- function() {
  list(cost.unit='€', effectiveness.unit='LY')
}

# Prior values of the age- and sex-specific transitions, women's six age groups
# then men's, in the order of MODEL.STRATA. The calibration starts from them.
PRIOR.LESION.ONSET <- c(0.028, 0.024, 0.012, 0.007, 0.005, 0.004,
                        0.022, 0.020, 0.012, 0.008, 0.006, 0.004)
PRIOR.LESION.CLEARANCE <- c(0.012, 0.010, 0.016, 0.016, 0.015, 0.014,
                            0.010, 0.009, 0.015, 0.015, 0.013, 0.012)
PRIOR.LESION.PROGRESSION <- rep(c(0.0021, 0.0030, 0.0039, 0.0039, 0.0039, 0.0039), 2)
PRIOR.HIGH.GRADE.REGRESSION <- rep(c(0.55, 0.15, 0.09, 0.05, 0.035, 0.03), 2)
PRIOR.INVASION <- rep(c(0.17, 0.10, 0.14, 0.14, 0.15, 0.16), 2)

stratified.parameter <- function(name, display.name, values, class) {
  lapply(seq_along(MODEL.STRATA), function(i) {
    list(name=name,
         display.name=display.name,
         base.value=values[i],
         stratum=MODEL.STRATA[i],
         class=class,
         distribution='beta',
         min.value=0,
         max.value=1)
  })
}

MODEL.PARAMETERS <- c(
  stratified.parameter('p.lesion.onset',
                       'Annual probability of acquiring a persistent infection with a low-grade lesion',
                       PRIOR.LESION.ONSET, 'Lesions'),
  stratified.parameter('p.lesion.clearance',
                       'Annual probability that a low-grade lesion clears',
                       PRIOR.LESION.CLEARANCE, 'Lesions'),
  stratified.parameter('p.lesion.progression',
                       'Annual probability that a low-grade lesion becomes high-grade',
                       PRIOR.LESION.PROGRESSION, 'Lesions'),
  stratified.parameter('p.high.grade.regression',
                       'Annual probability that a high-grade lesion regresses to low-grade',
                       PRIOR.HIGH.GRADE.REGRESSION, 'Lesions'),
  stratified.parameter('p.invasion',
                       'Annual probability that a high-grade lesion becomes a preclinical cancer',
                       PRIOR.INVASION, 'Cancer'),
  list(
    list(name='p.symptomatic.diagnosis',
         display.name='Annual probability that a preclinical cancer is diagnosed through symptoms',
         base.value=0.30, class='Cancer', distribution='beta', min.value=0, max.value=1),
    list(name='mortality.base.F',
         display.name='Other-cause mortality of women at age 20',
         base.value=0.0004, class='Mortality', distribution='beta', min.value=0, max.value=1),
    list(name='mortality.base.M',
         display.name='Other-cause mortality of men at age 20',
         base.value=0.0007, class='Mortality', distribution='beta', min.value=0, max.value=1),
    list(name='mortality.gompertz.rate',
         display.name='Gompertz rate of other-cause mortality per year of age',
         base.value=0.085, class='Mortality', distribution='lognormal', min.value=0),
    list(name='cost.cancer',
         display.name='Cost of diagnosing and treating a cancer',
         base.value=30000, class='Costs', distribution='cost', min.value=0)
  )
)

get.parameters <- function() {
  MODEL.PARAMETERS
}

# Base values in the shape run.simulation() receives them: stratified
# parameters as named lists over the strata.
default.parameters <- function() {
  pars <- list()
  for (descriptor in MODEL.PARAMETERS) {
    if (is.null(descriptor$stratum)) {
      pars[[descriptor$name]] <- descriptor$base.value
    } else {
      if (is.null(pars[[descriptor$name]])) pars[[descriptor$name]] <- list()
      pars[[descriptor$name]][[descriptor$stratum]] <- descriptor$base.value
    }
  }
  pars
}

run.simulation <- function(strategies, pars, strict=FALSE) {
  results <- simulate(pars, strict=strict)
  summary <- data.frame(strategy=strategies,
                        C=results$cancers * pars$cost.cancer,
                        E=results$life.years)
  list(summary=summary,
       outputs=setNames(rep(list(results$outputs), length(strategies)), strategies),
       traces=setNames(rep(list(results$traces), length(strategies)), strategies))
}

# Registry incidence and a prevalence survey, by sex and age group.
read.targets <- function() {
  rows <- read.csv('data/targets.csv', stringsAsFactors=FALSE)
  series <- unique(rows$series)
  setNames(lapply(series, function(s) {
    block <- rows[rows$series == s, ]
    as.list(setNames(block$value, block$stratum))
  }), series)
}
CALIB.TARGET <- read.targets()

CALIB.PARAMS <- CALIBRATED.TRANSITIONS
CALIB.STRATA <- MODEL.STRATA

# The calibration vector holds the parameters of CALIB.PARAMS in that order,
# each as a block of one value per stratum of CALIB.STRATA.
initial.guess.for <- function(params) {
  base <- default.parameters()
  unname(unlist(lapply(params, function(param) {
    sapply(CALIB.STRATA, function(stratum) base[[param]][[stratum]])
  })))
}

# Relative squared error, averaged within each target series and summed over
# series: incidence and prevalence are on scales orders of magnitude apart.
calibration.output.error <- function(output, target) {
  error <- 0
  for (series in names(target)) {
    target.values <- unlist(target[[series]])
    error <- error + mean(((output[[series]][names(target.values)] - target.values) / target.values)^2)
  }
  error
}

scheme.error <- function(pars, target, strict) {
  strategy <- 'natural_history'
  tryCatch({
    simulated <- run.simulation(strategy, pars, strict=strict)$outputs[[strategy]]
    output <- list()
    for (series in names(target))
      output[[series]] <- simulated[[series]][names(unlist(target[[series]]))]
    list(error=calibration.output.error(output, target), output=output)
  }, error=function(e) {
    output <- lapply(target, function(series) {
      setNames(rep(NA_real_, length(series)), names(series))
    })
    list(error=Inf, output=output)
  })
}

# Exits of a state that add up to more than 1 are scaled down to share it.
calibration.error <- function(pars, target) {
  scheme.error(pars, target, strict=FALSE)
}

# The same model with exits that add up to more than 1 treated as an error.
calibration.error.strict <- function(pars, target) {
  scheme.error(pars, target, strict=TRUE)
}

# Training set for the latent space methods: the initial guess scaled by a
# uniform random factor per value, kept inside [0, 1].
generate.training.dataset <- function(initial_guess, n, ...) {
  variation <- list(...)$variation
  if (is.null(variation)) variation <- 1
  k <- length(initial_guess)
  factors <- matrix(runif(n * k, min=max(0, 1 - variation), max=1 + variation), nrow=n)
  pmin(pmax(sweep(factors, 2, initial_guess, `*`), 0), 1)
}

# The largest change of each transition between neighbouring age groups of a
# sex that is considered credible: twice the largest such change of its prior.
AGE.STEP.LIMIT <- c(p.lesion.onset=0.024, p.lesion.clearance=0.012, p.lesion.progression=0.0018,
                    p.high.grade.regression=0.8, p.invasion=0.14)

# The neighbouring age groups of each sex, as pairs of positions in CALIB.STRATA.
age.step.pairs <- function(n.strata=length(CALIB.STRATA)) {
  pairs <- do.call(rbind, lapply(SEXES, function(sex) {
    at <- which(startsWith(CALIB.STRATA, paste0(sex, '.')))
    cbind(head(at, -1), tail(at, -1))
  }))
  pairs[pairs[, 2] <= n.strata, , drop=FALSE]
}

# Smoothness in age, in the c(x) <= 0 convention: for every transition and sex,
# |p(a + 1) - p(a)| <= limit, returned as |p(a + 1) - p(a)| / limit - 1, one value
# per transition, sex and pair of neighbouring age groups, in that order.
# With `partial` the vector covers the strata calibrated so far, the first of
# CALIB.STRATA, in the same layout, and only the pairs among them are answered.
# The leading constant is satisfied by construction and gives a step with no
# pair yet a constraint to model.
calibration.constraints <- function(x, partial=FALSE) {
  x <- as.numeric(x)
  n.strata <- length(x) %/% length(CALIB.PARAMS)
  pairs <- age.step.pairs(n.strata)
  values <- unlist(lapply(seq_along(CALIB.PARAMS), function(j) {
    block <- x[(j - 1) * n.strata + seq_len(n.strata)]
    abs(block[pairs[, 2]] - block[pairs[, 1]]) / AGE.STEP.LIMIT[[CALIB.PARAMS[j]]] - 1
  }))
  c(if (partial) -1, values)
}

calibration.constraints.llm <- function() {
  paste0(
    'The vector holds ', paste(CALIB.PARAMS, collapse=', '), ', in that order, each as a block of ',
    length(CALIB.STRATA), ' values running over the strata ', paste(CALIB.STRATA, collapse=', '), '.\n',
    'Each transition changes smoothly with age: within each sex, the values of neighbouring age groups ',
    'differ by at most ',
    paste(sprintf('%g for %s', AGE.STEP.LIMIT, names(AGE.STEP.LIMIT)), collapse=', '), '.'
  )
}

# Training set of the constrained scheme, feasible by construction: per
# transition and sex, the first age group is drawn as in
# generate.training.dataset(), and each next one uniformly among the values of
# its range that keep the step from the one before it within the limit.
generate.smooth.training.dataset <- function(initial_guess, n, ...) {
  variation <- list(...)$variation
  if (is.null(variation)) variation <- 1
  lower <- pmax(0, initial_guess * (1 - variation))
  upper <- pmin(1, initial_guess * (1 + variation))
  n.strata <- length(CALIB.STRATA)
  dataset <- matrix(NA_real_, nrow=n, ncol=length(initial_guess))
  for (j in seq_along(CALIB.PARAMS)) {
    # Drawn a hair inside the limit, so that rounding never puts a step past it.
    limit <- AGE.STEP.LIMIT[[CALIB.PARAMS[j]]] * (1 - 1e-9)
    for (sex in SEXES) {
      cols <- (j - 1) * n.strata + which(startsWith(CALIB.STRATA, paste0(sex, '.')))
      previous <- NULL
      for (col in cols) {
        if (is.null(previous)) {
          value <- runif(n, lower[col], upper[col])
        } else {
          lo <- pmax(lower[col], previous - limit)
          hi <- pmin(upper[col], previous + limit)
          value <- runif(n, pmin(lo, hi), pmax(lo, hi))
          # A range left empty by the bounds of the box: the nearest credible step.
          value <- pmin(pmax(value, previous - limit), previous + limit)
        }
        dataset[, col] <- value
        previous <- value
      }
    }
  }
  dataset
}

# Simulated against registry incidence, by sex and age group. The plot is given
# the values of the evaluation that belong to a stratum, named after it.
plot.incidence.fit <- function(data) {
  if (!requireNamespace('ggplot2', quietly=TRUE)) return(NULL)
  prefix <- 'sim_Cancer incidence.'
  simulated <- data[startsWith(names(data), prefix)]
  if (length(simulated) == 0) return(NULL)
  strata <- sub(prefix, '', names(simulated), fixed=TRUE)
  observed <- unlist(CALIB.TARGET$`Cancer incidence`)[strata]
  df <- data.frame(sex=ifelse(startsWith(strata, 'F.'), 'Women', 'Men'),
                   group=factor(sub('^[FM]\\.', '', strata), levels=AGE.GROUPS),
                   value=1e5 * c(as.numeric(simulated), as.numeric(observed)),
                   source=rep(c('Simulated', 'Target'), each=length(strata)))
  df$text <- sprintf('%s %s, %s: %.1f per 100,000', df$sex, df$group, df$source, df$value)
  ggplot2::ggplot(df, ggplot2::aes(x=group, y=value, colour=sex, linetype=source,
                                   group=interaction(sex, source), text=text)) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::labs(x='Age group', y='Cancer incidence per 100,000', colour='', linetype='') +
    ggplot2::theme_minimal()
}

get.calibration.schemes <- function() {
  scheme <- function(description, error_function) list(
    description=description,
    parameters=CALIB.PARAMS,
    target=CALIB.TARGET,
    strata=CALIB.STRATA,
    initial_guess=initial.guess.for(CALIB.PARAMS),
    error_function=error_function,
    output_error=calibration.output.error,
    latent_space_training_set=generate.training.dataset,
    evaluation.plots=list(`Cancer incidence` = plot.incidence.fit)
  )
  constrained <- scheme('Natural history by sex and age, smooth in age', calibration.error)
  constrained$constraints <- calibration.constraints
  constrained$constraints_description <- calibration.constraints.llm
  constrained$latent_space_training_set <- generate.smooth.training.dataset
  list(
    standard=scheme('Natural history by sex and age', calibration.error),
    failures=scheme('Natural history by sex and age, strict transitions', calibration.error.strict),
    constrained=constrained
  )
}

get.model.states <- function() {
  list(
    description='Annual transitions of each cohort. Every living state also leads to other-cause death, which is not drawn.',
    nodes=data.frame(
      id=c('healthy', 'low.grade', 'high.grade', 'preclinical', 'diagnosed'),
      label=c('Healthy', 'Low-grade\nlesion', 'High-grade\nlesion', 'Preclinical\ncancer', 'Diagnosed\ncancer'),
      description=c('No persistent infection and no lesion. Both cohorts start here at age 20.',
                    'Persistent infection with a low-grade lesion.',
                    'High-grade lesion, the precursor of cancer.',
                    'Undiagnosed invasive cancer.',
                    'Cancer diagnosed through symptoms.')
    ),
    edges=data.frame(
      from=c('healthy', 'low.grade', 'low.grade', 'high.grade', 'high.grade', 'preclinical'),
      to=c('low.grade', 'healthy', 'high.grade', 'low.grade', 'preclinical', 'diagnosed'),
      label=c('Onset', 'Clearance', 'Progression', 'Regression', 'Invasion', 'Symptoms'),
      description=c('p.lesion.onset', 'p.lesion.clearance', 'p.lesion.progression',
                    'p.high.grade.regression', 'p.invasion', 'p.symptomatic.diagnosis')
    )
  )
}

get.code.sample <- function() {
  list(`model.R`=paste(readLines('model.R'), collapse='\n'))
}
