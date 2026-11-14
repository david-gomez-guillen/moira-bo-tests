# A calibration scheme of a model written for THALASSA, as the problem
# moira_bo() runs: the folder's thalassa_interface.R read, and the scheme
# turned into the problem THALASSA's Calibration tab builds from it.
#
# The calibration vector holds the parameters of the scheme, each over the
# strata of the model in their order, as the scheme's constraints and sampler
# read it, and starts at the base values of the model. The box is the initial
# guess times 1 +/- scope, within [lower, upper]. Every call into the model is
# made from its folder, which its code reads its files relative to.

`%||%` <- function(a, b) if (is.null(a)) b else a

thalassa_problem <- function(dir, scheme = "constrained", scope = 1, lower = 0, upper = 1) {
  dir <- normalizePath(dir, mustWork = TRUE)
  in_dir <- function(expr) {
    old <- setwd(dir)
    on.exit(setwd(old))
    expr
  }

  # The model in an environment of its own, the files it sources included, so
  # that its functions neither mask nor are masked by those of this program.
  model <- new.env(parent = globalenv())
  model$source <- function(file, local = FALSE, chdir = FALSE, ...)
    invisible(sys.source(file, envir = model, chdir = chdir))
  in_dir(sys.source("thalassa_interface.R", envir = model))

  schemes <- in_dir(model$get.calibration.schemes())
  if (!scheme %in% names(schemes))
    stop(sprintf("the model has no scheme '%s'; its schemes: %s", scheme,
                 paste(names(schemes), collapse = ", ")), call. = FALSE)
  sch <- schemes[[scheme]]
  for (needed in c("error_function", "output_error", "latent_space_training_set"))
    if (!is.function(sch[[needed]]))
      stop(sprintf("the scheme '%s' gives no %s, which MOIRA-BO needs", scheme, needed),
           call. = FALSE)

  strata <- if (exists("get.strata", envir = model, inherits = FALSE))
    in_dir(model$get.strata()) else NULL
  values <- base_values(in_dir(model$get.parameters()))
  slots <- if (length(strata)) strata[strata %in% (sch$strata %||% strata)] else NULL

  # The vector, parameter by parameter and stratum by stratum, and back.
  x0 <- unlist(lapply(sch$parameters, function(p) {
    v <- values[[p]]
    if (is.null(slots)) return(stats::setNames(as.numeric(v), p))
    v <- if (is.list(v)) unlist(v[slots]) else rep(v, length(slots))
    stats::setNames(as.numeric(v), paste0(p, ".", slots))
  }))
  to_parameters <- function(x) {
    pars <- values
    i <- 0L
    for (p in sch$parameters) {
      if (is.null(slots)) {
        pars[[p]] <- x[[i <- i + 1L]]
        next
      }
      if (!is.list(pars[[p]])) pars[[p]] <- stats::setNames(as.list(rep(pars[[p]], length(strata))), strata)
      for (st in slots) pars[[p]][[st]] <- x[[i <- i + 1L]]
    }
    pars
  }

  near <- x0 * max(0, 1 - scope)
  far <- x0 * (1 + scope)
  target <- sch$target

  # Rows of flattened outputs back into the layout the error function returns,
  # learned from the first output it gives.
  layout <- NULL
  rebuild <- function(row) utils::relist(as.numeric(row), layout)

  constraints <- sch[["constraints", exact = TRUE]]
  sampler <- sch$latent_space_training_set

  list(
    name = basename(dir),
    scheme = scheme,
    x0 = x0,
    lower = pmax(lower, pmin(near, far)),
    upper = pmin(upper, pmax(near, far)),
    target = unlist(target),
    evaluate = function(x) {
      res <- in_dir(sch$error_function(to_parameters(x), target))
      if (is.null(layout) && !is.null(res$output)) layout <<- res$output
      res
    },
    score = function(M) {
      M <- as.matrix(M)
      in_dir(vapply(seq_len(nrow(M)), function(i) {
        v <- tryCatch(as.numeric(sch$output_error(rebuild(M[i, ]), target)),
                      error = function(e) NA_real_)
        if (length(v) == 1) v else NA_real_
      }, numeric(1)))
    },
    feasible = function(X) {
      X <- matrix(as.numeric(X), ncol = length(x0))
      if (!is.function(constraints)) return(rep(TRUE, nrow(X)))
      in_dir(apply(X, 1, function(x) {
        v <- tryCatch(as.numeric(constraints(x)), error = function(e) NA_real_)
        length(v) == 0 || all(is.finite(v) & v <= 0)
      }))
    },
    feasibility_sampler = function(center, n, width)
      matrix(as.numeric(in_dir(sampler(as.numeric(center), n, variation = width))), ncol = length(x0)),
    to_parameters = to_parameters
  )
}

# The base value of every parameter the model declares, by stratum for a
# stratified one, whether it declares an entry per stratum or one entry
# holding its values by stratum.
base_values <- function(declared) {
  values <- list()
  for (p in declared) {
    v <- p$base.value
    if (!is.null(p$stratum)) {
      if (is.null(values[[p$name]])) values[[p$name]] <- list()
      values[[p$name]][[p$stratum]] <- v
    } else {
      values[p$name] <- list(v)
    }
  }
  values
}
