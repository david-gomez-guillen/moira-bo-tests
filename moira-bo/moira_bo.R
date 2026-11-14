# MOIRA-BO, Multi-Output Inverse Region Acquisition: Bayesian optimization of
# the calibration of a simulation model whose error is a known function of its
# outputs.
#
# Near a fit each output is a smooth function of the parameters, while the
# error folds at every target. So a Gaussian process is fitted per output, and
# a candidate is scored by putting the predicted outputs through the error
# function. Each iteration gathers candidates in a trust region from a random
# sampler, a Latin hypercube of the region, the problem's feasibility sampler
# and, for a few iterations after each of its trainings, an inverse model that proposes parameters meeting every
# target. The candidate of highest expected improvement is simulated, and when
# none promises an improvement the predicted error is descended instead.
#
# The Gaussian processes and the inverse model are in moira_bo.py, which
# moira_bo() expects loaded as `py` (see load_moira_python()).
#
# A problem is a list of:
#   x0          the initial guess, a named vector
#   lower, upper  the box searched
#   evaluate(x)   runs the model: list(error =, output =), output a list or
#                 vector of the model's outputs
#   score(M)      the error of each row of M, rows of flattened outputs
#   target        the targets, flattened as the outputs are
#   feasible(X)   whether each row of X meets the constraints
#   feasibility_sampler(center, n, width)  n plausible points around center,
#                 a matrix, feasible by construction where it can be
#   guesses       optional further initial guesses, rows of a matrix

`%||%` <- function(a, b) if (is.null(a)) b else a

# The settings of the article, Table S1 of its supplementary material.
moira_settings <- function(...) {
  s <- list(
    n_init = 5L,              # n_0, further initial guesses drawn when the problem has none
    n_start = 5L,             # n_r, draws of the feasibility sampler around the best of those
    start_width = 1,          # w_r, their width
    full_every = 10L,         # T_f, iterations between full fits of the processes
    gp_starts = 5L,           # n_f, starts of a full fit
    gp_iter = 100L,           # L-BFGS-B iterations of a fit
    tr_init = 0.2, tr_min = 0.005, tr_max = 1.6,
    tr_success = 3L,          # improvements before the region doubles
    tr_fail = NULL,           # failures before it halves, max(4, ceiling(d / 2))
    n_candidates = 5000L,     # n_c, shared out among the samplers by their weights
    candidate_weights = c(random = 0.4,       # w_rand
                          feasibility = 0.4,  # w_feas
                          inverse = 0.2),     # w_inv, 0 trains no inverse model
    n_score = 256L,           # n_s, candidates whose expected improvement is estimated
    n_draws = 128L,           # n_d, joint draws of the outputs for it
    n_descent = 4L,           # n_gd, starts of the descent besides the best point
    inverse_every = 10L,      # T
    inverse_expire = 3L,      # E, 0 offers the samples until the next training
    latent_dim = NULL,        # q, max(2, d - K + 1)
    inverse_widths = c(0.05, 0.1, 0.2, 0.4),
    inverse_labelled = 2000L, # n_t, spread over the widths
    inverse_weight = 5,       # omega
    inverse_max_sd = 0.05,    # sigma_max
    inverse_epochs = 300L,    # e_0
    inverse_rounds = 3L,      # R
    inverse_round_samples = 1000L,  # n_rho
    inverse_round_epochs = 100L,    # e_1
    inverse_pool = 2000L      # n_q
  )
  overrides <- list(...)
  s[names(overrides)] <- overrides
  s
}

# Loads moira_bo.py through reticulate, as `py` in the global environment.
load_moira_python <- function(path) {
  py <<- reticulate::import_from_path("moira_bo", path = dirname(path), convert = TRUE)
  invisible(py)
}

moira_key <- function(x) paste(sprintf("%a", as.numeric(x)), collapse = ",")

lhs <- function(n, d) {
  m <- matrix(NA_real_, n, d)
  for (j in seq_len(d)) m[, j] <- sample((seq_len(n) - stats::runif(n)) / n)
  m
}

# The candidates of each sampler: n shared out in proportion to the weights, by
# largest remainder so that they add up to n. A weight of 0 leaves its sampler out.
candidate_counts <- function(n, weights) {
  w <- pmax(0, as.numeric(weights))
  if (sum(w) == 0) return(stats::setNames(rep(0L, length(w)), names(weights)))
  exact <- n * w / sum(w)
  counts <- floor(exact)
  top <- order(exact - counts, decreasing = TRUE)[seq_len(n - sum(counts))]
  counts[top] <- counts[top] + 1
  stats::setNames(as.integer(counts), names(weights))
}

clip_box <- function(X, lo, hi) {
  X <- matrix(as.numeric(X), ncol = length(lo))
  pmin(pmax(X, matrix(lo, nrow(X), length(lo), byrow = TRUE)),
       matrix(hi, nrow(X), length(hi), byrow = TRUE))
}

# Runs MOIRA-BO on `problem` until the run holds `budget` evaluations.
# Iteration t is seeded with seed + t, in R, numpy and torch, so that a run
# repeats exactly. `on_eval` is called with each new row of the history.
# Returns the history, a row per evaluation, and the best feasible one.
moira_bo <- function(problem, budget = 56L, seed = 1L, settings = moira_settings(),
                     on_eval = NULL) {
  s <- settings
  d <- length(problem$x0)
  lo <- as.numeric(problem$lower)
  hi <- as.numeric(problem$upper)
  span <- pmax(hi - lo, 1e-12)
  norm <- function(X) sweep(sweep(matrix(X, ncol = d), 2, lo, "-"), 2, span, "/")
  from_norm <- function(Z) sweep(sweep(matrix(Z, ncol = d), 2, span, "*"), 2, lo, "+")
  tr_fail <- s$tr_fail %||% max(4L, ceiling(d / 2))
  n_src <- candidate_counts(s$n_candidates, s$candidate_weights[c("random", "feasibility", "inverse")])

  run <- new.env()
  run$X <- matrix(numeric(0), 0, d, dimnames = list(NULL, names(problem$x0)))
  run$Y <- NULL
  run$error <- numeric(0)
  run$feasible <- logical(0)
  run$ok <- logical(0)
  run$source <- character(0)
  run$predicted <- numeric(0)

  # Simulates x and adds it to the history; a failed simulation is kept as one,
  # since it cost a simulation, and left out of what the processes are fitted on.
  evaluate <- function(x, source, predicted = NA_real_) {
    x <- as.numeric(x)
    res <- tryCatch(problem$evaluate(x), error = function(e) list(error = NA_real_))
    out <- if (!is.null(res$output)) unlist(res$output, use.names = TRUE)
    ok <- !is.null(res$error) && length(res$error) == 1 && is.finite(res$error) &&
      !is.null(out) && (is.null(run$Y) || length(out) == ncol(run$Y)) && all(is.finite(out))
    # The outputs are laid out by the first that succeeds, the failures before
    # it holding a row of NA as every failure does.
    if (is.null(run$Y) && ok)
      run$Y <- matrix(NA_real_, nrow(run$X), length(out), dimnames = list(NULL, names(out)))
    run$X <- rbind(run$X, x)
    run$Y <- if (!is.null(run$Y)) rbind(run$Y, if (ok) out else rep(NA_real_, ncol(run$Y)))
    run$error <- c(run$error, if (ok) res$error else NA_real_)
    run$feasible <- c(run$feasible, problem$feasible(matrix(x, 1)))
    run$ok <- c(run$ok, ok)
    run$source <- c(run$source, source)
    run$predicted <- c(run$predicted, predicted)
    if (is.function(on_eval)) on_eval(history_row(length(run$error)))
  }
  history_row <- function(i) {
    best <- best_so_far(i)
    data.frame(evaluation = i, source = run$source[i], error = run$error[i],
               feasible = run$feasible[i], predicted = run$predicted[i],
               best = if (is.null(best)) NA_real_ else run$error[best])
  }
  best_so_far <- function(n = length(run$error)) {
    e <- run$error[seq_len(n)]
    usable <- run$ok[seq_len(n)] & run$feasible[seq_len(n)]
    if (!any(usable)) usable <- run$ok[seq_len(n)]
    if (!any(usable)) return(NULL)
    which(usable)[which.min(e[usable])]
  }
  evaluated <- function() apply(run$X, 1, moira_key)

  # The start: the initial guess, the further guesses, then draws around the
  # best of them.
  set.seed(seed)
  py$set_seed(seed)
  evaluate(problem$x0, "initial guess")
  guesses <- problem$guesses
  if (is.null(guesses) && s$n_init > 0)
    guesses <- problem$feasibility_sampler(problem$x0, s$n_init, s$start_width)
  for (i in seq_len(NROW(guesses)))
    if (length(run$error) < budget) evaluate(clip_box(guesses[i, ], lo, hi), "guess")
  if (is.null(best_so_far())) stop("the model failed on every point of the start", call. = FALSE)
  start <- clip_box(problem$feasibility_sampler(run$X[best_so_far(), ], s$n_start, s$start_width), lo, hi)
  for (i in seq_len(nrow(start)))
    if (length(run$error) < budget) evaluate(start[i, ], "start")

  gps <- NULL
  logged_before <- NULL
  tr_length <- s$tr_init
  successes <- 0L
  failures <- 0L
  pool <- NULL
  trained_at <- NA_integer_
  iteration <- 0L

  while (length(run$error) < budget) {
    iteration <- iteration + 1L
    set.seed(seed + iteration)
    py$set_seed(seed + iteration)

    # The surrogate: a GP per output, on its logarithm where it is positive.
    ok <- run$ok
    Y <- run$Y[ok, , drop = FALSE]
    K <- ncol(Y)
    logged <- apply(Y > 0, 2, all)
    to_scale <- function(M) { M[, logged] <- log(M[, logged]); M }
    from_scale <- function(M) { M[, logged] <- exp(M[, logged]); colnames(M) <- colnames(run$Y); M }
    full <- is.null(gps) || (iteration - 1L) %% s$full_every == 0L || !identical(logged, logged_before)
    gps <- py$fit_gps(gps, norm(run$X[ok, , drop = FALSE]), to_scale(Y), full,
                      starts = s$gp_starts, n_iter = s$gp_iter)
    logged_before <- logged
    predict <- function(X) {
      p <- py$predict_gps(gps, norm(X))
      list(mu = matrix(p[[1]], ncol = K), sd = matrix(p[[2]], ncol = K))
    }
    score <- function(M) { f <- problem$score(from_scale(M)); f[!is.finite(f)] <- Inf; f }

    best <- best_so_far()
    center <- as.numeric(run$X[best, ])
    usable <- run$ok & run$feasible
    incumbent <- if (any(usable)) min(run$error[usable]) else Inf

    # The trust region, its sides following the lengthscales of the outputs.
    log_ls <- as.numeric(py$mean_log_lengthscales(gps))
    sides <- tr_length * pmin(pmax(exp(log_ls - mean(log_ls)), 0.1), 10)
    c_norm <- as.numeric(norm(center))
    tr_lo <- pmax(0, c_norm - sides / 2)
    tr_hi <- pmin(1, c_norm + sides / 2)
    region <- function(n) from_norm(sweep(sweep(lhs(n, d), 2, tr_hi - tr_lo, "*"), 2, tr_lo, "+"))

    cand <- if (n_src[["random"]] > 0) region(n_src[["random"]]) else matrix(numeric(0), 0, d)
    source <- rep("random", nrow(cand))
    if (n_src[["inverse"]] > 0 && s$inverse_every > 0 && (iteration - 1L) %% s$inverse_every == 0L) {
      pool <- inverse_pool(run, problem, predict, logged, lo, hi, s,
                           latent_dim = s$latent_dim %||% max(2L, d - K + 1L),
                           seed = 100000L + iteration)
      trained_at <- iteration
    }
    if (s$inverse_expire > 0 && !is.na(trained_at) && iteration - trained_at >= s$inverse_expire)
      pool <- NULL
    if (NROW(pool)) {
      pool <- pool[!apply(pool, 1, moira_key) %in% evaluated(), , drop = FALSE]
      take <- sample(nrow(pool), min(n_src[["inverse"]], nrow(pool)))
      cand <- rbind(cand, pool[take, , drop = FALSE])
      source <- c(source, rep("inverse", length(take)))
    }
    if (n_src[["feasibility"]] > 0) {
      drawn <- clip_box(problem$feasibility_sampler(center, n_src[["feasibility"]], min(1, tr_length)), lo, hi)
      inside <- apply(norm(drawn), 1, function(u) all(u >= tr_lo & u <= tr_hi))
      cand <- rbind(cand, drawn[inside, , drop = FALSE])
      source <- c(source, rep("feasibility", sum(inside)))
    }
    fresh <- !apply(cand, 1, moira_key) %in% evaluated()
    cand <- cand[fresh, , drop = FALSE]
    source <- source[fresh]
    if (!nrow(cand)) {
      cand <- region(1L)
      source <- "fallback"
    }
    feasible <- problem$feasible(cand)
    keep <- if (any(feasible)) which(feasible) else seq_len(nrow(cand))

    # The candidates of lowest predicted error, and the expected improvement of
    # each over common draws of the outputs.
    p <- predict(cand[keep, , drop = FALSE])
    predicted <- score(p$mu)
    top <- order(predicted)[seq_len(min(s$n_score, length(keep)))]
    keep <- keep[top]
    predicted <- predicted[top]
    mu <- p$mu[top, , drop = FALSE]
    sd <- p$sd[top, , drop = FALSE]
    ei <- 0
    if (is.finite(incumbent)) {
      eps <- matrix(stats::rnorm(s$n_draws * K), nrow = s$n_draws)
      draws <- do.call(rbind, lapply(seq_len(s$n_draws), function(j) mu + sweep(sd, 2, eps[j, ], "*")))
      f <- matrix(score(draws), nrow = length(keep))
      ei <- rowMeans(pmax(incumbent - f, 0))
      pick <- if (max(ei) > 0) which.max(ei) else which.min(apply(f, 1, min))
    } else {
      pick <- which.min(predicted)
    }
    x_new <- cand[keep[pick], ]
    picked <- source[keep[pick]]
    picked_predicted <- predicted[pick]

    # With no candidate promising an improvement, the predicted error is
    # descended within the region from the best point and the best candidates.
    if (is.finite(incumbent) && max(ei) <= 0) {
      predicted_at <- function(Z) score(predict(from_norm(Z))$mu)
      h <- 1e-7
      descend <- function(z0) tryCatch(stats::optim(
        z0, function(z) predicted_at(z),
        function(z) {
          step <- ifelse(z + h <= tr_hi, h, -h)
          v <- predicted_at(rbind(z, t(z + diag(step, d))))
          (v[-1] - v[1]) / step
        },
        method = "L-BFGS-B", lower = tr_lo, upper = tr_hi,
        control = list(maxit = 100, factr = 10))$par, error = function(e) z0)
      starts <- rbind(c_norm, norm(cand[keep[order(predicted)[seq_len(min(s$n_descent, length(keep)))]], ,
                                        drop = FALSE]))
      U <- from_norm(do.call(rbind, lapply(seq_len(nrow(starts)), function(i) descend(starts[i, ]))))
      U <- U[problem$feasible(U), , drop = FALSE]
      U <- U[!apply(U, 1, moira_key) %in% evaluated(), , drop = FALSE]
      if (nrow(U)) {
        descended <- score(predict(U)$mu)
        if (min(descended) < incumbent) {
          x_new <- U[which.min(descended), ]
          picked <- "descent"
          picked_predicted <- min(descended)
        }
      }
    }

    evaluate(x_new, picked, picked_predicted)
    n <- length(run$error)
    if (run$ok[n] && run$feasible[n] && run$error[n] < incumbent) {
      successes <- successes + 1L
      failures <- 0L
    } else {
      failures <- failures + 1L
      successes <- 0L
    }
    if (successes >= s$tr_success) {
      tr_length <- min(2 * tr_length, s$tr_max)
      successes <- 0L
    } else if (failures >= tr_fail) {
      tr_length <- tr_length / 2
      failures <- 0L
    }
    if (tr_length < s$tr_min) tr_length <- s$tr_init
  }

  history <- do.call(rbind, lapply(seq_along(run$error), history_row))
  history <- cbind(history, as.data.frame(run$X, optional = TRUE))
  if (!is.null(run$Y)) history <- cbind(history, as.data.frame(run$Y, optional = TRUE))
  best <- best_so_far()
  list(history = history, best = stats::setNames(as.numeric(run$X[best, ]), names(problem$x0)),
       best_error = run$error[best], best_feasible = run$feasible[best])
}

# The inverse model's pool: samples of a conditional VAE of the parameters
# given the misses of the outputs, drawn where every miss is 0. It is trained
# on the run's evaluations, weighted up, and on draws of the feasibility sampler
# around the best labelled with the outputs the surrogate predicts, kept where
# the surrogate is sure of them. Each of a few rounds adds its own samples,
# labelled alike, so that its capacity goes where the targets are met.
inverse_pool <- function(run, problem, predict, logged, lo, hi, s, latent_dim, seed) {
  ok <- run$ok
  X <- run$X[ok, , drop = FALSE]
  Y <- run$Y[ok, , drop = FALSE]
  # The targets in the order of the outputs: by name where they carry the
  # names of the outputs, by position otherwise.
  target <- problem$target
  target <- as.numeric(if (!is.null(names(target)) && all(colnames(Y) %in% names(target)))
                         target[colnames(Y)] else target)
  ratio <- logged & target > 0
  scale <- pmax(abs(target), 1e-12)
  # The miss of an output: log(m / t) where both are positive, the relative
  # miss otherwise; Mt is on the scale each output is modelled on.
  miss <- function(Mt) {
    Mt <- as.matrix(Mt)
    out <- Mt
    for (k in seq_along(target)) {
      m <- if (logged[k]) exp(Mt[, k]) else Mt[, k]
      out[, k] <- if (ratio[k]) Mt[, k] - log(target[k]) else (m - target[k]) / scale[k]
    }
    out
  }
  miss_sd <- function(Mt, Sd) {
    out <- as.matrix(Sd)
    for (k in seq_along(target)) if (!ratio[k])
      out[, k] <- (if (logged[k]) exp(Mt[, k]) * Sd[, k] else Sd[, k]) / scale[k]
    out
  }
  label <- function(M) {
    if (!NROW(M)) return(list(X = M, C = matrix(numeric(0), 0, length(target))))
    p <- predict(M)
    sure <- apply(abs(miss_sd(p$mu, p$sd)), 1, max) <= s$inverse_max_sd
    list(X = M[sure, , drop = FALSE], C = miss(p$mu[sure, , drop = FALSE]))
  }
  observed <- Y
  observed[, logged] <- log(observed[, logged])

  usable <- run$feasible[ok]
  f <- run$error[ok]
  center <- X[which.min(ifelse(usable | !any(usable), f, Inf)), ]
  per_width <- ceiling(s$inverse_labelled / length(s$inverse_widths))
  draws <- do.call(rbind, lapply(s$inverse_widths, function(w)
    clip_box(problem$feasibility_sampler(center, per_width, w), lo, hi)))
  lab <- label(draws)
  pool_X <- rbind(X, lab$X)
  pool_C <- rbind(miss(observed), lab$C)
  pool_W <- c(rep(s$inverse_weight, nrow(X)), rep(1, nrow(lab$X)))

  zero <- rep(0, length(target))
  model <- py$train_cvae(pool_X, pool_C, pool_W, latent_dim = as.integer(latent_dim),
                         epochs = as.integer(s$inverse_epochs), seed = as.integer(seed))
  for (r in seq_len(s$inverse_rounds)) {
    S <- clip_box(py$sample_cvae(model, zero, as.integer(s$inverse_round_samples),
                                 seed = as.integer(seed + r)), lo, hi)
    lab <- label(S[problem$feasible(S), , drop = FALSE])
    pool_X <- rbind(pool_X, lab$X)
    pool_C <- rbind(pool_C, lab$C)
    pool_W <- c(pool_W, rep(1, nrow(lab$X)))
    model <- py$train_cvae(pool_X, pool_C, pool_W, epochs = as.integer(s$inverse_round_epochs),
                           seed = as.integer(seed + r), model = model)
  }
  S <- clip_box(py$sample_cvae(model, zero, as.integer(s$inverse_pool), seed = as.integer(seed + 100)),
                lo, hi)
  S[problem$feasible(S), , drop = FALSE]
}
