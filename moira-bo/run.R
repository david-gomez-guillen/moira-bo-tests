# Runs MOIRA-BO on a problem of this repository, or on any folder holding a
# model written for THALASSA, and writes the run to a folder of results.
#
#   Rscript moira-bo/run.R <problem> [options]
#
#   <problem>          a problem of the repository, with or without its number
#                      (1-bimodal or bimodal), or the path of a folder holding a
#                      thalassa_interface.R
#   --scheme NAME      the calibration scheme (constrained)
#   --budget N         evaluations of the run, its start included (56)
#   --seed N           the seed; iteration t is seeded with seed + t (1)
#   --scope X          the box, the initial guess times 1 +/- scope (1)
#   --guesses FILE     further initial guesses, a CSV with a column per
#                      parameter named as in the history; without it, draws of
#                      the problem's feasibility sampler
#   --threads N        threads of PyTorch and the linear algebra (1)
#   --out DIR          where the results go (moira-bo/results/<problem>-seed<N>)
#   --quiet            no line per evaluation

args <- commandArgs(trailingOnly = TRUE)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
usage <- function() {
  lines <- readLines(file.path(here, "run.R"))
  cat(sub("^# ?", "", lines[seq_len(which(!startsWith(lines, "#"))[1] - 1)]), sep = "\n")
  quit(status = 1)
}
if (!length(args) || args[1] %in% c("-h", "--help")) usage()

opt <- list(scheme = "constrained", budget = 56L, seed = 1L, scope = 1, guesses = NULL,
            threads = 1L, out = NULL, quiet = FALSE)
problem_arg <- args[1]
i <- 2L
while (i <= length(args)) {
  key <- sub("^--", "", args[i])
  if (identical(key, "quiet")) {
    opt$quiet <- TRUE
    i <- i + 1L
    next
  }
  if (!key %in% names(opt) || i == length(args)) usage()
  opt[[key]] <- args[i + 1L]
  i <- i + 2L
}
opt$budget <- as.integer(opt$budget)
opt$seed <- as.integer(opt$seed)
opt$scope <- as.numeric(opt$scope)
opt$threads <- as.integer(opt$threads)

# A problem of the repository is known by its folder or by its name without
# the number the folder starts with.
root <- dirname(here)
numbered <- list.files(root, pattern = paste0("^[0-9]+-", problem_arg, "$"))
dir <- if (dir.exists(problem_arg)) {
  problem_arg
} else if (length(numbered) == 1) {
  file.path(root, numbered)
} else {
  file.path(root, problem_arg)
}
if (!file.exists(file.path(dir, "thalassa_interface.R")))
  stop(sprintf("'%s' is neither a problem of the repository nor a folder holding a thalassa_interface.R",
               problem_arg), call. = FALSE)
dir <- normalizePath(dir)
out <- opt$out %||% file.path(here, "results", sprintf("%s-seed%d", basename(dir), opt$seed))

# The threads have to be set before Python and its libraries load, which read
# them only then.
for (v in c("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"))
  do.call(Sys.setenv, stats::setNames(list(as.character(opt$threads)), v))
prefix <- Sys.getenv("CONDA_PREFIX")
if (nzchar(prefix)) {
  python <- if (.Platform$OS.type == "windows") file.path(prefix, "python.exe")
            else file.path(prefix, "bin", "python")
  if (file.exists(python)) reticulate::use_python(python, required = TRUE)
}

source(file.path(here, "moira_bo.R"))
source(file.path(here, "thalassa_problem.R"))
load_moira_python(file.path(here, "moira_bo.py"))
reticulate::import("torch")$set_num_threads(opt$threads)

problem <- thalassa_problem(dir, scheme = opt$scheme, scope = opt$scope)
if (!is.null(opt$guesses)) {
  guesses <- utils::read.csv(opt$guesses, check.names = FALSE)
  missing <- setdiff(names(problem$x0), names(guesses))
  if (length(missing))
    stop(sprintf("the guesses lack the parameters %s", paste(missing, collapse = ", ")), call. = FALSE)
  problem$guesses <- as.matrix(guesses[, names(problem$x0), drop = FALSE])
}

cat(sprintf("MOIRA-BO on %s, scheme '%s': %d parameters, %d targets, %d evaluations, seed %d\n",
            problem$name, opt$scheme, length(problem$x0), length(problem$target), opt$budget, opt$seed))
started <- Sys.time()
on_eval <- if (!opt$quiet) function(row)
  cat(sprintf("%4d  %-13s  error %-11s  %-10s  best %s  (%s)\n", row$evaluation, row$source,
              formatC(row$error, format = "e", digits = 3), if (row$feasible) "feasible" else "infeasible",
              formatC(row$best, format = "e", digits = 3),
              format(round(difftime(Sys.time(), started, units = "secs")))))
res <- moira_bo(problem, budget = opt$budget, seed = opt$seed, on_eval = on_eval)

dir.create(out, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(res$history, file.path(out, "history.csv"), row.names = FALSE)
utils::write.csv(data.frame(parameter = names(res$best), value = res$best),
                 file.path(out, "best.csv"), row.names = FALSE)
cat(sprintf("\nBest error %s (%s), after %s. Written to %s\n",
            formatC(res$best_error, format = "e", digits = 4),
            if (res$best_feasible) "feasible" else "infeasible",
            format(round(difftime(Sys.time(), started, units = "mins"), 1)), out))
