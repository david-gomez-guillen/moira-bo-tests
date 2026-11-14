# MOIRA-BO, standalone

A standalone version of MOIRA-BO (Multi-Output Inverse Region Acquisition), the Bayesian optimization method for calibrating expensive simulation models that the problems of this repository were built to test. It runs the problems from the command line, with no THALASSA and no calibration engine, and it can run any model written for THALASSA.

It follows the algorithm and the settings of the article. It is not the implementation the article's results come from, which is the calibration engine of THALASSA, so a run here does not reproduce one of those evaluation by evaluation (see [Differences from the article's runs](#differences-from-the-articles-runs)).

## Setting up

The method runs in R, the language of the models, and fits its Gaussian processes and its inverse model in Python with PyTorch, through reticulate. One conda environment holds both, [`environment.yml`](../environment.yml) at the root of the repository. It needs no GPU.

From the root of the repository, the script `run_moira_bo.sh` creates the environment the first time, which takes a few minutes, and runs the method:

```bash
./run_moira_bo.sh                      # 1-bimodal, to 56 evaluations
./run_moira_bo.sh 7-valley --seed 3    # another problem, with the options below
```

Without the script:

```bash
conda env create -f environment.yml
conda activate moira-bo
```

## Running

From the root of the repository, with the environment active:

```bash
Rscript moira-bo/run.R 1-bimodal
```

This calibrates the `constrained` scheme of `1-bimodal` to 56 evaluations, printing a line per evaluation:

```
MOIRA-BO on 1-bimodal, scheme 'constrained': 12 parameters, 12 targets, 56 evaluations, seed 1
   1  initial guess  error 9.491e-02    feasible    best 9.491e-02  (0 secs)
   ...
  12  inverse        error 2.995e-02    feasible    best 2.995e-02  (12 secs)
```

The second column says where the point came from: the start (`initial guess`, `guess`, `start`), or the candidate the iteration chose, from the `random` sampler, the `feasibility` sampler, the `inverse` model or the `descent`.

| Option | Meaning | Default |
|---|---|---|
| `--scheme NAME` | the calibration scheme | `constrained` |
| `--budget N` | evaluations of the run, its start included | 56 |
| `--seed N` | the seed; iteration *t* is seeded with seed + *t* in R, NumPy and PyTorch, so a run repeats exactly | 1 |
| `--scope X` | the box searched: the initial guess times 1 ± *X*, within [0, 1] | 1 |
| `--guesses FILE` | further initial guesses, a CSV with a column per parameter, named as in the history | draws of the feasibility sampler |
| `--threads N` | threads of PyTorch and the linear algebra | 1 |
| `--out DIR` | where the results go | `moira-bo/results/<problem>-seed<N>` |
| `--quiet` | no line per evaluation | |

A problem of the repository can be named with or without its number, `1-bimodal` or `bimodal`. It can also be any folder holding a `thalassa_interface.R`, given by its path.

### What it writes

- `history.csv`: a row per evaluation, with where the point came from, its error, whether it is feasible, the error the surrogate predicted for it, the best feasible error so far, its parameters and the outputs the model gave.
- `best.csv`: the parameters of the best feasible evaluation.

The parameters are named `<parameter>.<stratum>`, as in THALASSA's Calibration tab.

## How it works

1. **The start.** The initial guess, the further initial guesses, and five draws of the problem's feasibility sampler around the best of them: 11 evaluations, on which the surrogate is first fitted.
2. **The surrogate.** A Gaussian process per output of the model, on its logarithm where it is positive, refitted at every iteration. A candidate is scored by putting the predicted outputs through the error function of the problem, so only the outputs are modelled.
3. **The trust region**, TuRBO's, centred on the best feasible evaluation and shaped by the lengthscales of the processes.
4. **The candidates**, 5000 of them, shared out among three samplers by their weights: 40% from the random sampler, a Latin hypercube of the region; 40% from the feasibility sampler, which draws plausible, feasible points around its centre; and 20% from the inverse model. A weight of 0 leaves a sampler out. The inverse model is a conditional variational autoencoder of the parameters given how far each output misses its target. It is trained every 10 iterations, and sampled where every miss is 0. Its samples are offered for 3 iterations after each training and then withdrawn, its share of the candidates left empty.
5. **The acquisition.** Of the feasible candidates, the 256 of lowest predicted error are kept, and the one of highest expected improvement over the best is simulated, estimated from 128 joint draws of the outputs. When none promises an improvement, the predicted error is descended by L-BFGS-B inside the region instead.

The settings are those of the article, listed with their symbols in `moira_settings()` of [moira_bo.R](moira_bo.R).

## Differences from the article's runs

- **The further initial guesses.** In the article they were five guesses of an agentic calibration method, a language model. Here they are five draws of the feasibility sampler, or the guesses given with `--guesses`. In the article's ablations, the agent's guesses gave a better start than random draws, so a run here usually starts further from a fit.
- **The implementation.** This version is written for reading, not as the engine the article ran. It makes the same choices, but draws its random numbers in its own order, so the same seed gives a different run.
- **The time.** In the article, a run of 56 evaluations took 10 to 30 minutes on one core, most of it fitting a Gaussian process per output.

## The files

- [run.R](run.R): the command line.
- [moira_bo.R](moira_bo.R): the method, on a problem given as a list of functions. It knows nothing about THALASSA, so it can calibrate any model whose error is a known function of its outputs.
- [thalassa_problem.R](thalassa_problem.R): reads a model written for THALASSA and turns one of its calibration schemes into that problem.
- [moira_bo.py](moira_bo.py): the Gaussian processes and the inverse model.
- [environment.yml](../environment.yml) and [run_moira_bo.sh](../run_moira_bo.sh), at the root of the repository: the conda environment, and the script that creates it and runs the method.

## Running it on another model

A model written for THALASSA runs as it is if its calibration scheme gives:

- an `error_function(parameters, target)` that returns the error and the outputs it was computed from, as `list(error =, output =)`;
- an `output_error(output, target)` that computes the same error from the outputs alone;
- a `latent_space_training_set(center, n, variation)`, its feasibility sampler, that draws `n` plausible parameter vectors around `center`;
- optionally `constraints(x)`, whose values are all at most 0 where `x` is feasible.

The parameters are searched within [0, 1], which suits probabilities and proportions. A model calibrated on other scales can be given its own bounds through the `lower` and `upper` of `thalassa_problem()`, or written as a problem of its own for `moira_bo()`: the comment at the top of [moira_bo.R](moira_bo.R) lists what a problem holds.
