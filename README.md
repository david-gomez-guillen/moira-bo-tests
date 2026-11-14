# MOIRA-BO test problems

Eight constrained calibration problems for testing MOIRA-BO (Multi-Output Inverse Region Acquisition), a Bayesian optimization method for calibrating expensive simulation models. The problems run in [THALASSA](https://thalassa.iconcologia.net), a web application for the cost-effectiveness analysis of simulation models, which offers MOIRA-BO among its calibration methods.

Each problem is a disease simulation model of the kind used in cost-effectiveness analysis. Each one is built to present a single difficulty to a calibration method.

## The problems

| Problem | Parameters | Targets | Constraints | Difficulty |
|---|---|---|---|---|
| `1-bimodal` | 12 | 12 | 1 nonlinear, conditional on the fast-pathway share at ages 70 to 84 | two basins, with the starting point in the worse one |
| `2-highdim` | 60 | 24 | 50 linear: neighbouring age groups differ by a bounded step | 60 parameters |
| `3-misspec` | 12 | 14 | 8 bounds on the mean lesion duration | a misspecified model whose targets cannot all be met |
| `4-noisy` | 16 | 16 | 4 linear: regression plus progression at most 0.18 | a microsimulation of 100,000 people, whose outputs carry Monte Carlo noise |
| `5-nonidentifiable` | 16 | 12 | 8 nonlinear: fast-pathway share of new cancers between 10% and 70% | parameters identified only in combination, with fewer targets than parameters |
| `6-threshold` | 16 | 16 | 1 nonlinear: basic reproduction number at least 1.05 | a transmission model whose prevalence is 0 below the threshold |
| `7-valley` | 21 | 18 | 19 linear: each age block non-decreasing, fast-pathway onset at most 35% of the low-grade one | two pathways that only one series separates, with a set of fits of dimension 3 |
| `8-vertex` | 12 | 16 | 13 nonlinear and coupled: lesion duration, share progressing, lifetime risk | an optimum at a vertex of the feasible set |

All the targets are synthetic. On `3-misspec` they come from a process the model cannot reproduce, with measurement error added. On `4-noisy` they are the expected outputs at a set of true parameters, which the noisy outputs of the microsimulation can only approach. On the other problems they are the outputs of the model at a set of true parameters. Those parameters are not part of this repository.

## How the repository is organized

- **Folders**: each problem is a numbered folder, `1-bimodal` to `8-vertex`, holding its model (`model.R`), the interface THALASSA drives it through (`thalassa_interface.R`), the text of its Overview tab (`overview.md`) and, for all but `7-valley`, its targets (`data/targets.csv`).
- **Branches**: each problem also has a branch named after it, without the number: `bimodal` to `vertex`. A branch adds a `thalassa_interface.R` at the root that loads the folder of its problem. THALASSA offers the branches of a model repository as its versions, so picking a version picks a problem.
- **The default branch** is `bimodal`, so loading the repository starts on that problem.
- **`moira-bo/`** holds a standalone version of MOIRA-BO, which runs the problems without THALASSA. `environment.yml` is the conda environment it runs in, and `run_moira_bo.sh` creates that environment and runs it.

## Running MOIRA-BO in THALASSA

1. Open [thalassa.iconcologia.net](https://thalassa.iconcologia.net) and press **Explore the demo**. The problems are among the models the demo offers, so no account is needed.
2. In the model bar at the top of the page, select **MOIRA-BO tests**, and pick a problem from its list of versions. The first load builds the environment the model runs in, which takes a minute or two.
3. In the **Calibration** tab, choose the `constrained` scheme and **MOIRA-BO** under *Bayesian optimization* in the method dropdown, then press **Run**.

The **Overview** tab of each problem describes the model and what makes its calibration hard. To compare methods on one problem, run them one after another in the Calibration tab: each new run keeps the previous ones on the plots. Nelder–Mead, SCBO and the agentic method are the ones MOIRA-BO was compared with. **Export** saves every run on the plots to a spreadsheet, which **Import** restores later.

## Running MOIRA-BO without THALASSA

[`moira-bo/`](moira-bo/) is a standalone version of the method, in R and Python, which runs on any of the problems from the command line. With [conda](https://conda-forge.org/download/) installed, from the root of the repository:

```bash
./run_moira_bo.sh
```

The first time, this creates the conda environment of [`environment.yml`](environment.yml), named `moira-bo`, which takes a few minutes. It then calibrates `1-bimodal` to 56 evaluations, printing a line per evaluation and writing the run to `moira-bo/results/`. Another problem, and any option of the method, go after it:

```bash
./run_moira_bo.sh 7-valley --seed 3
```

Without the script, the same is:

```bash
conda env create -f environment.yml
conda activate moira-bo
Rscript moira-bo/run.R 1-bimodal
```

The [README of `moira-bo/`](moira-bo/README.md) explains the options, what a run writes and how this version differs from the runs of the article.
