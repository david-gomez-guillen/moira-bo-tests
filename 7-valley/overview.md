# Two-pathway cancer screening model

A cost-effectiveness model of screening for a generic cancer. Cancer arises
through **two parallel routes**, which is what makes the model difficult to calibrate: the targets constrain the *sum* of what the two
routes produce much better than they constrain the *split* between them, so the
error surface has several distinct local minima instead of one basin.

## Model structure

A Markov cohort model with annual cycles, followed from age 20 to age 89 (70
cycles), with no new entries.

**Slow pathway** (most cancers):

    Normal -> low-grade lesion -> high-grade lesion -> preclinical cancer (localized -> advanced)

**Fast pathway** (a minority of cancers):

    Normal -> fast-pathway lesion -> preclinical cancer (localized -> advanced)

Low-grade and fast-pathway lesions regress at a low annual rate; high-grade
lesions regress back to low-grade ones. Each pathway has its own pair of
preclinical (undiagnosed) cancer states, because fast-pathway cancers progress to
an advanced stage faster and stay asymptomatic longer than slow-pathway ones.
This is what makes the stage distribution at diagnosis carry information about
the pathway mix.

Preclinical cancers surface either through symptoms or through screening, and
become diagnosed cancer at the stage they had reached. Diagnosed patients are
treated and either recover (with a small annual risk of late recurrence and
death) or die of cancer. Everyone alive is exposed to a Gompertz other-cause
mortality that competes with all of the above.

The 13 states are `normal`, `lgl`, `hgl`, `fpl`, `pre.early.slow`,
`pre.late.slow`, `pre.early.fast`, `pre.late.fast`, `clin.early`, `clin.late`,
`survivor`, `dead.cancer` and `dead.other`.

Within a cycle, screening is applied first, then costs and utilities are
accrued, then the natural history transition. Rows of the transition matrix that
extreme parameter sets would push above 1 are rescaled, so an optimizer probing
the corners of the search box gets a valid Markov chain rather than an error.

`simulate()` also takes an artificial per-call delay
(`MODEL.SIMULATION.DELAY` in `model.R`), there to stand in for a slow model when
testing progress reporting or calibration budgets. It is 3 seconds on `main`
and 0 on the `no_delay` branch, where a run takes the few milliseconds the model
actually needs.

## Strategies

Two generic screening modalities are modelled: a non-invasive screening *test*,
whose positives go on to the *procedure*, and the *procedure* itself, which finds
and removes lesions.

| Strategy | Description |
|---|---|
| `no_screening` | Natural history only. This is the strategy the calibration is run on. |
| `test_biennial` | Screening test every two years from 50 to 74; positives go on to the procedure. |
| `procedure_10y` | Procedure every ten years from 50 to 79. |
| `procedure_45` | The same, starting at 45. |

Detected lesions are removed and the person returns to `normal`; detected
preclinical cancers are diagnosed at their current stage. Fast-pathway lesions
are harder to detect, so both modalities miss them far more often than
slow-pathway lesions, and fast-pathway cancers are detected somewhat less often
than slow-pathway ones (`rr.detection.fast`). Procedures carry a risk of a
serious complication, with its own cost and disutility. Surveillance after a
positive finding is not modelled.

## Parameters

Age-dependent parameters take one value per stratum. Everything else is a single
value. The `class` field only groups the parameters in the interface.

| Parameter | Base value | Class |
|---|---|---|
| `p.lgl.onset` | 0.0050 - 0.0370 by age | Natural history (slow pathway) |
| `p.lgl.progress` | 0.015 | Natural history (slow pathway) |
| `p.lgl.regress` | 0.015 | Natural history (slow pathway) |
| `p.hgl.progress` | 0.0053 - 0.0251 by age | Natural history (slow pathway) |
| `p.hgl.regress` | 0.005 | Natural history (slow pathway) |
| `p.fpl.onset` | 0.00018 - 0.00102 by age | Natural history (fast pathway) |
| `p.fpl.progress` | 0.040 | Natural history (fast pathway) |
| `p.fpl.regress` | 0.020 | Natural history (fast pathway) |
| `p.stage.progress.slow`, `p.stage.progress.fast` | 0.300, 0.450 | Cancer progression |
| `p.symptomatic.early.slow`, `p.symptomatic.late.slow` | 0.270, 0.550 | Cancer progression |
| `p.symptomatic.early.fast`, `p.symptomatic.late.fast` | 0.120, 0.450 | Cancer progression |
| `p.cure.early`, `p.cancer.death.early` | 0.550, 0.030 | Cancer survival |
| `p.cure.late`, `p.cancer.death.late` | 0.120, 0.250 | Cancer survival |
| `p.survivor.death` | 0.010 | Cancer survival |
| `mortality.other.base`, `mortality.other.rate` | 0.0005, 0.085 | General |
| `adherence.test`, `spec.test` | 0.650, 0.950 | Screening (test) |
| `sens.test.lgl`, `sens.test.hgl`, `sens.test.fpl`, `sens.test.cancer` | 0.040, 0.240, 0.080, 0.550 | Screening (test) |
| `adherence.procedure` | 0.550 | Screening (procedure) |
| `sens.procedure.lgl`, `.hgl`, `.fpl`, `.cancer` | 0.750, 0.920, 0.550, 0.950 | Screening (procedure) |
| `rr.detection.fast` | 0.850 | Screening (procedure) |
| `p.procedure.complication` | 0.002 | Screening (procedure) |
| `cost.test`, `cost.procedure`, `cost.removal`, `cost.complication` | 30, 800, 250, 3000 | Costs |
| `cost.treatment.early`, `cost.treatment.late`, `cost.followup` | 18000, 40000, 1200 | Costs |
| `utility.cancer.early`, `utility.cancer.late`, `utility.survivor` | 0.720, 0.500, 0.920 | Utilities |
| `disutility.complication` | 0.020 | Utilities |
| `discount` | 0.030 | General |

Sensitivity for advanced preclinical cancer is derived from the localized one
(times 1.15 for the test, 1.05 for the procedure, capped at 1) rather than being
a parameter of its own.

## Strata

Ten-year age groups from `20-29` to `80-89`.

## Outputs

- `summary`: one row per strategy with the total discounted cost (`C`) and the
  total discounted quality-adjusted life years (`E`) of the cohort.
- `outputs`: per strategy, the three series the calibration uses, each by
  stratum:
  - `Cancer incidence`, diagnosed cancers per person-year among those without a
    previous cancer diagnosis (screen-detected cancers included);
  - `Lesion prevalence`, the proportion of the same population carrying any
    lesion, of any type, as a study examining everyone for lesions would find it;
  - `Late-stage share`, the proportion of the cancers diagnosed in the stratum
    that were at an advanced stage.
- `incidence`: the incidence series alone, per strategy.
- `cohort.info`: the full cohort trace per strategy, one row per age.

## Calibration

Both schemes calibrate on `no_screening` against the same three target
series, which cover the six age groups from `30-39` to `80-89`. The `20-29`
group is burn-in: the cohort starts with no lesions, so its incidence is
essentially zero whatever the parameters. It is left out of the *targets* rather
than out of the scheme's strata, because the app builds the initial guess from
the base values over every stratum, so a scheme covering only some of them would
be handed a longer vector than it consumes.

| Scheme | Calibrated parameters | Dimension |
|---|---|---|
| `natural_history` | `p.lgl.onset`, `p.fpl.onset`, `p.hgl.progress` | 21 |
| `constrained` | the same three, with constraints | 21 |

The error is the mean squared *relative* deviation within each series, summed
over the three series. The series are on very different scales (incidence is per
mil, the other two are proportions), so a scale-free error is what makes all
three count.

The `constrained` scheme adds a `constraints` function, evaluated on the raw
calibration vector, for the constrained BO wrappers, and a
`constraints_description` function stating the same thing in words for the
agentic wrapper, which appends it to the agent's system prompt. Feasible
parameter sets are non-decreasing in age for all three parameters and keep total
fast-pathway lesion onset below 35% of total low-grade lesion onset; its training set generator only
produces feasible rows.
