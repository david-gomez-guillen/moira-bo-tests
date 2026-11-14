# Natural history of an HPV-related cancer by sex and age

A Markov cohort model of the natural history of a cancer caused by a sexually
acquired infection, in men and women without screening or vaccination. Two
closed cohorts, one of women and one of men, are followed from age 20 to 85 in
annual cycles. Both start healthy.

## States

- **Healthy**: no persistent infection and no lesion.
- **Low-grade lesion**: a persistent infection with a low-grade lesion. It can
  clear, back to healthy, or progress to a high-grade lesion.
- **High-grade lesion**: the precursor of cancer. It can regress to a low-grade
  lesion or become invasive.
- **Preclinical cancer**: an invasive cancer that has not caused symptoms yet.
  It is diagnosed through symptoms with an annual probability of 30%.
- **Diagnosed cancer**: where a cohort member leaves the population at risk.

Every living state is also exposed to other-cause death, which follows a
Gompertz curve in age, higher for men than for women.

## Strata

Twelve strata, sex by age group: `F.20_29`, `F.30_39`, `F.40_49`, `F.50_59`,
`F.60_69`, `F.70_85`, then the same six for men (`M.`). A stratified parameter
takes one value per stratum, and a cohort member uses the value of their sex and
of the age group they are in that year.

## Calibrated parameters

Five annual transition probabilities, each with a value per stratum, so 60
values in all:

| Parameter | Transition |
|---|---|
| `p.lesion.onset` | healthy to low-grade lesion |
| `p.lesion.clearance` | low-grade lesion to healthy |
| `p.lesion.progression` | low-grade to high-grade lesion |
| `p.high.grade.regression` | high-grade to low-grade lesion |
| `p.invasion` | high-grade lesion to preclinical cancer |

The calibration vector holds them in that order, each as a block of twelve
values in the order of the strata above.

## Targets

Two series, each by sex and age group, 24 values in all:

- **Cancer incidence**: diagnosed cancers per person-year among those alive and
  never diagnosed, as a cancer registry reports it.
- **Lesion prevalence**: the share of those alive and never diagnosed who carry
  a lesion of either grade, as a population survey would find it, averaged over
  the years of the age group.

The error is the relative squared error, averaged within each series and summed
over the two.

## Calibration schemes

All three schemes calibrate the same 60 values against the same targets.

- **standard**: when the probabilities of leaving a state in a cycle add up to
  more than 1, they are scaled down so that together they move the whole state.
- **failures**: the same model without that correction. Exit probabilities that
  add up to more than 1 stop the simulation with an error, and the evaluation
  has no outputs.
- **constrained**: the model of `standard`, with every transition required to
  change smoothly with age. Within each sex, the values of neighbouring age
  groups differ by at most a limit of the transition's own scale, twice the
  largest such step of its prior:

  | Parameter | Largest step between neighbouring age groups |
  |---|---|
  | `p.lesion.onset` | 0.024 |
  | `p.lesion.clearance` | 0.012 |
  | `p.lesion.progression` | 0.0018 |
  | `p.high.grade.regression` | 0.8 |
  | `p.invasion` | 0.14 |

  That is 5 transitions × 2 sexes × 5 pairs of age groups, 50 constraints.
