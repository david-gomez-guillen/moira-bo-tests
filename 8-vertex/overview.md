# High-grade lesion natural history

A cohort model of the natural history of a cancer that arises from a precancerous
high-grade lesion, in a population with no screening. It follows a single birth cohort
from age 20, lesion-free, to age 85 in annual cycles.

## States and transitions

- **Lesion-free**: each year a person may develop a high-grade lesion
  (`p.lesion.onset`).
- **High-grade lesion**: each year a lesion may regress, returning the person to the
  lesion-free state (`p.lesion.regress`), or become a preclinical cancer
  (`p.lesion.progress`). A person whose lesion regressed can develop another.
- **Preclinical cancer**: an undiagnosed cancer, diagnosed once it causes symptoms
  (`p.symptoms`, 0.25 a year).
- **Diagnosed cancer**: followed for costs and quality of life, with an annual risk of
  cancer death.

Everyone alive also faces other-cause mortality, which rises with age along a Gompertz
curve.

The three natural history probabilities depend on age, with one value for each of the
age groups 20-34, 35-49, 50-64 and 65-85.

## Calibration

The `constrained` scheme calibrates onset, regression and progression in every age
group, twelve values in all, against three sources:

- **Cancer incidence**: registry incidence of diagnosed cancer per 100,000 person-years
  among the undiagnosed, in eight age groups.
- **Lesion prevalence**: the share of undiagnosed people with a high-grade lesion, from a
  population study that examined people of every age group.
- **Lesion persistence**: the share of untreated lesions still present when examined
  again two years after they were found, from a follow-up cohort in each age group.

The error is the relative squared error, averaged within each source and summed over the
three.

The search is restricted to natural histories the clinical literature considers credible:

- in every age group, a lesion lasts between 3 and 6 years on average before it regresses
  or progresses;
- in every age group, at most a quarter of lesions progress rather than regress;
- the lifetime risk of developing a lesion between 20 and 85, in the absence of death, is
  at most 32%.
