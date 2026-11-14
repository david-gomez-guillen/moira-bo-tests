# Natural history of a progressive cancer

An individual-level microsimulation of a cancer that develops from a precancerous lesion. A birth cohort is followed in annual cycles from age 20 to 85.

## States

- **Healthy**: no lesion.
- **Lesion**: a precancerous lesion. Each year it can regress to healthy or progress to a preclinical cancer.
- **Preclinical cancer**: a cancer that has not yet caused symptoms. Each year it can cause symptoms and be diagnosed.
- **Clinical cancer**: a cancer diagnosed after symptoms. People stay here for the rest of the simulation.
- **Other-cause death**: from a fixed life table by five-year age band, the same in every state.

## Parameters

The model is stratified by age group: 20-39, 40-54, 55-69 and 70-84. Each natural-history probability takes one value per group, applied to every year of age in it.

- **Lesion onset**: the annual probability that a healthy person develops a lesion.
- **Lesion regression**: the annual probability that a lesion regresses.
- **Lesion progression**: the annual probability that a lesion becomes a preclinical cancer.
- **Symptoms**: the annual probability that a preclinical cancer causes symptoms and is diagnosed. Its inverse is roughly the mean time a cancer stays preclinical.

The whole cohort is healthy at 20.

## How it is simulated

Each evaluation simulates 100,000 people. Every year each living person draws one uniform number, which decides both death from other causes and the disease transition of that year. The outputs are Monte Carlo estimates: two runs with the same parameters give slightly different outputs unless the random state is the same.

## Calibration targets

All by age group.

- **Clinical incidence** per 100,000 person-years over the ages of the group, from the cancer registry. Person-years count everyone alive, diagnosed or not.
- **Lesion prevalence**: the average share of the living who have a lesion over the ages of the group, from a population study.
- **Preclinical cancer prevalence**: the average share of the living with a cancer that has not yet caused symptoms, from the same study.
- **Lesion persistence**: of the people found with a lesion at ages 25, 40, 55 and 70, the share of those still alive five years later who have a lesion then, from a follow-up of the study. A lesion that regressed and came back counts as present.

The error is the relative squared error, averaged within each of the four series and summed over them.

## Calibration schemes

- **standard**: the four natural-history probabilities in every age group, 16 values.
- **constrained**: the same, with the annual probability of leaving the lesion state, regression plus progression, at most 0.18 in every age group.

## Strategies

| Strategy | Description |
|---|---|
| `no_intervention` | Natural history only: cancers are diagnosed when they cause symptoms. The calibration is run on this strategy. |
| `screening_5y` | An examination at 50, 55, 60, 65 and 70. |

An examination comes at the start of the year. Each person has a draw for it, and with a lesion or a preclinical cancer is found when the draw is below the adherence (`screening.adherence`) times the sensitivity to what they have (`screening.sens.lesion`, `screening.sens.cancer`). A lesion found is removed, which leaves the person healthy, and a cancer found is diagnosed and counted in clinical incidence.

Every strategy of a run is simulated from the same random state, and the draws of the examinations come from a stream of their own, so the yearly draws of every person, and with them the deaths from other causes, are the same under every strategy. With these common random numbers the difference between two strategies is that of the strategies, not of the Monte Carlo noise, and the outputs before 50 are the same under both.

## Costs and effects

A diagnosis costs a one-off amount, and every year lived after it has a lower utility. An examination and the removal of a lesion each have a cost. All are discounted to age 20. They play no part in the calibration.
