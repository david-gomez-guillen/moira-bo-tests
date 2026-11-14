# Natural history of a lesion-derived cancer

A Markov cohort model of a cancer that develops from a precancerous lesion. It follows a
cohort from age 20 to age 85 in yearly cycles, with no screening: every cancer is
diagnosed when it gives symptoms.

## States

- **Healthy**: no lesion. The whole cohort starts here at age 20.
- **Lesion**: a precancerous lesion. Each year it can regress to healthy or become a
  preclinical cancer.
- **Preclinical localized** and **preclinical advanced**: an undiagnosed cancer. A
  localized cancer either gives symptoms and is diagnosed at the localized stage, or
  advances without being noticed. An advanced cancer is diagnosed when it gives symptoms.
- **Diagnosed**: a cancer diagnosed after symptoms.
- **Other-cause death**: everyone alive runs the same Gompertz risk of dying of something
  else.

## Parameters

Three parameters depend on age, with one value per age stratum (20-34, 35-49, 50-64 and
65-84):

- `p.lesion.onset`: annual probability that a healthy person develops a lesion;
- `p.lesion.progress`: annual probability that a lesion becomes a preclinical cancer;
- `p.stage.progress`: annual probability that a preclinical localized cancer advances.

The rest are held at values from the literature: lesion regression (4% a year), the
annual probability of symptoms of a localized (15%) and of an advanced (55%) preclinical
cancer, and other-cause mortality.

## Calibration

The `standard` scheme calibrates the three age-dependent parameters, 12 values in all,
against three series of data on the same population:

- **Cancer incidence**: symptomatic diagnoses per 100,000 person-years among those never
  diagnosed, in eight age groups, from the cancer registry;
- **Lesion prevalence**: the share of people with a lesion in each age stratum, from a
  population study that examined a sample of every age group;
- **Late-stage share**: the share of diagnoses made at the advanced stage, at ages 30-54
  and 55-84, from the registry's staging records.

The `constrained` scheme calibrates the same values against the same data, keeping the
mean time a lesion lasts before it regresses or progresses, 1 / (regression + progression),
between 15.5 and 24 years in every age stratum, the range follow-up studies of the lesion
report.

The data are estimates and carry sampling error. The error of a parameter set is the
relative squared error of each output against its target, averaged within each series and
summed over the three series.
