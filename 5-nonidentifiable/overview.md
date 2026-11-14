# Two-pathway cancer natural history

A Markov cohort model of the natural history of a cancer, in annual cycles. A
cohort enters at age 30 and is followed to age 89, with no screening or other
intervention.

## Disease states

- **Healthy**: no lesion.
- **Low-grade lesion** and **high-grade lesion**: the precursor lesions of the
  slow pathway. A low-grade lesion progresses to high grade or regresses to
  healthy. A high-grade lesion progresses to a preclinical cancer or regresses
  to low grade.
- **Preclinical cancer, slow pathway**: an undiagnosed cancer that arose from a
  high-grade lesion.
- **Preclinical cancer, fast pathway**: an undiagnosed cancer that arose
  directly from healthy tissue, with no detectable precursor lesion. It causes
  symptoms sooner than a slow-pathway cancer.
- **Diagnosed cancer**: a cancer found because of its symptoms. Its prognosis
  does not depend on the pathway it came from. Each year it is cured, kills, or
  stays under treatment.
- **Cured**, **cancer death** and **other-cause death**.

At age 30 the cohort already carries some low-grade and high-grade lesions.
Other-cause mortality follows a Gompertz curve and applies to every living
state.

## Age groups

The four natural-history probabilities that are calibrated take one value per
age group: 30-44, 45-59, 60-74 and 75-89.

- `p.lgl.onset`: healthy to low-grade lesion.
- `p.lgl.progress`: low-grade to high-grade lesion.
- `p.hgl.progress`: high-grade lesion to preclinical cancer.
- `p.fast.onset`: healthy to preclinical cancer, the fast pathway.

Everything else (regression, time to symptoms, survival after diagnosis, the
lesions at entry, other-cause mortality) is fixed at literature values.

## Calibration targets

All by age group, from `data/targets.csv`:

- **Clinical incidence**: diagnoses per person-year among those never
  diagnosed, as a cancer registry reports it.
- **Cancer mortality**: cancer deaths per person-year among the living.
- **Lesion prevalence**: the share of the never-diagnosed carrying a lesion of
  either grade, averaged over the ages of the group. It comes from a survey
  that examined people for lesions without recording their grade.

The error is the relative squared error, averaged within each series and
summed over the three. Cancer mortality is two orders of magnitude below
incidence and four below lesion prevalence, so only a scale-free error lets
every series count.

## Calibration schemes

- **standard**: the four natural-history probabilities in every age group.
- **constrained**: the same, keeping the fast pathway between 10% and 70% of new
  cancers in every age group, which is what is known of it apart from the targets.
  The share is that of the healthy and lesion states of the cohort, followed from
  age 30 with the calibrated probabilities.
