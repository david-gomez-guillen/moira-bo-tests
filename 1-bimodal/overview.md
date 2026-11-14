# Natural history of a cancer with two pathways

A Markov cohort model of a cancer that can develop along two pathways. A birth cohort is followed in annual cycles from age 20 to 85, and the whole cohort is free of lesions at 20.

## Pathways

- **Slow pathway**: a low-grade lesion that can become high-grade, and a high-grade lesion that can become a preclinical cancer. Both lesions can regress: a low-grade lesion to healthy, a high-grade lesion to low-grade.
- **Fast pathway**: a short-lived lesion that becomes a preclinical cancer within a few years or clears. Only part of the cohort carries the susceptibility it needs, and clearing a lesion of either pathway leaves a person no longer at risk of it.

Both pathways end in the same **preclinical cancer**, which each year can cause symptoms and be diagnosed as a **clinical cancer**. After diagnosis there is an annual risk of cancer death. Other-cause mortality follows a Gompertz curve and is the same in every state.

## Parameters

The model is stratified by age group: 20-39, 40-54, 55-69 and 70-84. Three natural-history probabilities take one value per group, applied to every year of age in it:

- **Low-grade lesion onset** (`p.lgl.onset`): the annual probability that a healthy person develops a low-grade lesion.
- **Fast-pathway lesion onset** (`p.fast.onset`): the annual probability that a susceptible healthy person develops a fast-pathway lesion.
- **High-grade lesion progression** (`p.hgl.progress`): the annual probability that a high-grade lesion becomes a preclinical cancer.

The other natural-history parameters take one value at every age: progression and regression of low-grade lesions, regression of high-grade lesions, the share of the cohort susceptible to the fast pathway, progression and clearance of fast-pathway lesions, and the annual probability that a preclinical cancer causes symptoms.

## Calibration targets

All by age group, each over the years of age in the group.

- **Clinical incidence** per 100,000 person-years of the undiagnosed, in two age groups per stratum, from the cancer registry. Cancer is too rare before 30 for the registry to give a rate.
- **Precursor prevalence** in one age group per stratum: the average share of the undiagnosed who have a lesion of either pathway, from a screening study that examined a sample of each age group.

The calibration scheme fits the three stratified probabilities. The error is the relative squared error, averaged within each series and summed over the two.

## Costs and effects

A clinical diagnosis costs a one-off amount, and every year lived after it has a lower utility. Both are discounted to age 20. They play no part in the calibration.
