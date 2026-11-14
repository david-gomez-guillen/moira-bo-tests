# Transmission model of a curable sexually transmitted infection

A deterministic compartmental model of a bacterial sexually transmitted
infection in a heterosexual population aged 15 to 64, used to estimate the
burden of infection that a testing or vaccination programme would act on.

## Structure

- **Population**: women and men in four age groups, 15-24, 25-34, 35-44 and
  45-64. The population is stationary: people enter at 15, move up an age group
  as they age, and leave the sexually active population at 65.
- **States**: susceptible, infected, and recovered. Recovery leaves antibodies and
  protection against reinfection, which wane over time and return the person
  to the susceptible state (an SIRS model).
- **Partnerships**: a fixed partnership matrix gives how many partnerships form
  each year between men and women of each pair of age groups at full sexual
  activity. Most are within an age group, and men's partners tend to be
  younger. The activity level of each age group scales the partnerships of
  that group, for its men and its women alike.
- **Transmission**: within a partnership with an infected partner, the
  infection passes on at each act with the per-act susceptibility of the
  susceptible partner. There are 60 acts a year, and men are 0.6 times as
  susceptible as women of the same age. A partnership passes the infection on
  only if it does so before it ends or the infection clears. A partnership ends
  when either partner moves on to a new one, so partnerships between very
  active people are short.
- **Clearance**: an infection ends by natural clearance, at a rate by age group,
  or by treatment, at 0.15 a year in women and 0.35 a year in men, who are more
  often symptomatic.

The model is integrated with forward-Euler steps of 0.05 years over 300
years, from 1% of every group infected, to reach its endemic equilibrium.

## Outputs

- Prevalence of current infection by age group and sex.
- Seroprevalence by age group and sex: the share carrying antibodies, which is
  the currently infected plus the recovered.

## Calibration

The `standard` scheme calibrates four parameters by age group, sixteen values
in all, against prevalence and seroprevalence surveys by age group and sex:

- `p.susceptibility`: per-act probability of acquiring the infection;
- `activity.level`: sexual activity relative to the partnership matrix;
- `r.clearance`: annual rate of natural clearance;
- `r.waning`: annual rate of losing antibodies and immunity.

The error is the relative squared error of each output, averaged within each
series and summed over the four series.

The `constrained` scheme calibrates the same values against the same targets,
restricted to parameters under which the infection is endemic: a basic
reproduction number of at least 1.05. R0 is the spectral radius of the
next-generation matrix at the infection-free state, from the susceptibility,
activity and clearance of every age group (`basic.reproduction.number()` in
`model.R`). Below R0 = 1 the infection dies out and every prevalence is 0.
