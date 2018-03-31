import scipy.stats as stat
import numpy as np

weight = stat.binom.pmf(
    k=400,
    p=0.5,
    n=573
)
print("The likelihood that a clinical study reports 400 of 573 participants survived at the end of the 5-year study \n"
      "period if 50% of the patients in the simulated cohort survived beyond 5 years", weight)


