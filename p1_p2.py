import p1_SurvivalModel as Cls
#import scr.FigureSupport as Fig

MORTALITY_PROB = 0.1    # annual probability of mortality
TIME_STEPS = 100       # simulation length
REAL_POP_SIZE = 573    # size of the real cohort to make the projections for
NUM_SIM_COHORTS = 1000   # number of simulated cohorts used for making projections
ALPHA = 0.05            # significance level



simcohort= Cls.Cohort(
    id=1,
    pop_size=REAL_POP_SIZE,
    mortality_prob=MORTALITY_PROB
)

simcohort.simulate(TIME_STEPS)
print("The five year survival percentage is :",simcohort.get_survival_percentage(), "; if the annual mortality rate is :", MORTALITY_PROB)



# calculating prediction interval for mean survival time
# create multiple cohorts
multiCohort = Cls.MultiCohort(
    ids=range(NUM_SIM_COHORTS),   # [0, 1, 2 ..., NUM_SIM_COHORTS-1]
    pop_sizes=[REAL_POP_SIZE] * NUM_SIM_COHORTS,  # [REAL_POP_SIZE, REAL_POP_SIZE, ..., REAL_POP_SIZE]
    mortality_probs=[MORTALITY_PROB]*NUM_SIM_COHORTS  # [p, p, ....]
)

# simulate all cohorts
multiCohort.simulate(TIME_STEPS)

print("Multicohort 4's (id=3) 5-y survival percentage is", multiCohort.get_cohort_survival_percentages(3))
print("All multicohorts' 5-y survival percentage are:", multiCohort.get_survival_percentages())



# Problem 2:
# "The number of participants that survived beyond 5 years in a cohort of 𝑁 participants” follows binomial distribution
# Since p represent 5-year survival probability, and each one's survival status is independent from other
# Thus the problem here is actually to get the number of success (survived beyond 5 years) in N times independent experiment
