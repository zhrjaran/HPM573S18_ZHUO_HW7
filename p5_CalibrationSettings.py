import scipy.stats as stat


SIM_POP_SIZE = 1000       # population size of simulated cohorts
TIME_STEPS = 1000        # length of simulation
ALPHA = 0.05             # significance level for calculating confidence intervals
NUM_SIM_COHORTS = 500   # number of simulated cohorts used to calculate prediction intervals

# HR: why do  we need multiple cohort? each cohort represent one value of the parameter theta

# details of a clinical study estimating the mean survival time
OBS_N = 1146      # number of patients involved in the study
OBS_SUR= 800        # number of patients survived after 5 years in the study
OBS_ALPHA = 0.05   # significance level
# the standard deviation of the mean survival time reported in the clinical study
# assumes that the reported confidence interval in this study is a t-confidence interval
#OBS_STDEV = OBS_HL / stat.t.ppf(1 - OBS_ALPHA / 2, OBS_N-1)

# HR: OBS_N-1 is the value of degree freedom
# HR: here, we should pay attention to the formula which calculate the standard deviation
# we need 1-alpha/2 instead of alpha/2 to get the t-distribution

# how to sample the posterior distribution of mortality probability
# minimum, maximum and the number of samples for the mortality probability
POST_L, POST_U, POST_N = 0.05, 0.20, 1000
# HR: POST_N here is the number of parameters we want to simulate