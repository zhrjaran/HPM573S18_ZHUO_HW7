from enum import Enum
import scipy.stats as stat
import numpy as np
import scr.InOutFunctions as InOutSupport
import scr.StatisticalClasses as StatSupport
import scr.FormatFunctions as FormatSupport
import p1_SurvivalModel as SurvivalCls
import p4_CalibrationSettings as CalibSets


class CalibrationColIndex(Enum):
    """ indices of columns in the calibration results cvs file  """
    ID = 0          # cohort ID
    W = 1  # likelihood weight
    MORT_PROB = 2   # mortality probability


class Calibration:
    def __init__(self):
        """ initializes the calibration object"""
        np.random.seed(1)   # specifying the seed of the numpy random number generator
        self._cohortIDs = range(CalibSets.POST_N)   # IDs of cohorts to simulate
        self._mortalitySamples = []      # values of mortality probability at which the posterior should be sampled
        self._mortalityResamples = []    # resampled values for constructing posterior estimate and interval
        self._weights = []               # likelihood weights of sampled mortality probabilities
        self._normalizedWeights = []     # normalized likelihood weights (sums to 1)
        self._csvRows = \
            [['Cohort ID', 'Likelihood Weights' ,'Mortality Prob']]  # list containing the calibration results

    def sample_posterior(self):
        """ sample the posterior distribution of the mortality probability """

        # find values of mortality probability at which the posterior should be evaluated
        # HR: this is a prior p, since we just simulate the parameter without having any observations
        self._mortalitySamples = np.random.uniform(
            low=CalibSets.POST_L,
            high=CalibSets.POST_U,
            size=CalibSets.POST_N)

        # create a multi cohort
        multiCohort = SurvivalCls.MultiCohort(
            ids=self._cohortIDs,
            mortality_probs=self._mortalitySamples,
            pop_sizes=[CalibSets.SIM_POP_SIZE]*CalibSets.POST_N
        )

        # simulate the multi cohort
        multiCohort.simulate(CalibSets.TIME_STEPS)

        # calculate the likelihood of each simulated cohort
        for cohort_id in self._cohortIDs:

            # get the survival percentage for this cohort
            survivalrate= multiCohort.get_cohort_survival_percentages(cohort_id)

            # construct a gaussian likelihood
            # with mean calculated from the simulated data and standard deviation from the clinical study.
            # evaluate this pdf (probability density function) at the mean reported in the clinical study.
            weight = stat.binom.pmf(
                k=CalibSets.OBS_SUR,
                p=survivalrate*0.01,
                n=CalibSets.OBS_N
            )

            # store the weight
            self._weights.append(weight)

        # normalize the likelihood weights
        sum_weights = np.sum(self._weights)
        self._normalizedWeights = np.divide(self._weights, sum_weights)

        # re-sample mortality probability (with replacement) according to likelihood weights
        self._mortalityResamples = np.random.choice(
            a=self._mortalitySamples,
            size=CalibSets.NUM_SIM_COHORTS,
            # HR: size is the number of resampled cohorts we would get eventually, this is not the number of cohorts we simulate first time
            replace=True,
            # HR: values with higher weights get sampled more often
            p=self._normalizedWeights)

        # produce the list to report the results
        for i in range(0, len(self._mortalitySamples)):
            self._csvRows.append(
                [self._cohortIDs[i], self._normalizedWeights[i], self._mortalitySamples[i]])

        # write the calibration result into a csv file
        InOutSupport.write_csv('CalibrationResults.csv', self._csvRows)

    def get_mortality_resamples(self):
        """
        :return: mortality resamples
        """
        return self._mortalityResamples

    def get_mortality_estimate_credible_interval(self, alpha, deci):
        """
        :param alpha: the significance level
        :param deci: decimal places
        :returns text in the form of 'mean (lower, upper)' of the posterior distribution"""

        # calculate the credible interval
        sum_stat = StatSupport.SummaryStat('Posterior samples', self._mortalityResamples)

        estimate = sum_stat.get_mean()  # estimated mortality probability
        credible_interval = sum_stat.get_PI(alpha)  # credible interval

        return FormatSupport.format_estimate_interval(estimate, credible_interval, deci)

    def get_effective_sample_size(self):
        """
        :returns: the effective sample size
        """
        return 1 / np.sum(self._normalizedWeights ** 2)
        # HR: formula here can check the wikipedia for "weighted effective sample size"
        # HR: https://en.wikipedia.org/wiki/Effective_sample_size
        # HR: In statistics, effective sample size is a notion defined for
        # a sample from a distribution when the observations in the sample are correlated or weighted

