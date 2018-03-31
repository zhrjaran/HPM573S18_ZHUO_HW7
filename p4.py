import p4_CalibrateClasses as Cls
import p4_CalibrationSettings as CalibSets
import scr.FigureSupport as Fig

# create a calibration object
calibration = Cls.Calibration()

# sample the posterior of the mortality probability
calibration.sample_posterior()

# create the histogram of the resampled mortality probabilities
Fig.graph_histogram(
    data=calibration.get_mortality_resamples(),
    title='Histogram of Resampled Mortality Probabilities',
    x_label='Mortality Probability',
    y_label='Counts',
    x_range=[CalibSets.POST_L, CalibSets.POST_U])

# Estimate of mortality probability and the posterior interval
print('Estimate of mortality probability ({:.{prec}%} credible interval):'.format(1-CalibSets.ALPHA, prec=0),
      calibration.get_mortality_estimate_credible_interval(CalibSets.ALPHA, 4))

# effective sample size
txtEff = 'Effective sample size: {:.1f}'.format(calibration.get_effective_sample_size())
print(txtEff)