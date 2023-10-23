# Analytic Common Spatial Pattern (aCSP) analysis on TMS-EEG data.
The main analysis codes are [here](cspAnalysis). The main script is [aCSPmain.m](cspAnalysis/aCSPmain.m), it calls the other finctions in the [cspAnalysis](cspAnalysis) subfolder. The data used as input is published in: 10.5281/zenodo.10034964 (not yet publicly available). Only the data from 12 out of 20 analysed subjects are published due to data protection agreement.

The codes for the creating of figures from the manuscript are [here](ManuscriptFigures). Since the training-test data assignment is random and can lead to different results, the dataset, the results of which are reported in the manuscript, is shared [here](ManuscriptFigures/fig_PSD_TC_data.mat).

The code for calculating spatial patterns of the aCSP components is [here](Patterns/calculatePatternsFun.m).

The code for temporal and spectral analysis of the aCSP component time courses is [here](TimeFreq/timecourse.m).

The code for calculating empirical chance level for classification is [here](Statistics/aCSP_Null.m).

The codes for calculating spatial pattern similarity across subjects are [here](Statistics/CorrelationChannelPermutation.m) and [here](Statistics/CorrelationPatternPermutation.m).

The code for extracting prediction accuracy from the output of the `aCSPmain.m` is [here](interpretation/extract_accuracies.m).

The code for calculating correlation coefficient between the power of real CSP and analytic CSP components is [here](interpretation/realVSanalyticCorrelation.m). Since the training-test data assignment is random and can lead to different results, the dataset, the results of which are reported in the manuscript, is shared [here](interpretation/realVSanalytic_data.mat).

The code for creating condition labels (i.e. excitability labels) from MEP amplitudes is [here](label_meps.m).
