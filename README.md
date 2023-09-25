Analytic Common Spatial Pattern (aCSP) analysis on TMS-EEG data.
The main analysis codes are in "cspAnalysis" folder. The main script is [here] (cspAnalysis/aCSPmain.m) which calls the other finctions in the cspAnalysis subfolder. The data used as input is published in: 10.5281/zenodo.8370584 (not yet publicly available). Only the data from 12 out of 20 analysed subjects are published here due to data protection agreement.
The codes for the creating of figures from the manuscript are in: ManuscriptFigures . Since the training-test data assignment is random and can lead to different results, the dataset, the results of which are reported in the manuscript, is included in: ManuscriptFigures/fig_PSD_TC_data.mat
The code for calculating spatial patterns of the aCSP components is in Patterns/calculatePatternsFun.m.
The code for temporal and spectral analysis of the aCSP component time courses is in: TimeFreq/timecourse.m
The code for calculating empirical chance level for classification is in Statistics/aCSP_Null.m
The codes for calculating spatial pattern similarity across subjects are in: Statistics/CorrelationChannelPermutation.m and Statistics/CorrelationPatternPermutation.m
The code for extracting prediction accuracy from the output of the aCSPmain.m is in: interpretation/extract_accuracies.m
The code for calculating correlation coefficient between the power of real CSP and analytic CSP components is in: interpretation/realVSanalyticCorrelation.m . Since the training-test data assignment is random and can lead to different results, the dataset, the results of which are reported in the manuscript, is included in: interpretation/realVSanalytic_data.mat
The code for creating condition labels (i.e.e excitability labels) from MEP amplitudes is in: label_meps.m