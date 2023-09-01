function [XTrainPred,XTestPred] = combineFeatures(XTrainPhase,XTrainBP,XTestPhase,XTestBP)
%Combine Variance and Phase into a feature matrix for LDA

XTrainPred = [XTrainPhase, XTrainBP];
XTestPred = [XTestPhase, XTestBP];
end