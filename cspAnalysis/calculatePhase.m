function [XTrainPhase,YTrain,XTestPhase,YTest] = calculatePhase(XTrain1,XTrain2, XTest1,XTest2, C)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
nTrain1 = size(XTrain1,3);
nTrain2 = size(XTrain2,3);
nTrain = nTrain1 + nTrain2;

%predictor features
XTrain1TC = pagemtimes(ctranspose(C),XTrain1);%filter High data with both filters
XTrain2TC = pagemtimes(ctranspose(C),XTrain2);%filter Low data with both filters
XTrain1Angle = squeeze(angle(XTrain1TC(:,end,:)));
XTrain2Angle = squeeze(angle(XTrain2TC(:,end,:)));
XTrainPhase = [cos(XTrain1Angle),sin(XTrain1Angle);cos(XTrain2Angle),sin(XTrain2Angle)]';

%class labels
YTrain = ones(nTrain, 1);
YTrain(nTrain1+1:end) = 2;

nTest1 = size(XTest1,3);
nTest2 = size(XTest2,3);
nTest = nTest1 + nTest2;

%predictor features
XTest1TC = pagemtimes(ctranspose(C),XTest1);%filter High data with both filters
XTest2TC = pagemtimes(ctranspose(C),XTest2);%filter Low data with both filters
XTest1Angle = squeeze(angle(XTest1TC(:,end,:)));
XTest2Angle = squeeze(angle(XTest2TC(:,end,:)));
XTestPhase = [cos(XTest1Angle),sin(XTest1Angle);cos(XTest2Angle),sin(XTest2Angle)]';

%class labels
YTest = ones(nTest, 1);
YTest(nTest1+1:end) = 2;