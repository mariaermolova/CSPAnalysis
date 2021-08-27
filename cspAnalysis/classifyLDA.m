function [model,cm,acc] = classifyLDA(XTrainBP,XTestBP,YTrain,YTest)
%perform classification with LDA

% model=fitcdiscr(XTrainBP,YTrain);
model=fitcdiscr(XTrainBP,YTrain,'OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',struct('Verbose',0,'ShowPlots',0));
% model=fitcdiscr(XTrainBP,YTrain,'CrossVal','on');
YPredicted = predict(model,XTestBP);
% YPredicted = kfoldPredict(model,XTestBP);
cm = confusionmat(YTest,YPredicted);                 
acc = sum(YPredicted==YTest)/length(YTest);

end

