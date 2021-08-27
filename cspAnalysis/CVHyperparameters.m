function [accCVmean,bestParamIdx] = CVHyperparameters(param,data1,data2,covN)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%create indices for 10 folds
[nCh,~,nTr1] = size(data1);
[~,~,nTr2] = size(data2);

indClass1 = crossvalind('Kfold',nTr1,param.nFolds);
indClass2 = crossvalind('Kfold',nTr2,param.nFolds);

accCV = zeros(size(param.hyperParamList,1),param.nFolds);

for idxFold = 1:param.nFolds
    
    parfor idxHParam = 1:size(param.hyperParamList,1)
        
        nChCSPIdx = param.hyperParamList(idxHParam,1);
        regulIdx = param.hyperParamList(idxHParam,2);
        
        %divide into training and test data
        [XTrain1,XTrain2,XTest1,XTest2] = divideTrainTest(data1,data2,indClass1,indClass2,idxFold);
        
        %calculate csp filters from training data
        [C,~,~] = calculateCSP(param,XTrain1,XTrain2,nCh,regulIdx,nChCSPIdx,covN);
        
        %calculate variance from CSP-filtered covariance
        [XTrainBP,YTrain,XTestBP,YTest] = calculateVar(param,XTrain1,XTrain2,XTest1,XTest2,C,nCh,nChCSPIdx);
        
        % Classification
        [~,~,acc] = classifyLDA(XTrainBP,XTestBP,YTrain,YTest);
        
        accCV(idxHParam,idxFold)=acc;
    end
    
end
accCVmean = mean(accCV,2);
[~,bestParamIdx] = max(accCVmean);

end

