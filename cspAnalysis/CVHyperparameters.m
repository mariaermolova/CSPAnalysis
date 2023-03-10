function [accCVmean,idxnChCSP,idxRegulCoef] = CVHyperparameters(param,data1,data2,covN)
%Select hyperparameters via CV

%create indices for 10 folds
[nCh,~,nTr] = size(data1);

indClass1 = crossvalind('Kfold',nTr,param.nFolds);
indClass2 = crossvalind('Kfold',nTr,param.nFolds);

accCV = zeros(size(param.hyperParamList,1),param.nFolds);

for idxFold = 1:param.nFolds
    
    parfor idxHParam = 1:size(param.hyperParamList,1)
        
        idxnChCSP = param.hyperParamList(idxHParam,1);
        idxRegulCoef = param.hyperParamList(idxHParam,2);
        
        %divide into training and test data
        [XTrain1,XTrain2,XTest1,XTest2] = divideTrainTest(data1,data2,indClass1,indClass2,idxFold);
        
        %calculate csp filters from training data
        [C,~,~] = calculateCSP(param,XTrain1,XTrain2,nCh,idxRegulCoef,idxnChCSP,covN);
        
        %calculate variance from CSP-filtered covariance
        [XTrainBP,YTrain,XTestBP,YTest] = calculateVar(param,XTrain1,XTrain2,XTest1,XTest2,C,nCh,idxnChCSP);
        
        % Classification
        [~,~,acc] = classifyLDA(XTrainBP,XTestBP,YTrain,YTest);
        
        accCV(idxHParam,idxFold)=acc;
    end
    
end
accCVmean = mean(accCV,2);
[~,bestParamIdx] = max(accCVmean);

idxnChCSP = param.hyperParamList(bestParamIdx,1);
idxRegulCoef = param.hyperParamList(bestParamIdx,2);

end

