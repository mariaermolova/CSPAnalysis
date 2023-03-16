function [V1inv,V2inv,V1,V2,C] = calculatePatternsFun(subnum,allSubOut,toiIdx, freqIdx)
%Calculate predictive CSP patterns
%   Pattern are calculated on all 400 samples of data. The hyperparameters
%   with the highest average classification accuracy across folds is taken
    
    [~,ind] = max(mean(mean(allSubOut{1, subnum}.accCVmeans{toiIdx, freqIdx},2),3));
    
    cspOut = allSubOut{1,subnum};
    param = cspOut.param;
    regulIdx = allSubOut{1, subnum}.param.hyperParamList(ind,2);
    nChCSPIdx = allSubOut{1, subnum}.param.hyperParamList(ind,1);
    
    [XTrain1,XTrain2] = prepDataforCSP(param,cspOut.dataCV,freqIdx,toiIdx,0,1);
    %     [XTest1,XTest2] = prepDataforCSP(param,cspOut.dataVal,freqIdx,toiIdx,0,2);
    
    [nCh,~,~] = size(XTrain1);
    %     [covN] = noiseCovariance(param,cspOut.dataCV);
    [covN] = [];
    [C,V1,V2,~,~] = calculateCSP(param,XTrain1,XTrain2,nCh,regulIdx,nChCSPIdx,covN);
    
    V1inv = inv(ctranspose(V1));
    V2inv = inv(ctranspose(V2));
    
end

