function [V1inv,V2inv,V1,V2,C] = calculatePatternsFun(subnum,allSubOut,toiIdx, freqIdx)
%Calculate predictive CSP patterns
%   Pattern are calculated on all 400 samples of data. The hyperparameters
%   with the highest average classification accuracy across folds is taken
    
    % regularisation coefficient and number of CSP components with the
    % highest classification accuracy
    [~,ind] = max(mean(mean(allSubOut{1, subnum}.accCVmeans{toiIdx, freqIdx},2),3));
    
    cspOut = allSubOut{1,subnum};
    param = cspOut.param;
    regulIdx = allSubOut{1, subnum}.param.hyperParamList(ind,2);
    nChCSPIdx = allSubOut{1, subnum}.param.hyperParamList(ind,1);
    
    % bandpass filter the raw data
    [XTrain1,XTrain2] = prepDataforCSP(param,cspOut.dataCV,freqIdx,toiIdx,0,1);
    %     [XTest1,XTest2] = prepDataforCSP(param,cspOut.dataVal,freqIdx,toiIdx,0,2);
    
    % calculate CSP filters
    [nCh,~,~] = size(XTrain1);
    %     [covN] = noiseCovariance(param,cspOut.dataCV);
    [covN] = [];
    [C,V1,V2,~,~] = calculateCSP(param,XTrain1,XTrain2,nCh,regulIdx,nChCSPIdx,covN);
    
    % calculate patterns from the filters
    V1inv = inv(ctranspose(V1));
    V2inv = inv(ctranspose(V2));

    %normalise the patterns to unit norm
    for cspIdx = 1:size(V1,2)
        V1inv(:,cspIdx) = V1inv(:,cspIdx)./norm(V1inv(:,cspIdx),2);
        V2inv(:,cspIdx) = V2inv(:,cspIdx)./norm(V2inv(:,cspIdx),2);
    end
    
end

