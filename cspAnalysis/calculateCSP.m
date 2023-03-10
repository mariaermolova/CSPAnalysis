function [C,V1,V2,D1,D2] = calculateCSP(param,XTrain1,XTrain2,nCh,regulIdx,nChCSPIdx,covN)
%calculate CSP filters

nTrain1 = size(XTrain1,3);
nTrain2 = size(XTrain2,3);

%demean trials across time
XTrain1 = XTrain1 - mean(XTrain1,2);
XTrain2 = XTrain2 - mean(XTrain2,2);

%covariance of the training data
cov1Tr = zeros(nCh,nCh);
cov2Tr = zeros(nCh,nCh);
if nTrain1==nTrain2
    for trial = 1:nTrain1
        trial1 = squeeze(XTrain1(:,:,trial));
        trial2 = squeeze(XTrain2(:,:,trial));
        cov1Tr = cov1Tr + (trial1*trial1')/trace(trial1*trial1');
        cov2Tr = cov2Tr +(trial2*trial2')/trace(trial2*trial2');
    end
else
    for trial = 1:nTrain1
        trial1 = squeeze(XTrain1(:,:,trial));
        cov1Tr = cov1Tr + (trial1*trial1')/trace(trial1*trial1');
    end
    for trial = 1:nTrain2
        trial2 = squeeze(XTrain2(:,:,trial));
        cov2Tr = cov2Tr + (trial2*trial2')/trace(trial2*trial2');
    end
end
cov1 = cov1Tr/nTrain1;
cov2 = cov2Tr/nTrain2;
covJoint = (cov1Tr+cov2Tr)/(nTrain1+nTrain2);

%csp
%TODO: decide between mean variance or separate for each channel
% [V1,D,~] = eig(cov1, cov1 + cov2 + mean(diag(covJoint))*(eye(nCh)*param.regul(regulIdx)),'qz');   % Mixing matrix V (spatial filters are columns)
[V1,D,~] = eig(cov1, cov1 + cov2 + diag(covJoint)'.*(eye(nCh)*param.regulCoef(regulIdx)),'qz');   % Mixing matrix V (spatial filters are columns)

%sort eigenvalues and vectors
[D1,d_order] = sort(D(D~=0), 'descend');
V1 = V1(:,d_order);

% [V2,D,~] = eig(cov2, cov1 + cov2 + mean(diag(covJoint))*(eye(nCh)*param.regul(regulIdx)),'qz');   % Mixing matrix V (spatial filters are columns)
[V2,D,~] = eig(cov2, cov1 + cov2 + diag(covJoint)'.*(eye(nCh)*param.regulCoef(regulIdx)),'qz');   % Mixing matrix V (spatial filters are columns)

%sort eigenvalues and vectors
[D2,d_order] = sort(D(D~=0), 'descend');
V2 = V2(:,d_order);

% matrix of selected filters
C = zeros(size(V2,1),param.nChCSP(nChCSPIdx));
for cspIdx = 1:param.nChCSP(nChCSPIdx)/2
    C(:,cspIdx) = V1(:,cspIdx);
    C(:,param.nChCSP(nChCSPIdx)-cspIdx+1) = V2(:,cspIdx);
end

end

