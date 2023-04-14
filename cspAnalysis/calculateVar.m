function [XTrainBP,YTrain,XTestBP,YTest] = calculateVar(param,XTrain1,XTrain2,XTest1,XTest2,C,nCh,nChCSPIdx)
%filter trial-wise covariance matrix with CSP, get variance
%TODO: no  need to separate covariances of classes?

nTrain1 = size(XTrain1,3);
nTrain2 = size(XTrain2,3);
nTrain = nTrain1 + nTrain2;

%covariance of the training data
covTrain1Tr = zeros(nCh,nCh);
covTrain2Tr = zeros(nCh,nCh);
XTrain1BP = zeros(nTrain1,param.nChCSP(nChCSPIdx));
XTrain2BP = zeros(nTrain2,param.nChCSP(nChCSPIdx));
% XTrain1BPlog = zeros(nTrain1,param.nChCSP(nChCSPIdx));
% XTrain2BPlog = zeros(nTrain2,param.nChCSP(nChCSPIdx));
if nTrain1==nTrain2
    for trial = 1:nTrain1
        trial1 = squeeze(XTrain1(:,:,trial));
        trial2 = squeeze(XTrain2(:,:,trial));                                                     
        covTrain1Tr = (trial1*trial1')/trace(trial1*trial1');   % S1~[C x C]
        covTrain2Tr = (trial2*trial2')/trace(trial2*trial2');   % S2~[C x C]
        XTrain1BP(trial,:) = abs(diag(C'*covTrain1Tr*C));
%         XTrain1BPlog(trial,:) = log(XTrain1BP(trial,:)/sum(XTrain1BP(trial,:)));
        XTrain2BP(trial,:) = abs(diag(C'*covTrain2Tr*C));
%         XTrain2BPlog(trial,:) = log(XTrain2BP(trial,:)/sum(XTrain2BP(trial,:)));
                                                                                                                                     
    end
else
    for trial = 1:nTrain1
        trial1 = squeeze(XTrain1(:,:,trial));
        covTrain1Tr(:,:) = (trial1*trial1')/trace(trial1*trial1');   % S1~[C x C]
        XTrain1BP(trial,:) = abs(diag(C'*covTrain1Tr*C));
%         XTrain1BPlog(trial,:) = log(XTrain1BP(trial,:)/sum(XTrain1BP(trial,:)));
    end
    for trial = 1:nTrain2
        trial2 = squeeze(XTrain2(:,:,trial));
        covTrain2Tr(:,:) = (trial2*trial2')/trace(trial2*trial2');   % S2~[C x C]
        XTrain2BP(trial,:) = abs(diag(C'*covTrain2Tr*C));
%         XTrain2BPlog(trial,:) = log(XTrain2BP(trial,:)/sum(XTrain2BP(trial,:)));
    end
end

%predictor features
XTrainBP = [XTrain1BP;XTrain2BP];
% XTrainBP = [XTrain1BPlog;XTrain2BPlog];

%class labels
YTrain = ones(nTrain, 1);
YTrain(nTrain1+1:end) = 2;



nTest1 = size(XTest1,3);
nTest2 = size(XTest2,3);
nTest = nTest1 + nTest2;

%covariance of the test data
covTest1Tr = zeros(nCh,nCh);
covTest2Tr = zeros(nCh,nCh);
XTest1BP = zeros(nTest1,param.nChCSP(nChCSPIdx));
XTest2BP = zeros(nTest2,param.nChCSP(nChCSPIdx));
% XTest1BPlog = zeros(nTest1,param.nChCSP(nChCSPIdx));
% XTest2BPlog = zeros(nTest2,param.nChCSP(nChCSPIdx));
if nTest1==nTest2
    for trial = 1:nTest1
        trial1 = squeeze(XTest1(:,:,trial));
        trial2 = squeeze(XTest2(:,:,trial));
        covTest1Tr(:,:) = (trial1*trial1')/trace(trial1*trial1');   % S1~[C x C]
        covTest2Tr(:,:) = (trial2*trial2')/trace(trial2*trial2');   % S2~[C x C]
        XTest1BP(trial,:) = abs(diag(C'*covTest1Tr*C));
        XTest2BP(trial,:) = abs(diag(C'*covTest2Tr*C));
%         XTest1BPlog(trial,:) = log(XTest1BP(trial,:)/sum(XTest1BP(trial,:)));
%         XTest2BPlog(trial,:) = log(XTest2BP(trial,:)/sum(XTest2BP(trial,:)));
    end
else
    for trial = 1:nTest1
        trial1 = squeeze(XTest1(:,:,trial));
        covTest1Tr(:,:) = (trial1*trial1')/trace(trial1*trial1');   % S1~[C x C]
        XTest1BP(trial,:) = abs(diag(C'*covTest1Tr*C));
%         XTest1BPlog(trial,:) = log(XTest1BP(trial,:)/sum(XTest1BP(trial,:)));
    end
    for trial = 1:nTest2
        trial2 = squeeze(XTest2(:,:,trial));
        covTest2Tr(:,:) = (trial2*trial2')/trace(trial2*trial2');   % S2~[C x C]
        XTest2BP(trial,:) = abs(diag(C'*covTest2Tr*C));
%         XTest2BPlog(trial,:) = log(XTest2BP(trial,:)/sum(XTest2BP(trial,:)));
    end
end

%predictor features
XTestBP = [XTest1BP;XTest2BP];
% XTestBP = [XTest1BPlog;XTest2BPlog];

%class labels
YTest = ones(nTest, 1);
YTest(nTest1+1:end) = 2;

end
