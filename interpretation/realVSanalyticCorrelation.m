%% calculate correlation between real and analytic CSP signal variances
clear
projectPath = 'W:\Projects\2018-12 POSTHOCSOURCE Project\analysis_maria\CSPRepo';
addpath(fullfile(projectPath,'CSPAnalysis'))
allSubOutAna = load(fullfile(projectPath,'output','reftep_15-Mar-2023.mat'));
allSubOutReal = load(fullfile(projectPath,'output','reftep_real_20-Aug-2023.mat'));

%% calculate analytic CSP and real CSP time courses and variances
clear XTestBPAna YTestAna XTestBPReal
subjects = 1:20;

for subnum = subjects

    %prepare data for analysis
    cspOut = allSubOutAna.allSubOut{1,subnum};
    param = cspOut.param;
    data = cspOut.dataCV;

    %bandpass filter and downsample
    [data1,data2] = prepDataforCSP(param,data,1,1,0,1);

    % Get channel, time, and trial dimensions
    [nCh,~,nTr] = size(data1);

    % Generate CV indices
    idcClass1 = crossvalind('Kfold',nTr,param.nFolds);
    idcClass2 = crossvalind('Kfold',nTr,param.nFolds);

    % Divide data into training and test sets
    [XTrain1,XTrain2,XTest1,XTest2] = divideTrainTest(data1,data2,idcClass1,idcClass2,1);

    %select the hyperparameters 
    [~,ind] = max(mean(mean(cspOut.accCVmeans{1, 1},2),3));
    idxRegulCoef = param.hyperParamList(ind,2);
    idxnChCSP = 1;

    % Calculate analytic csp filters from training data
    [C,~,~] = calculateCSP(param,XTrain1,XTrain2,nCh,idxRegulCoef,idxnChCSP,[]);

    [~,~,XTestBP,YTest] = calculateVar(param,XTrain1,XTrain2,XTest1,XTest2,C,nCh,idxnChCSP);

    XTestBPAna(:,:,subnum) = XTestBP;
    YTestAna(:,subnum) = YTest;

    %prepare data for analysis
    cspOut = allSubOutReal.allSubOut{1,subnum};
    param = cspOut.param;
    data = cspOut.dataCV;

    %bandpass filter and downsample
    [data1,data2] = prepDataforCSPreal(param,data,1,1,0,1);

    % Divide data into training and test sets
    [XTrain1,XTrain2,XTest1,XTest2] = divideTrainTest(data1,data2,idcClass1,idcClass2,1);

    %select the hyperparameters 
    [~,ind] = max(mean(mean(cspOut.accCVmeans{1, 1},2),3));
    idxRegulCoef = param.hyperParamList(ind,2);
    idxnChCSP = 1;

    % Calculate csp filters from training data
    [C,~,~] = calculateCSPreal(param,XTrain1,XTrain2,nCh,idxRegulCoef,idxnChCSP,[]);

    [~,~,XTestBP,~] = calculateVar(param,XTrain1,XTrain2,XTest1,XTest2,C,nCh,idxnChCSP);

    XTestBPReal(:,:,subnum) = XTestBP;
end

%% plot real vs analytic CSP variance in each test trial for the High component
figure,tiledlayout('flow')
for subnum = subjects

    %normalise variance by log-transform
    VarAna = log(XTestBPAna(:,1,subnum));
    VarReal = log(XTestBPReal(:,1,subnum));

    %divide into two conditions by labels
    labels = YTestAna(:,subnum);

    nexttile
    scatter(VarAna(labels == 1),VarReal(labels == 1),'red')
    hold on
    scatter(VarAna(labels == 2),VarReal(labels == 2),'blue')
    legend('High trials','Low trials','Location','northwest')
    title('High comp, Sub ',num2str(subnum))
    xlabel('aCSP variance (log)')
    ylabel('CSP variance (log)')

end
%% plot real vs analytic CSP variance in each test trial for the Low component
figure,tiledlayout('flow')
for subnum = subjects

    %normalise variance by log-transform
    VarAna = log(XTestBPAna(:,2,subnum));
    VarReal = log(XTestBPReal(:,2,subnum));

    %divide into two conditions by labels
    labels = YTestAna(:,subnum);

    nexttile
    scatter(VarAna(labels == 1),VarReal(labels == 1),'red')
    hold on
    scatter(VarAna(labels == 2),VarReal(labels == 2),'blue')
    legend('High trials','Low trials','Location','northwest')
    title('Low comp, Sub ',num2str(subnum))
    xlabel('aCSP variance (log)')
    ylabel('CSP variance (log)')
    
end
%% calculate Pearson correlation between CSP and aCSP
clear r_High r_Low p_High p_Low

%normalise variance by log-transform
BPAna = log(XTestBPAna);
BPReal = log(XTestBPReal);

%calculate Pearson corr coefs and p-values 
for subnum = subjects
    [r,p] = corrcoef(BPAna(:,1,subnum),BPReal(:,1,subnum));
    r_High(subnum) = r(1,2);
    p_High(subnum) = p(1,2);
    [r,p] = corrcoef(BPAna(:,2,subnum),BPReal(:,2,subnum));
    r_Low(subnum) = r(1,2);
    p_Low(subnum) = p(1,2);
end

% %% calculate Spearman correlation between CSP and aCSP
% for subnum = subjects
%     [r,p] = corr(BPAna(:,1,subnum),BPReal(:,1,subnum),'type','Spearman');
%     r_High(subnum) = r;
%     p_High(subnum) = p;
%     [r,p] = corr(BPAna(:,2,subnum),BPReal(:,2,subnum),'type','Spearman');
%     r_Low(subnum) = r;
%     p_Low(subnum) = p;
% end



