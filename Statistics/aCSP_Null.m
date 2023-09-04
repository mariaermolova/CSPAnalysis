%the script to rule them all
%dependencies: Signal Processing Toolbox, Image Processing Toolbox, Statistics and Machine Learning Toolbox, BioInformatics Toolbox
addpath('W:\Projects\2018-12 POSTHOCSOURCE Project\analysis_maria\CSPRepo\cspAnalysis')
clear

%% Preparation

dataTable = readtable('REFTEP_list.xlsx', 'Basic', 1); %load data path info

subjects = [1:20]; %select subjects for analysis

allSubOut = cell(1,length(subjects)); %main output

rng('shuffle')

hyperParamList(:,1) = [repmat(1,5,1)];%; repmat(2,5,1); repmat(3,5,1)]; %create all combos of hyperparameters for parallel processing
hyperParamList(:,2) = [repmat([1 2 3 4 5]',1,1)];

%% CSP analysis

for subnum = subjects

    tic

    fprintf('Subject %d ',subnum)

    %% Load the data

    dataTableSub = dataTable(subnum,:);

    [eeg,~,~,~,mepsize,labels] = loadData(dataTableSub);
    [eeg,classId] = labelData(eeg,labels,mepsize,200); %select number of trials in each class: 200


    nTheory = 1000;
    cspOut = zeros(nTheory,1);

    parfor theorIdx = 1:nTheory

        %% Label data into 2 classes and select sample size
        classIdPerm = classId(randperm(length(classId)));

        %% Set analysis parameters

        param = [];
        param.subID = dataTableSub.ID; %dataset id
        param.freq = [8]; %lowest frequencies of analysed freq bands: [8] [7 13 22 31] [4 8 13 30]
        param.freqWindow = [22]; %widths of analysed freq bands (freq:(freq+freqband)), vector of the same length as param.freq: [22] [6 9 9 10] [4 5 17 10]
        param.toiWindow = 499; %time window in samples, 499 for reftep analysis
        param.toi = [size(eeg,2)-param.toiWindow]; %toi: [size(dataCV,2)-param.toiWindow]. [1:250:751], if several time windows
        param.nFolds = 5; %number of folds
        param.nTimes = 1; %number of times, 1 or 5
        param.bpfiltparam.FilterOrder = 6; %order of the bandpass filter
        param.bpfiltparam.DesignMethod = 'butter';%type of the bandpass filter
        param.SampleRate = 1000; %SF of input eeg data
        param.class.CV = classIdPerm; %store class indices
        param.nChCSP = [2 4 6]; %hyperparameter: number of CSP channels
        param.regulCoef = [1e-8 1e-6 1e-4 1e-2 1e-1]; %hyperparameter: regularization coefficient for csp
        param.hyperParamList = hyperParamList;
        param.nTheory = nTheory;

        %% Clean up and prepare storage space for results

        %         clearvars -except eeg param subnum subjects dataTable allSubOut theorIdx
        idxFreq = 1;
        idxToi = 1;

        %% Prepare for CSP

        % Prepare data for analysis
        % filter in foi, hilbert transform, cut to toi, equalize trials if needed
        [data1,data2] = prepDataforCSP(param,eeg,idxFreq,idxToi,0,1); % trialFlag(0=equal N of trials,1=non-equal),classFlag(1=CV,2=Val)

        % Prepare spatial noise covariance for CSP regularization
        [covN] = [];

        % Get channel, time, and trial dimensions
        [nCh,~,nTr] = size(data1);

        %% Loop through times and folds of CV (only 1 time and 1 fold here)

        % Generate CV indices
        idcClass1 = crossvalind('Kfold',nTr,param.nFolds);
        idcClass2 = crossvalind('Kfold',nTr,param.nFolds);

        for idxFold = 1

            % Divide data into training and test sets
            [XTrain1,XTrain2,XTest1,XTest2] = divideTrainTest(data1,data2,idcClass1,idcClass2,idxFold);

            % Hyperparameter selection
            [accCVmean,idxnChCSP,idxRegulCoef] = CVHyperparameters(param,XTrain1,XTrain2,covN);

            % Calculate csp filters from training data
            [C,~,~] = calculateCSP(param,XTrain1,XTrain2,nCh,idxRegulCoef,idxnChCSP,covN);

            % Calculate classification features
            [XTrainBP,YTrain,XTestBP,YTest] = calculateVar(param,XTrain1,XTrain2,XTest1,XTest2,C,nCh,idxnChCSP);

            [~,~,acc] = classifyLDA(XTrainBP,XTestBP,YTrain,YTest);
            cspOut(theorIdx) = acc;
        end
    end

    save(append('nullDistribution_reftep_sub', '_',num2str(subnum),'_', date),'cspOut','-v7.3')

    fprintf('\n alpha = %d \n', quantile(cspOut,0.95))

    toc
end


