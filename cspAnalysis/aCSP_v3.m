%the script to rule them all
%dependencies: Signal Processing Toolbox, Image Processing Toolbox, Statistics and Machine Learning Toolbox, BioInformatics Toolbox

clear

%% Preparation

dataTable = readtable('SCREEN3_list.xlsx', 'Basic', 1); %load data path info

subjects = [1:11]; %select subjects for analysis

allSubOut = cell(1,length(subjects)); %main output

rng(0) %fix the random seed, for reproducibility (doesn't work)

%% CSP analysis

for subnum = subjects

    tic

    fprintf(num2str(subnum))

    %% Load the data

    dataTableSub = dataTable(subnum,:);

    [eeg,~,~,~,mepSorter,labels] = loadData(dataTableSub);

    %% Label data into 2 classes and select sample size

    [eeg,classId] = labelData(eeg,labels,mepSorter,200); %select number of trials in each class: 200

    %% Set analysis parameters

    param = [];
    param.subID = dataTableSub.ID; %dataset id
    param.freq = [8]; %lowest frequencies of analysed freq bands: [8] [7 13 22 31] [4 8 13 30] 
    param.freqWindow = [22]; %widths of analysed freq bands (freq:(freq+freqband)), vector of the same length as param.freq: [22] [6 9 9 10] [4 5 17 10] 
    param.toiWindow = 499; %time window in samples, 499 for reftep analysis
    param.toi = [size(eeg,2)-toiWindow]; %toi: [size(dataCV,2)-toiWindow]. [1:250:751], if several time windows
    param.nFolds = 5; %number of folds
    param.nTimes = 5; %number of times, 1 or 5
    param.bpfiltparam.FilterOrder = 6; %order of the bandpass filter
    param.bpfiltparam.DesignMethod = 'butter';%type of the bandpass filter
    param.SampleRate = 1000; %SF of input eeg data
    param.class.CV = classId; %store class indices
    param.nChCSP = [2 4 6]; %hyperparameter: number of CSP channels
    param.regulCoef = [1e-8 1e-6 1e-4 1e-2 1e-1]; %hyperparameter: regularization coefficient for csp
    hyperParamList(:,1) = [repmat(1,5,1); repmat(2,5,1); repmat(3,5,1)]; %create all combos of hyperparameters for parallel processing
    hyperParamList(:,2) = [repmat([1 2 3 4 5]',3,1)];
    param.hyperParamList = hyperParamList;

    %% Clean up and prepare storage space for results

    clearvars -except dataCV dataVal param subnum subjects datTable allSub

    cspOut = struct; %subject's output
    ACC = cell(length(param.toi),length(param.freq)); % classification accuracy
    KAPPA = cell(length(param.toi),length(param.freq)); % kappa
    Model = cell(param.nFolds,param.nTimes,length(param.toi),length(param.freq)); % LDA model
    accCVmeans = cell(length(param.toi),length(param.freq)); % results of hyperparameter selection

    %% Loop through frequency bands and time windows of analysis, if several

    for idxFreq = 1:length(param.freq)

        for idxToi = 1:length(param.toi)

            %% Prepare for CSP

            % Pre-allocate space
            CM{idxToi,idxFreq} = zeros(2);
            ACC{idxToi,idxFreq}.all = zeros(param.nFolds,param.nTimes);
            accCVmeans{idxToi,idxFreq} = zeros(size(param.hyperParamList,1),param.nFolds,param.nTimes);

            % Prepare data for analysis
            % filter in foi, hilbert transform, cut to toi, equalize trials if needed
            [data1,data2] = prepDataforCSP(param,eeg,idxFreq,idxToi,0,1); % trialFlag(0=equal N of trials,1=non-equal),classFlag(1=CV,2=Val)

            % Prepare spatial noise covariance for CSP regularization
            % [covN] = noiseCovariance(param,dataCV);
            [covN] = [];

            % Get channel, time, and trial dimensions
            [nCh,~,nTr] = size(data1);

            %% Loop through times and folds of CV

            for idxTime = 1:param.nTimes

                % Generate CV indices
                idcClass1 = crossvalind('Kfold',nTr,param.nFolds);
                idcClass2 = crossvalind('Kfold',nTr,param.nFolds);

                for idxFold = 1:param.nFolds

                    Model{idxFold,idxTime,idxToi,idxFreq} = [];

                    % Divide data into training and test sets
                    [XTrain1,XTrain2,XTest1,XTest2] = divideTrainTest(data1,data2,idcClass1,idcClass2,idxFold);

                    % Hyperparameter selection
                    [accCVmean,idxnChCSP,idxRegulCoef] = CVHyperparameters(param,XTrain1,XTrain2,covN);

                    % Calculate csp filters from training data
                    [C,~,~] = calculateCSP(param,XTrain1,XTrain2,nCh,idxnChCSP,idxRegulCoef,covN);

                    % Calculate classification features
                    [XTrainBP,YTrain,XTestBP,YTest] = calculateVar(param,XTrain1,XTrain2,XTest1,XTest2,C,nCh,idxnChCSP);

                    % Perform LDA
                    [model,cm,acc] = classifyLDA(XTrainBP,XTestBP,YTrain,YTest);

                    % Store output
                    KAPPA{idxToi,idxFreq}.all(idxFold,idxTime) = kappa(cm);
                    ACC{idxToi,idxFreq}.all(idxFold,idxTime) = acc;
                    Model{idxFold,idxTime,idxToi,idxFreq} = model;
                    accCVmeans{idxToi,idxFreq}(:,idxFold,idxTime) = accCVmean;

                end
            end

            %% Store classification results

            % Mean accuracies
            ACC{idxToi,idxFreq}.mean = (mean(mean(ACC{idxToi,idxFreq}.all),1));
            ACC{idxToi,idxFreq}.std = std(ACC{idxToi,idxFreq}.all(:));
            ACC{idxToi,idxFreq}.sem = ACC{idxToi,idxFreq}.std/sqrt(param.nTimes*param.nFolds);
            KAPPA{idxToi,idxFreq}.mean = (mean(mean(KAPPA{idxToi,idxFreq}.all),1));
            KAPPA{idxToi,idxFreq}.std = std(KAPPA{idxToi,idxFreq}.all(:));
            KAPPA{idxToi,idxFreq}.sem = KAPPA{idxToi,idxFreq}.std/sqrt(param.nTimes*param.nFolds);

            % get a sneakpeek of the result
            fprintf(num2str(ACC{idxToi,idxFreq}.mean))

        end
    end

    %% Save subject's output
    cspOut.ACC = ACC;
    cspOut.KAPPA = KAPPA;
    cspOut.LDA = Model;
    cspOut.param = param;
    cspOut.dataCV = eeg;
    cspOut.accCVmeans = accCVmeans;

    allSubOut{1,subnum} = cspOut;

    toc

end

%% Save all output
save(append('W:\Projects\2018-12 POSTHOCSOURCE Project\analysis_maria\CSPRepo\screen3\output\screen3_',date),'allSubOut','-v7.3')


