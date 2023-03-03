%TODO: different times between subjects, automatize cutting
addpath('W:\Projects\2018-12 POSTHOCSOURCE Project\analysis_maria\CSPRepo\cspAnalysis')
clear
% clearvars -except sub_table phases

%% load subject info
datTable = readtable('subsTue.xlsx', 'Basic', 1);

subjects = [3,9,16,19,20];

allSub = cell(1,length(subjects));

rng(0)

for subnum = subjects
    tic

    fprintf('%d\n', subnum)

    datTableSub = datTable(subnum,:);

    [Xclean2,chanlocs0,badCh,goodTrials,isort,y] = loadData(datTableSub);


    %%

    Xclean2 = Xclean2(:,:,isort([1:200,end-199:end])); % matrix of all data (chan x time x trial)
    classId = double(y(isort([1:200,end-199:end])));
    classId(classId == 0) = 2;

    dataCV = Xclean2;
    classIdCV = classId;

    %%

    %set parameters
    param = [];
    param.subNum = datTableSub.ID; %dataset id
    param.freq = [8]; %set all start frequencies [7 13 22 31]
    param.freqband = [22]; %set all sizes of freq bands, same length as freq [6 9 9 10]
    param.nChCSP = [2 4 6]; %set all total numbers of CSP channels
    param.toi = [size(dataCV,2)-499];; %set toi 991 491
    %     param.toi = [size(dataCV,2)-506]; %set toi 991 491
    %     param.toi = 191:200:591; %991:100:1490 491:100:990 691:200:1091 191:200:591
    param.toiWindow = 499; %set total time window for sliding, 499 for single window
    param.regul = [1e-8 1e-6 1e-4 1e-2 1e-1]; % regularization of covariance matrix
    param.nFolds = 5;
    param.nTimes = 1;
    param.bpfiltparam.FilterOrder = 6;
    param.bpfiltparam.DesignMethod = 'butter';
    param.SampleRate = 1000;
    param.class.CV = classIdCV;
    param.ndel = 500;
    hyperParamList(:,1) = [repmat(1,5,1); repmat(2,5,1); repmat(3,5,1)];
    hyperParamList(:,2) = [repmat([1 2 3 4 5]',3,1)];
    param.hyperParamList = hyperParamList;
    param.nTheory = 5000;

    clearvars -except dataCV dataVal param subnum subjects datTable allSub phases

    cspOut = zeros(param.nTheory,1);

    freqIdx = 1;
    toiIdx = 1;



    clearvars -except theorIdx cspOut freqIdx toiIdx dataCV dataVal param subnum subjects datTable allSub phases

    %tic
    %         fprintf('%d\n',theorIdx)
    %         param.class.CV = param.class.CV(randperm(length(param.class.CV)));

    % prepare data: filter in foi, hilbert transform, cut
    % to toi, equalize trials if needed
    [data1,data2] = prepDataforCSP(param,dataCV,freqIdx,toiIdx,0,1); % trialFlag(0=equal trials,1=non-equal),classFlag(1=CV,2=Val)

    [covN] = [];

    [nCh,nTm,nTr1] = size(data1);
    [~,~,nTr2] = size(data2);

    ACC_all_folds = zeros(param.nFolds,param.nTimes);

    for idxTime = 1%:param.nTimes

        %create indices for 10 folds
        indClass1 = crossvalind('Kfold',nTr1,param.nFolds);
        indClass2 = crossvalind('Kfold',nTr2,param.nFolds);

        for idxFold = 1%:param.nFolds

            %divide into training and test data
            [XTrain1,XTrain2,XTest1,XTest2] = divideTrainTest(data1,data2,indClass1,indClass2,idxFold);

            [accCVmean,bestParamIdx] = CVHyperparameters(param,XTrain1,XTrain2,covN);

            nChCSPIdx = param.hyperParamList(bestParamIdx,1);
            regulIdx = param.hyperParamList(bestParamIdx,2);

            %calculate csp filters from training data
            [C,~,~] = calculateCSP(param,XTrain1,XTrain2,nCh,regulIdx,nChCSPIdx,covN);

            %calculate variance from CSP-filtered covariance
            [XTrainBP,YTrain,XTestBP,YTest] = calculateVar(param,XTrain1,XTrain2,XTest1,XTest2,C,nCh,nChCSPIdx);

            parfor theorIdx = 1:nTheory
                YTestPerm = YTest(randperm(length(YTest)));
                [~,~,acc] = classifyLDA(XTrainBP,XTestBP,YTrain,YTestPerm);
                cspOut(theorIdx) = acc;

            end
        end
    end

    save(append('nullDistribution_sub-',num2str(subnum), '_', date),'cspOut','-v7.3')
    
    fprintf('\n alpha = %d \n', quantile(cspOut,0.95))

    toc
end




