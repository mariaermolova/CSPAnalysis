%TODO: different times between subjects, automatize cutting
% addpath('C:\Users\BNPPC08\Desktop\Maria\matlab\Projects\CSP')
%dependencies: BioInformatics Toolbox

clear
% clearvars -except sub_table phases

%% load subject info
datTable = readtable('SCREEN3_list.xlsx', 'Basic', 1);

subjects = [1:11];

allSub = cell(1,length(subjects));

rng(0)

for subnum = subjects
    tic

    fprintf(num2str(subnum))

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
    param.freq = [8]; %set all start frequencies [7 13 22 31] [4 8 13 30] [8]
    param.freqband = [22]; %set all sizes of freq bands, same length as freq [6 9 9 10] [4 5 17 10] [22]
    param.nChCSP = [2 4 6]; %set all total numbers of CSP channels
    param.toi = [size(dataCV,2)-494]; %set toi
%     param.toi =  1:250:751; %991:100:1490 491:100:990 691:200:1091 191:200:591 1:250:751
    param.toiWindow = 494; %set total time window for sliding, 499 for single window
    param.regul = [1e-8 1e-6 1e-4 1e-2 1e-1]; % regularization of covariance matrix
    param.nFolds = 5;
    param.nTimes = 5; %1 5
    param.bpfiltparam.FilterOrder = 6;
    param.bpfiltparam.DesignMethod = 'butter';
    param.SampleRate = 1000;
    param.class.CV = classIdCV;
    param.ndel = 500;
    hyperParamList(:,1) = [repmat(1,5,1); repmat(2,5,1); repmat(3,5,1)];
    hyperParamList(:,2) = [repmat([1 2 3 4 5]',3,1)];
    param.hyperParamList = hyperParamList;


    clearvars -except dataCV dataVal param subnum subjects datTable allSub phases


    % Prepare storage for results
    cspOut = struct;
    ACC = cell(length(param.toi),length(param.freq));
    KAPPA = cell(length(param.toi),length(param.freq));
    Model = cell(param.nFolds,param.nTimes,length(param.toi),length(param.freq));
    accCVmeans = cell(length(param.toi),length(param.freq));

    for freqIdx = 1:length(param.freq)
        for toiIdx = 1:length(param.toi)

            CM{toiIdx,freqIdx} = zeros(2);
            ACC{toiIdx,freqIdx}.all = zeros(param.nFolds,param.nTimes);
            accCVmeans{toiIdx,freqIdx} = zeros(size(param.hyperParamList,1),param.nFolds,param.nTimes);

            % prepare data: filter in foi, hilbert transform, cut
            % to toi, equalize trials if needed
            [data1,data2] = prepDataforCSP(param,dataCV,freqIdx,toiIdx,0,1); % trialFlag(0=equal trials,1=non-equal),classFlag(1=CV,2=Val)

            % [covN] = noiseCovariance(param,dataCV);
            [covN] = [];

            [nCh,nTm,nTr1] = size(data1);
            [~,~,nTr2] = size(data2);

            for idxTime = 1:param.nTimes

                %create indices for 10 folds
                indClass1 = crossvalind('Kfold',nTr1,param.nFolds);
                indClass2 = crossvalind('Kfold',nTr2,param.nFolds);

                for idxFold = 1:param.nFolds

                    Model{idxFold,idxTime,toiIdx,freqIdx} = [];

                    %divide into training and test data
                    [XTrain1,XTrain2,XTest1,XTest2] = divideTrainTest(data1,data2,indClass1,indClass2,idxFold);

                    [accCVmean,bestParamIdx] = CVHyperparameters(param,XTrain1,XTrain2,covN);

                    nChCSPIdx = param.hyperParamList(bestParamIdx,1);
                    regulIdx = param.hyperParamList(bestParamIdx,2);

                    %calculate csp filters from training data
                    [C,~,~] = calculateCSP(param,XTrain1,XTrain2,nCh,regulIdx,nChCSPIdx,covN);

                    %calculate variance from CSP-filtered covariance
                    [XTrainBP,YTrain,XTestBP,YTest] = calculateVar(param,XTrain1,XTrain2,XTest1,XTest2,C,nCh,nChCSPIdx);
                    [model,cm,acc] = classifyLDA(XTrainBP,XTestBP,YTrain,YTest);

                    KAPPA{toiIdx,freqIdx}.all(idxFold,idxTime) = kappa(cm);
                    ACC{toiIdx,freqIdx}.all(idxFold,idxTime) = acc;
                    Model{idxFold,idxTime,toiIdx,freqIdx} = model;
                    accCVmeans{toiIdx,freqIdx}(:,idxFold,idxTime) = accCVmean;

                end
            end


            % Mean accuracies
            ACC{toiIdx,freqIdx}.mean = (mean(mean(ACC{toiIdx,freqIdx}.all),1));
            ACC{toiIdx,freqIdx}.std = std(ACC{toiIdx,freqIdx}.all(:));
            ACC{toiIdx,freqIdx}.sem = ACC{toiIdx,freqIdx}.std/sqrt(param.nTimes*param.nFolds);
            KAPPA{toiIdx,freqIdx}.mean = (mean(mean(KAPPA{toiIdx,freqIdx}.all),1));
            KAPPA{toiIdx,freqIdx}.std = std(KAPPA{toiIdx,freqIdx}.all(:));
            KAPPA{toiIdx,freqIdx}.sem = KAPPA{toiIdx,freqIdx}.std/sqrt(param.nTimes*param.nFolds);

            fprintf(num2str(ACC{toiIdx,freqIdx}.mean))
            
        end
    end

    %Save output
    cspOut.ACC = ACC;
    cspOut.KAPPA = KAPPA;
    cspOut.LDA = Model;
    cspOut.param = param;
    cspOut.dataCV = dataCV;
    cspOut.accCVmeans = accCVmeans;
    allSub{1,subnum} = cspOut;


    toc


end
%%
save(append('W:\Projects\2018-12 POSTHOCSOURCE Project\analysis_maria\CSPRepo\screen3\output\screen3_',date),'allSub','-v7.3')


