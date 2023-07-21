%% calculate CSP time courses, their PSD and envelope
clear
addpath('W:\Projects\2018-12 POSTHOCSOURCE Project\toolboxes\bbci')
startup_bbci_toolbox
projectPath = 'W:\Projects\2018-12 POSTHOCSOURCE Project\analysis_maria\CSPRepo';
addpath(fullfile(projectPath,'CSPAnalysis'))
load(fullfile(projectPath,'output','reftep_15-Mar-2023.mat'));

%% fcalculate CSP time courses
clear AllTimeCrsHDat AllTimeCrsLDat HPatHDatSpctr LPatHDatSpctr HPatLDatSpctr LPatLDatSpctr
subjects = 1:20;
for subnum = subjects

    %prepare data for analysis
    cspOut = allSubOut{1,subnum};
    param = cspOut.param;
    data = cspOut.dataCV;
    param.toi = 251; 
    param.toiWindow = 1244; %select -1250:-7 msec before the pulse

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
    [~,ind] = max(mean(mean(allSubOut{1, subnum}.accCVmeans{1, 1},2),3));
    idxRegulCoef = param.hyperParamList(ind,2);
    idxnChCSP = 1;

    %select 0.5 sec of signal before the pulse
    XTrain1 = XTrain1(:,(311-125+1):311,:);
    XTrain2 = XTrain2(:,(311-125+1):311,:);

    % Calculate csp filters from training data
    [C,~,~] = calculateCSP(param,XTrain1,XTrain2,nCh,idxRegulCoef,idxnChCSP,[]);

    % Apply the csp filters to test data
    AllTimeCrsHDat{subnum} = pagemtimes(ctranspose(C),XTest1);%filter High data with both filters
    AllTimeCrsLDat{subnum} = pagemtimes(ctranspose(C),XTest2);%filter Low data with both filters
end

%% calculate Power Spectrum
for subnum = subjects

    %set analysis parameters
    filteredData = [];
    filteredData.fs = 250;
    opt.Scaling = 'db';
    
    %power spectrum of High component on High data
    filteredData.x = permute(AllTimeCrsHDat{subnum}(1,:,:),[2 1 3]);
    filteredData_TEMP = proc_spectrum(filteredData, [8 30], opt);
    F = filteredData_TEMP.t;
    spectrum = filteredData_TEMP.x;
    HPatHDatSpctr(:,subnum) = mean(spectrum,3);

    %power spectrum of Low component on High data
    filteredData.x = permute(AllTimeCrsHDat{subnum}(end,:,:),[2 1 3]);
    filteredData_TEMP = proc_spectrum(filteredData, [8 30], opt);
    F = filteredData_TEMP.t;
    spectrum = filteredData_TEMP.x;
    LPatHDatSpctr(:,subnum) = mean(spectrum,3);

    %power spectrum of High component on Low data
    filteredData.x = permute(AllTimeCrsLDat{subnum}(1,:,:),[2 1 3]);
    filteredData_TEMP = proc_spectrum(filteredData, [8 30], opt);
    F = filteredData_TEMP.t;
    spectrum = filteredData_TEMP.x;
    HPatLDatSpctr(:,subnum) = mean(spectrum,3);

    %power spectrum of Low component on Low data
    filteredData.x = permute(AllTimeCrsLDat{subnum}(end,:,:),[2 1 3]);
    filteredData_TEMP = proc_spectrum(filteredData, [8 30], opt);
    F = filteredData_TEMP.t;
    spectrum = filteredData_TEMP.x;
    LPatLDatSpctr(:,subnum) = mean(spectrum,3);

end


%% Calculate envelope
clear normTimeCrsHPatHDat normTimeCrsHPatLDat normTimeCrsLPatHDat normTimeCrsLPatLDat

[~,nTm,~] = size(AllTimeCrsHDat{subnum});

for subnum = subjects

    TimeCrsHPatHDat = mean(abs(squeeze(AllTimeCrsHDat{subnum}(1,:,:))),2);
    TimeCrsHPatLDat = mean(abs(squeeze(AllTimeCrsLDat{subnum}(1,:,:))),2);
    
    normHPat = [TimeCrsHPatHDat;TimeCrsHPatLDat];
    normHPat = normHPat/norm(normHPat,2);

    normTimeCrsHPatHDat(:,subnum) = normHPat(1:nTm);
    normTimeCrsHPatLDat(:,subnum) = normHPat(nTm+1:end);
   
    TimeCrsLPatHDat = mean(abs(squeeze(AllTimeCrsHDat{subnum}(end,:,:))),2);
    TimeCrsLPatLDat = mean(abs(squeeze(AllTimeCrsLDat{subnum}(end,:,:))),2);

    normLPat = [TimeCrsLPatHDat;TimeCrsLPatLDat];
    normLPat = normLPat/norm(normLPat,2);

    normTimeCrsLPatHDat(:,subnum) = normLPat(1:nTm);
    normTimeCrsLPatLDat(:,subnum) = normLPat(nTm+1:end);

end


