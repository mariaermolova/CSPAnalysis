%% Calculate spatial similarity between patterns via correlation analysis
% Interpolate missing channels, reorder channels to a template structure,
% calculate spatial correlation between every pair of patterns and average
% across pairs. Derive upper limit of the CI by permuting patterns in every subject
% and repeating the analysis 10 000 times.

clear
projectPath = 'W:\Projects\2018-12 POSTHOCSOURCE Project\analysis_maria\CSPRepo';
addpath(fullfile(projectPath,'cspAnalysis'))
addpath(fullfile(projectPath,'Patterns'))
addpath('C:\Users\BNPPC08\Desktop\Maria\matlab\toolboxes\eeglab14_1_2b')
eeglab
load(fullfile(projectPath,'output','reftep_15-Mar-2023.mat'), 'allSubOut')
load(fullfile(projectPath,'Patterns','chanlocs.mat'))
datTable = readtable(fullfile(projectPath,'CSPAnalysis','REFTEP_list.xlsx'), 'Basic', 1);

%% interpolate missing channels and reorder channels to a common order
clear allPatterns Vinv

subjects = [1:7,9:20]; %select significant subjects

for idxSub = 1:length(subjects)

    idSub = subjects(idxSub);

    %calculate patterns (either High or Low)
    %[V1inv,V2inv,V1,V2,C] = calculatePatternsFun(iSub,allSubOut,toiIdx,freqIdx);
    [Vinv,~,~,~,~] = calculatePatternsFun(idSub,allSubOut,1,1);

    %load channel info
    load(char(datTable.Data(idSub)), 'chanlocs0');
    if ismember(datTable.Struct(idSub),'Xclean2')
        load(char(datTable.Data(idSub)), 'rmch');
    elseif ismember(datTable.Struct(idSub),'XAl')
        load(char(datTable.Data(idSub)), 'badCh');
        rmch = badCh;
    end

    %interpolate missing channels
    EEG.data = abs(Vinv(:,1)); %select first pattern, compute magnitude
    EEG.chanlocs = chanlocs0(~rmch);
    EEG.pnts = 1;
    EEG.trials = 1;
    EEG.nbchan = sum(~rmch);
    EEG = pop_interp( EEG, chanlocs,'spherical');

    %reorder data by a standard channel order
    [~,reorderingIdx] = ismember(lower({chanlocs.labels}), lower({EEG.chanlocs.labels}));
    EEG.chanlocs = EEG.chanlocs(reorderingIdx);
    EEG.data     = EEG.data(reorderingIdx);

    allPatterns(:,idxSub) = EEG.data;

end
%% calculate spatial correlation of each pattern with an average pattern
clear coefSim

%calculate average pattern
avgPat = mean(allPatterns,2);

%calculate correlation coefficient for each pattern with the average pattern
for iPat = 1 : size(allPatterns,2)
    testPat = allPatterns(:,iPat);
    r = corrcoef(testPat,avgPat);
    coefSim(iPat) = r(1,2);
end

%average coefficients across all pairs
corrTrue = mean(coefSim)
%% repeat channel interpolation for all patterns of the subject
clear allPatInterpl

for idxSub = 1:length(subjects)

    idSub = subjects(idxSub);

    %calculate patterns (either High or Low)
    %[V1inv,V2inv,V1,V2,C] = calculatePatternsFun(iSub,allSubOut,toiIdx,freqIdx);
    [Vinv,~,~,~,~] = calculatePatternsFun(idSub,allSubOut,1,1); 

    %load channel info
    load(char(datTable.Data(idSub)), 'chanlocs0');
    if ismember(datTable.Struct(idSub),'Xclean2')
        load(char(datTable.Data(idSub)), 'rmch');
    elseif ismember(datTable.Struct(idSub),'XAl')
        load(char(datTable.Data(idSub)), 'badCh');
        rmch = badCh;
    end

    %interpolate missing channels
    EEG.data = abs(Vinv); %sic! including all patterns now
    EEG.chanlocs = chanlocs0(~rmch);
    EEG.pnts = size(Vinv,2);
    EEG.trials = 1;
    EEG.nbchan = sum(~rmch);
    EEG = pop_interp( EEG, chanlocs,'spherical');

    %reorder data by a standard channel order
    [~,reorderingIdx] = ismember(lower({chanlocs.labels}), lower({EEG.chanlocs.labels}));
    EEG.chanlocs = EEG.chanlocs(reorderingIdx);
    EEG.data     = EEG.data(reorderingIdx,:);

    allPatInterpl{idxSub} = EEG.data;

end
%% calculate upper CI limit from null distribution
clear cosSimNull
for iNull = 1:10000

    clear allPatterns coefSim

    %randomly permute patterns in each subject
    for iSub = 1:length(subjects)
        iPat = randi(size(allPatInterpl{iSub},2),1);
        allPatterns(:,iSub) = allPatInterpl{iSub}(:,iPat);
    end

    %calculate average pattern
    avgPat = mean(allPatterns,2);

    %calculate correlation coefficient for each pattern with the average pattern
    for iPat = 1 : size(allPatterns,2)
        testPat = allPatterns(:,iPat);
        r = corrcoef(testPat,avgPat); 
        coefSim(iPat) = r(1,2);
    end

    %average coefficients across all pairs
    cosSimNull(iNull) = mean(coefSim);

end

upperCI = quantile(cosSimNull,0.95)

