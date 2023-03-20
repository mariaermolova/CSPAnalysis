%% Plot average spatial pattern magnitude
% Load spatial patterns, take magnitude part, interpolate missing channels,
% reorder channels to a common order, optionally normalise by subject, plot the
% average across subjects

clear
projectPath = 'W:\Projects\2018-12 POSTHOCSOURCE Project\analysis_maria\CSPRepo';
addpath('C:\Users\BNPPC08\Desktop\Maria\matlab\toolboxes\eeglab14_1_2b')
addpath('C:\Users\BNPPC08\Desktop\Maria\matlab\toolboxes\MatlabFns\Colourmaps')
eeglab
load(fullfile(projectPath,'Patterns','allPatterns_15-Mar-2023.mat'));
datTable = readtable(fullfile(projectPath,'cspAnalysis','REFTEP_list.xlsx'), 'Basic', 1);

%% interpolate missing channels and reorder channels to a common order
clear allPatternsHigh allPatternsLow

subjects = [1:7,9:20]; %select significant subjects

load(fullpath(projectPath,'Patterns','chanlocs.mat')) %load template channel location structure

for idxSub = 1:length(subjects)
    
    idSub = subjects(idxSub);

    %get the list of channels of the subject
    load(char(datTable.Data(idSub)), 'chanlocs0');
    if ismember(datTable.Struct(idSub),'Xclean2')
        load(char(datTable.Data(idSub)), 'rmch');
    elseif ismember(datTable.Struct(idSub),'XAl')
        load(char(datTable.Data(idSub)), 'badCh');
        rmch = badCh;
    end

    %get the patterns of the subject
    VHinv = allV1inv{idSub};
    VLinv = allV2inv{idSub};

    %prepare patterns of High condition
    EEG.data = abs(VHinv(:,1)); %take magnitude of the first pattern
    EEG.chanlocs = chanlocs0(~rmch);
    EEG.pnts = 1;
    EEG.trials = 1;
    EEG.nbchan = sum(~rmch);
    EEG = pop_interp( EEG, chanlocs,'spherical'); %interpolate missing channels

    %reorder data by a standard channel order
    [~,reorderingIdx] = ismember(lower({chanlocs.labels}), lower({EEG.chanlocs.labels}));
    EEG.chanlocs = EEG.chanlocs(reorderingIdx);
    EEG.data     = EEG.data(reorderingIdx);

    allPatternsHigh(:,idxSub) = EEG.data;

    %prepare patterns of Low condition
    EEG.data = abs(VLinv(:,1));
    EEG.chanlocs = chanlocs0(~rmch);
    EEG.pnts = 1;
    EEG.trials = 1;
    EEG.nbchan = sum(~rmch);
    EEG = pop_interp( EEG, chanlocs,'spherical');

    [~,reorderingIdx] = ismember(lower({chanlocs.labels}), lower({EEG.chanlocs.labels}));
    EEG.chanlocs = EEG.chanlocs(reorderingIdx);
    EEG.data     = EEG.data(reorderingIdx);

    allPatternsLow(:,idxSub) = EEG.data;

end

%% normalise by 1-norm
for idxSub = 1:length(subjects)
    allPatternsHigh(:,idxSub) = allPatternsHigh(:,idxSub)./norm(allPatternsHigh(:,idxSub),1);
    allPatternsLow(:,idxSub) = allPatternsLow(:,idxSub)./norm(allPatternsLow(:,idxSub),1);
end
%% normalise by 2-norm
for idxSub = 1:length(subjects)
    allPatternsHigh(:,idxSub) = allPatternsHigh(:,idxSub)./norm(allPatternsHigh(:,idxSub),2);
    allPatternsLow(:,idxSub) = allPatternsLow(:,idxSub)./norm(allPatternsLow(:,idxSub),2);
end
%% normalise by z-scoring
for idxSub = 1:length(subjects)
    allPatternsHigh(:,idxSub) = (allPatternsHigh(:,idxSub)-mean(allPatternsHigh(:,idxSub)))./norm(allPatternsHigh(:,idxSub),2);
    allPatternsLow(:,idxSub) = (allPatternsLow(:,idxSub)-mean(allPatternsHigh(:,idxSub)))./norm(allPatternsHigh(:,idxSub),2);
end
%% plot average pattern
avgPatternHigh = mean(allPatternsHigh,2);
avgPatternLow = mean(allPatternsLow,2);

figure; tiledlayout(1,2)
nexttile
topoplot(avgPatternHigh, chanlocs,'style','map','maplimits','minmax','electrodes','off');
colorcet('L3','reverse',1)
title("Mean High pattern")

nexttile
topoplot(avgPatternLow, chanlocs,'style','map','maplimits','minmax','electrodes','off');
colorcet('L3','reverse',1)
title("Mean Low pattern")




