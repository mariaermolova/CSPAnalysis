%% Plot average spatial pattern phase lag distribution
% Load spatial patterns (complex), interpolate missing channels,
% reorder channels to a common order, re-reference, compute the angle,
% plot the average across subjects

clear
projectPath = 'W:\Projects\2018-12 POSTHOCSOURCE Project\analysis_maria\CSPRepo';
addpath('C:\Users\BNPPC08\Desktop\Maria\matlab\toolboxes\eeglab14_1_2b')
addpath('C:\Users\BNPPC08\Desktop\Maria\matlab\toolboxes\MatlabFns\Colourmaps')
eeglab
load(fullfile(projectPath,'Patterns','allPatterns_15-Mar-2023.mat'));
datTable = readtable(fullfile(projectPath,'cspAnalysis','REFTEP_list.xlsx'), 'Basic', 1);

%% Interpolate missing channels (on complex patterns)
clear allPatternsLow allPatternsHigh

subjects = [1,2,4,9,18]; %select subjects to plot

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
    EEG.data = VHinv(:,1);
    EEG.chanlocs = chanlocs0(~rmch);
    EEG.pnts = 1;
    EEG.trials = 1;
    EEG.nbchan = sum(~rmch);
    EEG = pop_interp( EEG, chanlocs,'spherical');

    %reorder data by a standard channel order
    [~,reorderingIdx] = ismember(lower({chanlocs.labels}), lower({EEG.chanlocs.labels}));
    EEG.chanlocs = EEG.chanlocs(reorderingIdx);
    EEG.data     = EEG.data(reorderingIdx);

    allPatternsHigh(:,idxSub) = EEG.data;

    %prepare patterns of Low condition
    EEG.data = VLinv(:,1);
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
%% Re-reference phases

clear topoAnglesHigh topoAnglesLow

for idxSub = 1:size(allPatternsHigh,2)

    %find the index of a new reference channel
    refind = find(ismember({chanlocs.labels},'C3'));

    %get the phase pattern of High condition
    subPattern = allPatternsHigh(:,idxSub);
    compDiff = conj(subPattern).*subPattern(refind); %re-reference
    topoAnglesHigh(:,idxSub) = angle(compDiff)*180/pi;   %get the phase

    %get the phase pattern of Low condition
    subPattern = allPatternsLow(:,idxSub);
    compDiff = conj(subPattern).*subPattern(refind);
    topoAnglesLow(:,idxSub) = angle(compDiff)*180/pi;
end


%% Plot average phase patterns

avgTopoAnglesHigh = mean(topoAnglesHigh,2);
avgTopoAnglesLow = mean(topoAnglesLow,2);

figure; tiledlayout(1,2)
nexttile
topoplot(avgTopoAnglesHigh, chanlocs,'maplimits','minmax','electrodes','off','emarker2',{[refind],'o','k',6});
colormap(parula)
title("Avg phase pattern for High")
caxis([-180 180])
colorbar

nexttile
topoplot(avgTopoAnglesLow, chanlocs,'maplimits','minmax','electrodes','off','emarker2',{[refind],'o','k',6});
colormap(parula)
title("Avg phase pattern for Low")
caxis([-180 180])
colorbar
