%% Plot average spatial pattern magnitude
% Load spatial patterns, take magnitude part, interpolate missing channels,
% reorder channels to a common order, normalise by subject, plot the
% average across subjects
clear
addpath('C:\Users\BNPPC08\Desktop\Maria\matlab\toolboxes\eeglab14_1_2b')
addpath 'C:\Users\BNPPC08\Desktop\Maria\matlab\toolboxes\MatlabFns\Colourmaps'
eeglab
load("allPatterns_15-03-2023.mat");
datTable = readtable('W:\Projects\2018-12 POSTHOCSOURCE Project\analysis_maria\CSPRepo\cspAnalysis\REFTEP_list.xlsx', 'Basic', 1);
%% interpolate missing channels and reorder channels to a common order
clear allPatternsHigh allPatternsLow
subjects = [1:7,9:20];
timeIdx = 1;
freqIdx = 1;


for idxSub = 1:length(subjects)
    clear chanlocs0 rmch
    
    iSub = subjects(idxSub);

    load(char(datTable.Data(iSub)), 'chanlocs0');
%     load(char(datTable.Data(iSub)), 'EEG');
    if ismember(datTable.Struct(iSub),'Xclean2')
        load(char(datTable.Data(iSub)), 'rmch');
    elseif ismember(datTable.Struct(iSub),'XAl')
        load(char(datTable.Data(iSub)), 'badCh');
        rmch = badCh;
    end

    V1inv = allV1inv{iSub};
    V2inv = allV2inv{iSub};


    EEG.data = abs(V1inv(:,1));
    EEG.chanlocs = chanlocs0(~rmch);
    EEG.pnts = 1;
    EEG.trials = 1;
    EEG.nbchan = sum(~rmch);
    EEG = pop_interp( EEG, chanlocs,'spherical');

    [~,reorderingIdx] = ismember(lower({chanlocs.labels}), lower({EEG.chanlocs.labels}));
    EEG.chanlocs = EEG.chanlocs(reorderingIdx);
    EEG.data     = EEG.data(reorderingIdx);

    allPatternsHigh(:,idxSub) = EEG.data;

    EEG.data = abs(V2inv(:,1));
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
%% normalise by sum
for idxSub = 1:length(subjects)
    allPatternsHigh(:,idxSub) = allPatternsHigh(:,idxSub)./sum(allPatternsHigh(:,idxSub),'all');
    allPatternsLow(:,idxSub) = allPatternsLow(:,idxSub)./sum(allPatternsLow(:,idxSub),'all');
end
%% normalise by norm
for idxSub = 1:length(subjects)
    allPatternsHigh(:,idxSub) = allPatternsHigh(:,idxSub)./norm(allPatternsHigh(:,idxSub));
    allPatternsLow(:,idxSub) = allPatternsLow(:,idxSub)./norm(allPatternsLow(:,idxSub));
end
%% normalise by z-scoring
for idxSub = 1:length(subjects)
    allPatternsHigh(:,idxSub) = (allPatternsHigh(:,idxSub)-mean(allPatternsHigh(:,idxSub))./norm(allPatternsHigh(:,idxSub)));
    allPatternsLow(:,idxSub) = (allPatternsLow(:,idxSub)-mean(allPatternsHigh(:,idxSub))./norm(allPatternsHigh(:,idxSub)));
end
%% plot average pattern
figure; tiledlayout(1,2)
nexttile
topoplot(mean(allPatternsHigh,2), chanlocs,'style','map','maplimits','minmax','electrodes','off');
colorcet('L3','reverse',1)
title("Avg pattern for High")
%caxis([0 0.07])

nexttile
topoplot(mean(allPatternsLow,2), chanlocs,'style','map','maplimits','minmax','electrodes','off');
colorcet('L3','reverse',1)
title("Avg pattern for Low")
%caxis([0.15 0.51])




