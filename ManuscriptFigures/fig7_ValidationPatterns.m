clear

projectPath = 'W:\Projects\2018-12 POSTHOCSOURCE Project\analysis_maria\CSPRepo';
datTable = readtable(fullfile(projectPath,'cspAnalysis','SCREEN3_list.xlsx'), 'Basic', 1); %insert the name of the subject spreadsheet with paths
addpath('C:\Users\BNPPC08\Desktop\Maria\matlab\toolboxes\MatlabFns\Colourmaps')
addpath('C:\Users\BNPPC08\Desktop\Maria\matlab\toolboxes\eeglab2021.0')
addpath(fullfile(projectPath,'interpretation'))
load(fullfile(projectPath,'Patterns','screenPatterns_14-Jul-2023.mat')) %insert the name of the CSP pattern file
eeglab

%% interpolate missing channels and reorder channels to a common order
clear allPatternsHigh allPatternsLow

subjects = [1:11]; %select all subjects

% load(fullfile(projectPath,'Patterns','chanlocs.mat')) %load template channel location structure

for idxSub = 1:length(subjects)

    idSub = subjects(idxSub);

    %get the list of channels of the subject
    load(char(datTable.Data(idSub)), 'chanlocs');
    if ismember(datTable.Struct(idSub),'Xclean2')
        load(char(datTable.Data(idSub)), 'rmch');
    elseif ismember(datTable.Struct(idSub),'XAl')
        load(char(datTable.Data(idSub)), 'badCh');
        rmch = badCh;
    end

    %get the patterns of the subject
    VHinv = allV1inv{idSub};
    VLinv = allV2inv{idSub};

    %interpolate real patterns of High condition
    EEG.data = real(VHinv(:,1)); %take real part of the first pattern
    EEG.chanlocs = chanlocs0(~rmch);
    EEG.pnts = 1;
    EEG.trials = 1;
    EEG.nbchan = sum(~rmch);
    EEG = pop_interp( EEG, chanlocs,'spherical'); %interpolate missing channels

    %reorder data by a standard channel order
    [~,reorderingIdx] = ismember(lower({chanlocs.labels}), lower({EEG.chanlocs.labels}));
    EEG.chanlocs = EEG.chanlocs(reorderingIdx);
    EEG.data     = EEG.data(reorderingIdx);

    realPattern = EEG.data;

    %interpolate imaginary patterns of High condition
    EEG.data = imag(VHinv(:,1)); %take imaginary part of the first pattern
    EEG.chanlocs = chanlocs0(~rmch);
    EEG.pnts = 1;
    EEG.trials = 1;
    EEG.nbchan = sum(~rmch);
    EEG = pop_interp( EEG, chanlocs,'spherical'); %interpolate missing channels

    %reorder data by a standard channel order
    [~,reorderingIdx] = ismember(lower({chanlocs.labels}), lower({EEG.chanlocs.labels}));
    EEG.chanlocs = EEG.chanlocs(reorderingIdx);
    EEG.data     = EEG.data(reorderingIdx);

    imagPattern = EEG.data;

    PatternsHigh(:,idxSub) = realPattern + 1i*imagPattern; %combine back into complex values

    %interpolate real patterns of Low condition
    EEG.data = real(VLinv(:,1)); %take real part of the first pattern
    EEG.chanlocs = chanlocs0(~rmch);
    EEG.pnts = 1;
    EEG.trials = 1;
    EEG.nbchan = sum(~rmch);
    EEG = pop_interp( EEG, chanlocs,'spherical'); %interpolate missing channels

    %reorder data by a standard channel order
    [~,reorderingIdx] = ismember(lower({chanlocs.labels}), lower({EEG.chanlocs.labels}));
    EEG.chanlocs = EEG.chanlocs(reorderingIdx);
    EEG.data     = EEG.data(reorderingIdx);

    realPattern = EEG.data;

    %interpolate imaginary patterns of High condition
    EEG.data = imag(VLinv(:,1)); %take imaginary part of the first pattern
    EEG.chanlocs = chanlocs0(~rmch);
    EEG.pnts = 1;
    EEG.trials = 1;
    EEG.nbchan = sum(~rmch);
    EEG = pop_interp( EEG, chanlocs,'spherical'); %interpolate missing channels

    %reorder data by a standard channel order
    [~,reorderingIdx] = ismember(lower({chanlocs.labels}), lower({EEG.chanlocs.labels}));
    EEG.chanlocs = EEG.chanlocs(reorderingIdx);
    EEG.data     = EEG.data(reorderingIdx);

    imagPattern = EEG.data;

    PatternsLow(:,idxSub) = realPattern + 1i*imagPattern; %combine back into complex values

end

%% plot average magnitude pattern
subjectsSign = [1:10]; %select significant subjects
avgPatternHigh = mean(abs(PatternsHigh(:,subjectsSign)),2);
avgPatternLow = mean(abs(PatternsLow(:,subjectsSign)),2);

figure;
tiledlayout(1,2)
ax1 = nexttile;
topoplot(avgPatternHigh, chanlocs,'maplimits','minmax','electrodes','off');
% colorbar

title("High magnitude map (group-average)")
f=gca;
f.FontSize = 14;
f.FontName = 'arial';

ax2 = nexttile;
topoplot(avgPatternLow, chanlocs,'maplimits','minmax','electrodes','off');
title("Low magnitude map (group-average)")
% colorbar
mapMagn = colorcet('L3','reverse',1);
mapPhase = colorcet('L9');
colormap(ax1, mapMagn);
colormap(ax2, mapMagn);

f=gca;
f.FontSize = 14;
f.FontName = 'arial';

%%
exportgraphics(gcf,['C:\Users\BNPPC08\Dropbox\CSP_manuscript\Figures for paper\fig7_Validation','fig7_ValidationPatterns.pdf'],...
    'BackgroundColor','none','ContentType','vector' )



