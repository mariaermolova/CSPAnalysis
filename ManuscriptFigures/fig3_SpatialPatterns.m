clear

projectPath = 'W:\Projects\2018-12 POSTHOCSOURCE Project\analysis_maria\CSPRepo';
datTable = readtable(fullfile(projectPath,'cspAnalysis','REFTEP_list.xlsx'), 'Basic', 1); %insert the name of the subject spreadsheet with paths
addpath('C:\Users\BNPPC08\Desktop\Maria\matlab\toolboxes\MatlabFns\Colourmaps')
addpath('C:\Users\BNPPC08\Desktop\Maria\matlab\toolboxes\eeglab2021.0\')
addpath(fullfile(projectPath,'interpretation'))
load(fullfile(projectPath,'Patterns','allPatterns_13-Jul-2023.mat')) %insert the name of the CSP pattern file
eeglab
%% interpolate missing channels and reorder channels to a common order
 

subjects = [1:20]; %select significant subjects

load(fullfile(projectPath,'Patterns','chanlocs.mat')) %load template channel location structure
%% interpolate complex topogarphies
clear PatternsHigh PatternsLow

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

    PatternsHigh(:,idxSub) = realPattern + 1i*imagPattern; %combine them back into complex values

    %interpolate real patterns of Low condition
    EEG.data = real(VLinv(:,1)); %take real of the first pattern
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

    %interpolate imaginary patterns of Low condition
    EEG.data = imag(VLinv(:,1)); %take imaginary of the first pattern
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

    PatternsLow(:,idxSub) = realPattern + 1i*imagPattern; %combine them back into complex values

end

%% Re-reference phases

clear topoAnglesHigh topoAnglesLow

for idxSub = 1:size(PatternsHigh,2)

    %find the index of a new reference channel
    refIndHigh = find(ismember({chanlocs.labels},'FCC3h'));

    %get the phase pattern of High condition
    subPattern = PatternsHigh(:,idxSub);
    compDiff = conj(subPattern).*subPattern(refIndHigh); %re-reference
    topoAnglesHigh(:,idxSub) = angle(compDiff)*180/pi;   %get the phase
end

for idxSub = 1:size(PatternsLow,2)

    %find the index of a new reference channel
    refIndLow = find(ismember({chanlocs.labels},'Cz'));

    %get the phase pattern of Low condition
    subPattern = PatternsLow(:,idxSub);
    compDiff = conj(subPattern).*subPattern(refIndLow); %re-reference
    topoAnglesLow(:,idxSub) = angle(compDiff)*180/pi; %get the phase
end

%% plot average magnitude pattern
subjectsSign = [1:7,9:20]; %select significant subjects
avgPatternHigh = mean(abs(PatternsHigh(:,subjectsSign)),2); 
avgPatternLow = mean(abs(PatternsLow(:,subjectsSign)),2);

figure; 
tiledlayout('flow')

ax1 = nexttile;
topoplot(avgPatternHigh, chanlocs,'maplimits','minmax','electrodes','off');
% colorbar
title("A. High magnitude map (group-average)")
f=gca;
f.FontSize = 14;
f.FontName = 'arial';

ax2 = nexttile;
topoplot(avgPatternLow, chanlocs,'maplimits','minmax','electrodes','off');
title("B. Low magnitude map (group-average)")
% colorbar
f=gca;
f.FontSize = 14;
f.FontName = 'arial';
%% plot subset average magnitude pattern 
subjectsHigh = [1,2,4,9,18,19]; %select sugjects with the same High pattern
subjectsLow = [1,12,16,17,18]; %select sugjects with the same Low pattern
avgSubsetPatternHigh = mean(abs(PatternsHigh(:,subjectsHigh)),2); 
avgSubsetPatternLow = mean(abs(PatternsLow(:,subjectsLow)),2); 

ax3 = nexttile;
topoplot(avgSubsetPatternHigh, chanlocs,'maplimits','minmax','electrodes','off');
title("C. High magnitude map (subset-average)")
% colorbar
f=gca;
f.FontSize = 14;
f.FontName = 'arial';

ax4 = nexttile;
topoplot(avgSubsetPatternLow, chanlocs,'maplimits','minmax','electrodes','off');
title("D. Low magnitude map (subset-average)")
% colorbar
f=gca;
f.FontSize = 14;
f.FontName = 'arial';


%% Plot subset average phase pattern

avgTopoAnglesHigh = mean(topoAnglesHigh(:,subjectsHigh),2); 
avgTopoAnglesLow = mean(topoAnglesLow(:,subjectsLow),2); 

%masking by median of average magnitude
maskHigh = avgSubsetPatternHigh>=quantile(avgSubsetPatternHigh,0.50);
maskLow = avgSubsetPatternLow>=quantile(avgSubsetPatternLow,0.50);

ax5 = nexttile;
topoplot_with_method(avgTopoAnglesHigh, chanlocs,'maplimits','minmax','electrodes','off','emarker2',{[refIndHigh],'o','k',6},...
    'pmask',maskHigh,'method','nearest','style','both');
% colorbar
title("E. High phase map (subset-average)")
caxis([-180 180])
f=gca;
f.FontSize = 14;
f.FontName = 'arial';

ax6 = nexttile;
topoplot_with_method(avgTopoAnglesLow, chanlocs,'maplimits','minmax','electrodes','off','emarker2',{[refIndLow],'o','k',6},...
    'pmask',maskLow,'method','nearest','style','both');
% colorbar
title("F. Low phase map (subset-average)")
caxis([-180 180])
f=gca;
f.FontSize = 14;
f.FontName = 'arial';
%%

mapMagn = colorcet('L3','reverse',1);
mapPhase = colorcet('C3','shift',0.5);
colormap(ax1, mapMagn);
colormap(ax2, mapMagn);
colormap(ax3, mapMagn);
colormap(ax4, mapMagn);
colormap(ax5, mapPhase)
colormap(ax6, mapPhase)

%%
exportgraphics(gcf,['C:\Users\BNPPC08\Dropbox\CSP_manuscript\Figures for paper\fig3_SpatialPatterns\','fig3_SpatialPatterns.pdf'],...
    'BackgroundColor','none','ContentType','vector' )

