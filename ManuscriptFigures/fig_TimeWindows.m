%% Plot average patterns from the results of a post-hoc time-varying analysis
clear

addpath 'W:\Projects\2018-12 POSTHOCSOURCE Project\analysis_maria\CSPRepo\interpretation'
projectPath = 'W:\Projects\2018-12 POSTHOCSOURCE Project\analysis_maria\CSPRepo';
datTable = readtable(fullfile(projectPath,'cspAnalysis','REFTEP_list.xlsx'), 'Basic', 1); %insert the name of the subject spreadsheet with paths
addpath('C:\Users\BNPPC08\Desktop\Maria\matlab\toolboxes\MatlabFns\Colourmaps')
addpath('C:\Users\BNPPC08\Desktop\Maria\matlab\toolboxes\eeglab2021.0')
addpath(fullfile(projectPath,'interpretation'))
load(fullfile(projectPath,'Patterns','timePatterns_14-Jul-2023.mat')) %insert the name of the CSP pattern file
eeglab
%% load classification accuracy
[accs] = extract_accuracies('W:\Projects\2018-12 POSTHOCSOURCE Project\analysis_maria\CSPRepo\output\timeWindows_reftep_17-Mar-2023.mat');
[accsEnd] = extract_accuracies('W:\Projects\2018-12 POSTHOCSOURCE Project\analysis_maria\CSPRepo\output\reftep_15-Mar-2023.mat');
accs(:,5) = accsEnd;
accs = accs*100;
meansSign = accs(:,2:5);
%%
figure;
tiledlayout(3,4,"TileSpacing","compact")
%% plot average time/accuracy with error bars
y = mean(meansSign,1); % your mean vector;
erbars = std(meansSign,[],1)/sqrt(size(meansSign,1));
curve1 = y + erbars;
curve2 = y - erbars;
x2 = [1:4, fliplr(1:4)];
inBetween = [curve1, fliplr(curve2)];
nexttile(1,[1,4])
fill(x2, inBetween, 'r','FaceAlpha',0.1,'EdgeColor','w');
hold on;
plot(1:4, y, 'k', 'LineWidth', 2);
f=gca;
f.XTick = 1:4;
f.XTickLabel = {'-1250:-750 ms','-1000:-500 ms','-750:-250 ms','-500:0 ms'};
f.YTick = [64,66,68,70];
f.YTickLabel = {'64%','66%','68%','70%'};
f.FontSize = 14;
f.FontName = 'arial';

xlim([0.6,4.4])
xlabel('Time window')
ylabel('Accuracy')
ylim([63,70])

%% Interpolate missing channels, plot average patterns for each time window or frequency band
clear PatternsHigh PatternsLow

subs = [1:20]; %select all subjects

load(fullfile(projectPath,'Patterns','chanlocs.mat')) %load template channel location structure

for timeIdx = 2:4 %loop either through time windows
    freqIdx = 1; %or through frequency bands

    allV1inv = allPatterns{timeIdx,freqIdx}.allV1inv;
    allV2inv = allPatterns{timeIdx,freqIdx}.allV2inv;

    for idxSub = 1:length(subs)

        iSub = subs(idxSub);

        %get the list of channels of the subject
        load(char(datTable.Data(iSub)), 'chanlocs0');
        if ismember(datTable.Struct(iSub),'Xclean2')
            load(char(datTable.Data(iSub)), 'rmch');
        elseif ismember(datTable.Struct(iSub),'XAl')
            load(char(datTable.Data(iSub)), 'badCh');
            rmch = badCh;
        end

        %get the patterns of the subject
        VHinv = allV1inv{iSub};
        VLinv = allV2inv{iSub};

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

    %plot average magnitude patterns
    avgPatternHigh = mean(abs(PatternsHigh),2);
    avgPatternLow = mean(abs(PatternsLow),2);

    nexttile(timeIdx+3,[1,1])
    topoplot(avgPatternHigh, chanlocs,'maplimits','minmax','electrodes','off');
    colorcet('L3','reverse',1)
    title("High",'FontSize',14)
    nexttile(timeIdx+7,[1,1])
    topoplot(avgPatternLow, chanlocs,'maplimits','minmax','electrodes','off');
    colorcet('L3','reverse',1)
    title("Low",'FontSize',14)

end

%% get patterns from the main analysis for the last time window

load(fullfile(projectPath,'Patterns','allPatterns_13-Jul-2023.mat')) %insert the name of the CSP pattern file

%% Interpolate missing channels, plot average patterns

clear PatternsHigh PatternsLow

subs = [1:20]; %select all subjects
timeIdx = 5;
load(fullfile(projectPath,'Patterns','chanlocs.mat')) %load template channel location structure

for idxSub = 1:length(subs)

    iSub = subs(idxSub);

    %get the list of channels of the subject
    load(char(datTable.Data(iSub)), 'chanlocs0');
    if ismember(datTable.Struct(iSub),'Xclean2')
        load(char(datTable.Data(iSub)), 'rmch');
    elseif ismember(datTable.Struct(iSub),'XAl')
        load(char(datTable.Data(iSub)), 'badCh');
        rmch = badCh;
    end

    %get the patterns of the subject
    VHinv = allV1inv{iSub};
    VLinv = allV2inv{iSub};

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

%plot average magnitude patterns
avgPatternHigh = mean(abs(PatternsHigh),2);
avgPatternLow = mean(abs(PatternsLow),2);

nexttile(timeIdx+3,[1,1])
topoplot(avgPatternHigh, chanlocs,'maplimits','minmax','electrodes','off');
colorcet('L3','reverse',1)
title("High",'FontSize',14)
nexttile(timeIdx+7,[1,1])
topoplot(avgPatternLow, chanlocs,'maplimits','minmax','electrodes','off');
colorcet('L3','reverse',1)
title("Low",'FontSize',14)

%%
exportgraphics(gcf,['C:\Users\BNPPC08\Dropbox\CSP_manuscript\Figures for paper\fig6_TimeWindows\','fig6_TimeWindows.pdf'],...
    'BackgroundColor','none','ContentType','vector' )
