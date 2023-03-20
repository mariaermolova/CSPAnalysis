%% Plot phase lags of spatial patterns for each subject individually
projectPath = 'W:\Projects\2018-12 POSTHOCSOURCE Project\analysis_maria\CSPRepo';
datTable = readtable(fullfile(projectPath,'cspAnalysis','REFTEP_list.xlsx'), 'Basic', 1); %insert the name of the subject spreadsheet with paths
load(fullfile(projectPath,'Patterns','allPatterns_15-Mar-2023.mat')); %insert the name of the CSP pattern file
addpath(fullfile(projectPath,'interpretation'))
addpath('C:\Users\BNPPC08\Desktop\Maria\matlab\toolboxes\MatlabFns\Colourmaps')
addpath('C:\Users\BNPPC08\Desktop\Maria\matlab\toolboxes\eeglab14_1_2b')
eeglab
%% get accuracies for each subject
[accs] = extract_accuracies(fullfile(projectPath,'output','reftep_15-Mar-2023.mat')); %insert the name of the CSP output file
accs = accs*100;
%% select subjects to plot
subjects = [1:20]; %select subjects to plot
%% Interpolate missing channels (on complex patterns)
clear allPatternsLow allPatternsHigh

for idxSub = 1:length(subjects)
    
    iSub = subjects(idxSub);

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

    %prepare patterns of High condition
    EEG.data = VHinv(:,1);
    EEG.chanlocs = chanlocs0(~rmch);
    EEG.pnts = 1;
    EEG.trials = 1;
    EEG.nbchan = sum(~rmch);
    EEG = pop_interp( EEG, chanlocs0,'spherical');

    allPatternsHigh(:,idxSub) = EEG.data;

    %prepare patterns of Low condition
    EEG.data = VLinv(:,1);
    EEG.chanlocs = chanlocs0(~rmch);
    EEG.pnts = 1;
    EEG.trials = 1;
    EEG.nbchan = sum(~rmch);
    EEG = pop_interp( EEG, chanlocs0,'spherical');

    allPatternsLow(:,idxSub) = EEG.data;

end
%% Plot phases

figure, tiledlayout("flow")

for idxSub = 1:size(allPatternsHigh,2)

    nexttile

    %get the complex pattern
    subPattern = allPatternsHigh(:,idxSub);

    %find the index of a new reference channel
    refInd = find(ismember({chanlocs0.labels},'FC3'));

    %get the phase pattern
    compDiff = conj(subPattern).*subPattern(refInd);
    compDiff = angle(compDiff)*180/pi;
    
    topoplot(compDiff,chanlocs0,'electrodes','off','emarker2',{[refInd],'o','k',6});
    colormap(parula)
    title(append('High Sub ',num2str(idxSub),' Acc ',num2str(round(accs(idxSub)*100)),'%'))

end

