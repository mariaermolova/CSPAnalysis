%% plot magnitides of spatial patterns for each subject individually

projectPath = 'W:\Projects\2018-12 POSTHOCSOURCE Project\analysis_maria\CSPRepo';
datTable = readtable(fullfile(projectPath,'cspAnalysis','REFTEP_list.xlsx'), 'Basic', 1); %insert the name of the subject spreadsheet with paths
addpath('C:\Users\BNPPC08\Desktop\Maria\matlab\toolboxes\MatlabFns\Colourmaps')
addpath('C:\Users\BNPPC08\Desktop\Maria\matlab\toolboxes\eeglab14_1_2b')
addpath(fullfile(projectPath,'interpretation'))
load(fullfile(projectPath,'Patterns','allPatterns_15-Mar-2023.mat')) %insert the name of the CSP pattern file
eeglab
%% get accuracies for each subject
[accs] = extract_accuracies(fullfile(projectPath,'output','reftep_15-Mar-2023.mat')); %insert the name of the CSP output file
accs = accs*100;
%% select subjects to plot
subjects = [1:20]; 
%% plot High patterns
figure
t = tiledlayout('flow');
t.TileSpacing = 'compact';
for subnum = subjects

    %get channel locations of the subject
    load(char(datTable.Data(subnum)), 'chanlocs0');
    if ismember(datTable.Struct(subnum),'Xclean2')
        load(char(datTable.Data(subnum)), 'rmch');
    elseif ismember(datTable.Struct(subnum),'XAl')
        load(char(datTable.Data(subnum)), 'badCh');
        rmch = badCh;
    end

    %get spatial pattern of the subject
    VHinv = allV1inv{subnum};
    VHinv = abs(VHinv(:,1)); %magnitude of the first pattern

    nexttile
    topoplot(VHinv, chanlocs0(~rmch),'electrodes','off','maplimits','maxmin');
    colorcet('L3','reverse',1)
    title(append('High Sub ',num2str(subnum),' Acc ',num2str(round(accs(subnum))),'%'))

end
%% plot Low patterns
figure
t = tiledlayout('flow');
t.TileSpacing = 'compact';
for subnum = subjects

    %get channel locations of the subject
    load(char(datTable.Data(subnum)), 'chanlocs0');
    if ismember(datTable.Struct(subnum),'Xclean2')
        load(char(datTable.Data(subnum)), 'rmch');
    elseif ismember(datTable.Struct(subnum),'XAl')
        load(char(datTable.Data(subnum)), 'badCh');
        rmch = badCh;
    end

    %get spatial pattern of the subject
    VLinv = allV2inv{subnum};
    VLinv = abs(VLinv(:,1)); %magnitude of the first pattern

    nexttile
    topoplot(VLinv, chanlocs0(~rmch),'electrodes','off','maplimits','maxmin');
    colorcet('L3','reverse',1)
    title(append('Low Sub ',num2str(subnum),' Acc ',num2str(round(accs(subnum))),'%'))

end