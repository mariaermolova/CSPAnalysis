%% filter CSP timecourses, dependencies: bbci toolbox, ...
clear
addpath('W:\Projects\2018-12 POSTHOCSOURCE Project\toolboxes\bbci')
startup_bbci_toolbox
projectPath = 'W:\Projects\2018-12 POSTHOCSOURCE Project\analysis_maria\CSPRepo';
addpath(fullfile(projectPath,'CSPAnalysis'))
load(fullfile(projectPath,'output','reftep_15-Mar-2023.mat'));

%% filter data with CSP filters
clear AllTimeCrsHDat AllTimeCrsLDat HPatHDatSpctr LPatHDatSpctr HPatLDatSpctr LPatLDatSpctr
subjects = 1:20;
for subnum = subjects

    cspOut = allSubOut{1,subnum};
    param = cspOut.param;
    data = cspOut.dataCV;

    param.toi = 251;
    param.toiWindow = 1244;

    [data1,data2] = prepDataforCSP(param,data,1,1,0,1);

    [covN] = [];

    % Get channel, time, and trial dimensions
    [nCh,~,nTr] = size(data1);

    % Generate CV indices
    idcClass1 = crossvalind('Kfold',nTr,param.nFolds);
    idcClass2 = crossvalind('Kfold',nTr,param.nFolds);

    % Divide data into training and test sets
    [XTrain1,XTrain2,XTest1,XTest2] = divideTrainTest(data1,data2,idcClass1,idcClass2,1);

    [~,ind] = max(mean(mean(allSubOut{1, subnum}.accCVmeans{1, 1},2),3));

    idxRegulCoef = param.hyperParamList(ind,2);
    idxnChCSP = param.hyperParamList(ind,1);

    XTrain1 = XTrain1(:,(311-125+1):311,:);
    XTrain2 = XTrain2(:,(311-125+1):311,:);

    % Calculate csp filters from training data
    [C,~,~] = calculateCSP(param,XTrain1,XTrain2,nCh,idxRegulCoef,idxnChCSP,covN);

    AllTimeCrsHDat{subnum} = pagemtimes(ctranspose(C),XTest1);%filter High data with both filters
    AllTimeCrsLDat{subnum} = pagemtimes(ctranspose(C),XTest2);%filter Low data with both filters
end

%% calculate Power Spectrum
clear HPatHDatSpctr LPatHDatSpctr HPatLDatSpctr LPatLDatSpctr
subjects = 1:20;
for subnum = subjects
    filteredData = [];
    filteredData.fs = 250;
    filteredData.x = permute(AllTimeCrsHDat{subnum}(1,:,:),[2 1 3]);

    opt.Scaling = 'db';

    filteredData_TEMP = proc_spectrum(filteredData, [8 30], opt);
    F = filteredData_TEMP.t;
    spectrum = filteredData_TEMP.x;
    HPatHDatSpctr(:,subnum) = mean(spectrum,3);

    filteredData.x = permute(AllTimeCrsHDat{subnum}(end,:,:),[2 1 3]);

    filteredData_TEMP = proc_spectrum(filteredData, [8 30], opt);
    F = filteredData_TEMP.t;
    spectrum = filteredData_TEMP.x;
    LPatHDatSpctr(:,subnum) = mean(spectrum,3);

    filteredData.x = permute(AllTimeCrsLDat{subnum}(1,:,:),[2 1 3]);

    filteredData_TEMP = proc_spectrum(filteredData, [8 30], opt);
    F = filteredData_TEMP.t;
    spectrum = filteredData_TEMP.x;
    HPatLDatSpctr(:,subnum) = mean(spectrum,3);

    filteredData.x = permute(AllTimeCrsLDat{subnum}(end,:,:),[2 1 3]);

    filteredData_TEMP = proc_spectrum(filteredData, [8 30], opt);
    F = filteredData_TEMP.t;
    spectrum = filteredData_TEMP.x;
    LPatLDatSpctr(:,subnum) = mean(spectrum,3);

end

%% Plot power spectrum
subjectsSign = [1:7,9:20];

SpectrtoPlot = HPatHDatSpctr(:,subjectsSign);
y = mean(SpectrtoPlot,2);
erbars = std(SpectrtoPlot,[],2)/sqrt(size(SpectrtoPlot,2));
y = y';
erbars = erbars';
curve1 = y + erbars;
curve2 = y - erbars;
x2 = [F, fliplr(F)];
inBetween = [curve1, fliplr(curve2)];

figure, tiledlayout(2,2)
nexttile
hold on
fill(x2, inBetween,[180 30 48]/255, 'FaceAlpha',0.1,'EdgeColor','w');
plot(F,y,'LineWidth',2,'Color',[180 30 48]/255)
xlabel('Frequency (Hz)')
ylabel('Power Spectrum (dB)')
ylim([6,21])
xlim([7,31])

SpectrtoPlot = HPatLDatSpctr(:,subjectsSign);
y = mean(SpectrtoPlot,2);
erbars = std(SpectrtoPlot,[],2)/sqrt(size(SpectrtoPlot,2));
y = y';
erbars = erbars';
curve1 = y + erbars;
curve2 = y - erbars;
x2 = [F, fliplr(F)];
inBetween = [curve1, fliplr(curve2)];

fill(x2, inBetween, [9 43 120]/255, 'FaceAlpha',0.1,'EdgeColor','w');
plot(F,y,':','LineWidth',1.4,'Color',[9 43 120]/255)

legend({'','High trials','','Low trials'})
title('High component')

f=gca;
f.FontSize = 14;
f.FontName = 'arial';
f.YTick = [5, 10, 15, 20];

nexttile
hold on

SpectrtoPlot = LPatLDatSpctr(:,subjectsSign);
y = mean(SpectrtoPlot,2);
erbars = std(SpectrtoPlot,[],2)/sqrt(size(SpectrtoPlot,2));
y = y';
erbars = erbars';
curve1 = y + erbars;
curve2 = y - erbars;
x2 = [F, fliplr(F)];
inBetween = [curve1, fliplr(curve2)];

fill(x2, inBetween, [9 43 120]/255, 'FaceAlpha',0.1,'EdgeColor','w');
plot(F,y,'LineWidth',2,'Color',[9 43 120]/255)

SpectrtoPlot = LPatHDatSpctr(:,subjectsSign);
y = mean(SpectrtoPlot,2);
erbars = std(SpectrtoPlot,[],2)/sqrt(size(SpectrtoPlot,2));
y = y';
erbars = erbars';
curve1 = y + erbars;
curve2 = y - erbars;
x2 = [F, fliplr(F)];
inBetween = [curve1, fliplr(curve2)];

fill(x2, inBetween,[180 30 48]/255, 'FaceAlpha',0.1,'EdgeColor','w');
plot(F,y,':','LineWidth',1.4,'Color',[180 30 48]/255)

legend({'','Low trials','','High trials'})
title('Low component')
xlabel('Frequency (Hz)')
ylabel('Power Spectrum (dB)')
ylim([6,21])
xlim([7,31])

f=gca;
f.FontSize = 14;
f.FontName = 'arial';
f.YTick = [5, 10, 15, 20];

%% Calculate time courses
zTimeHPatHDat = [];
zTimeHPatLDat = [];
zTimeLPatHDat = [];
zTimeLPatLDat = [];

[nCh,nTm,nTr] = size(AllTimeCrsHDat{subnum});
subjects = [1:20];

for subnum = subjects
    TimeCrsHPatHDat = mean(abs(squeeze(AllTimeCrsHDat{subnum}(1,:,:))),2);
    TimeCrsHPatLDat = mean(abs(squeeze(AllTimeCrsLDat{subnum}(1,:,:))),2);
    
    normHPat = [TimeCrsHPatHDat;TimeCrsHPatLDat];
    normHPat = normHPat/norm(normHPat,2);

    zTimeHPatHDat(:,subnum) = normHPat(1:nTm);
    zTimeHPatLDat(:,subnum) = normHPat(nTm+1:end);
   
    TimeCrsLPatHDat =  mean(abs(squeeze(AllTimeCrsHDat{subnum}(end,:,:))),2);
    TimeCrsLPatLDat = mean(abs(squeeze(AllTimeCrsLDat{subnum}(end,:,:))),2);

    normLPat = [TimeCrsLPatHDat;TimeCrsLPatLDat];
    normLPat = normLPat/norm(normLPat,2);

    zTimeLPatHDat(:,subnum) = normLPat(1:nTm);
    zTimeLPatLDat(:,subnum) = normLPat(nTm+1:end);

end

%% Plot time courses
nexttile

subjectsSign = [1:7,9:20];
times = -1252:4:-12;

TCoursetoPlot = zTimeHPatHDat(:,subjectsSign);
y = mean(TCoursetoPlot,2);
erbars = std(TCoursetoPlot,[],2)/sqrt(size(TCoursetoPlot,2));
y = y';
erbars = erbars';
curve1 = y + erbars;
curve2 = y - erbars;
x2 = [times, fliplr(times)];
inBetween = [curve1, fliplr(curve2)];

hold on
fill(x2, inBetween, [180 30 48]/255, 'FaceAlpha',0.1,'EdgeColor','w');
plot(times,y,'LineWidth',2,'Color',[180 30 48]/255)

TCoursetoPlot = zTimeHPatLDat(:,subjectsSign);
y = mean(TCoursetoPlot,2);
erbars = std(TCoursetoPlot,[],2)/sqrt(size(TCoursetoPlot,2));
y = y';
erbars = erbars';
curve1 = y + erbars;
curve2 = y - erbars;
x2 = [times, fliplr(times)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, [9 43 120]/255, 'FaceAlpha',0.1,'EdgeColor','w');
plot(times,y,':','LineWidth',1.4,'Color',[9 43 120]/255)

xline(-500,'--')
title('High component')
legend({'','High trials','','Low trials'},'Location','northeast')
ylabel('Signal envelope (normalised)')
xlabel('Time (ms)')
xlim([-1300 0])
ylim([0.03 0.05])

f=gca;
f.FontSize = 14;
f.FontName = 'arial';
f.YTick = [0.03, 0.04, 0.05];

nexttile
hold on
TCoursetoPlot = zTimeLPatLDat(:,subjectsSign);
y = mean(TCoursetoPlot,2);
erbars = std(TCoursetoPlot,[],2)/sqrt(size(TCoursetoPlot,2));
y = y';
erbars = erbars';
curve1 = y + erbars;
curve2 = y - erbars;
x2 = [times, fliplr(times)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, [9 43 120]/255, 'FaceAlpha',0.1,'EdgeColor','w');
plot(times,y,'LineWidth',2,'Color',[9 43 120]/255)

TCoursetoPlot = zTimeLPatHDat(:,subjectsSign);
y = mean(TCoursetoPlot,2);
erbars = std(TCoursetoPlot,[],2)/sqrt(size(TCoursetoPlot,2));
y = y';
erbars = erbars';
curve1 = y + erbars;
curve2 = y - erbars;
x2 = [times, fliplr(times)];
inBetween = [curve1, fliplr(curve2)];

fill(x2, inBetween, [180 30 48]/255, 'FaceAlpha',0.1,'EdgeColor','w');
plot(times,y,':','LineWidth',1.4,'Color',[180 30 48]/255)

xline(-500,'--')
title('Low component')
legend({'','Low trials','','High trials'},'Location','northeast')
ylabel('Signal envelope (normalised)')
xlabel('Time (ms)')
xlim([-1300 0])
ylim([0.03 0.05])

f=gca;
f.FontSize = 14;
f.FontName = 'arial';
f.YTick = [0.03, 0.04, 0.05];
%%
exportgraphics(gcf,['C:\Users\BNPPC08\Dropbox\CSP_manuscript\Figures for paper\fig4_PSD_TC\','fig4_PSD_TC.pdf'],...
    'BackgroundColor','none','ContentType','vector' )
