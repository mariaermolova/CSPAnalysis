function [eeg,chanlocs0,badCh,goodTrials,mepsize,labels] = loadData(dataTableSub)
%load EEG data, channel locations, removed channel indices, retained trial
%indices, MEP sorting indices, class labels
load(char(dataTableSub.Data),char(dataTableSub.Struct));
load(char(dataTableSub.Data),'chanlocs0');
load(char(dataTableSub.Data),'badCh');
if exist('XAl','var')
    eeg = XAl;
    load(char(dataTableSub.Data),'badTrEMG');
    load(char(dataTableSub.Data),'badTr');
    goodTrials=logical(~badTr);
    goodTrials(badTrEMG)=0;
else
    eeg = Xclean2;
    load(char(dataTableSub.Data),'badTrInds');
    load(char(dataTableSub.Data),'goodTrials');
    goodTrials = logical(goodTrials);
    goodTrials(badTrInds)=0;
end

load(char(dataTableSub.Data),'mepsize')
load(char(dataTableSub.Data),'labels')


end