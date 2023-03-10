function [eeg,chanlocs0,badCh,goodTrials,isort,y] = loadData(datTableSub)
%load EEG data, channel locations, removed channel indices, retained trial
%indices, MEP sorting indices, class labels
load(char(datTableSub.Data),char(datTableSub.Struct));
load(char(datTableSub.Data),'chanlocs0');
load(char(datTableSub.Data),'badCh');
if exist('XAl','var')
    eeg = XAl;
    load(char(datTableSub.Data),'badTrEMG');
    load(char(datTableSub.Data),'badTr');
    goodTrials=logical(~badTr);
    goodTrials(badTrEMG)=0;
else
    load(char(datTableSub.Data),'badTrInds');
    load(char(datTableSub.Data),'goodTrials');
    goodTrials = logical(goodTrials);
    goodTrials(badTrInds)=0;
end

load(char(datTableSub.Data),'isort')
load(char(datTableSub.Data),'y')


end