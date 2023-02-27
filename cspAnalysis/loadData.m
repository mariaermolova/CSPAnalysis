function [Xclean2,chanlocs0,badCh,goodTrials,isort,y] = loadData(datTableSub)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
load(char(datTableSub.Data),char(datTableSub.Struct));
load(char(datTableSub.Data),'chanlocs0');
load(char(datTableSub.Data),'badCh');
if exist('XAl','var')
    Xclean2 = XAl;
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

if ismember(datTableSub.Subnum,[14,15,16])
    chanlocs = chanlocs0(~badCh);
    badCh2 = ismember({chanlocs.labels},{'P3','P4'});
    Xclean2 = Xclean2(~badCh2,:,:);
elseif ismember(datTableSub.Subnum,[17])
    chanlocs = chanlocs0(~badCh);
    badCh2 = ismember({chanlocs.labels},{'P3','P4','P8','C4'});
    Xclean2 = Xclean2(~badCh2,:,:);
end

end

