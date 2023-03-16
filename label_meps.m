% label data as High or Low condition based on MEP amplitude relative to a
% moving median MEP. Select EMG channel with higher median MEP amplitude.

dataTable = readtable('REFTEP_list.xlsx', 'Basic', 1); %load data path info

subjects = [1:20]; %select subjects for analysis
%%

for subnum = subjects

    dataTableSub = dataTable(subnum,:);
    subject = dataTableSub.Data;

    load(subject{:}, 'AmpsAl')

    if median(AmpsAl(:,1))>median(AmpsAl(:,2))
        mepsize=AmpsAl(:,1);
    elseif median(AmpsAl(:,1))<median(AmpsAl(:,2))
        mepsize=AmpsAl(:,2);
    end

    mepmedian = movmedian(mepsize,150);

    labels = mepsize > mepmedian;

    save(subject{:},'labels','mepsize','-append')

end

