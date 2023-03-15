function [eegLabeled,classId] = labelData(eeg,labels,mepSize,nSamples)
%Sort EEG trials according to MEP amplitudes and divide trials into 2 classes

mepHigh = mepSize(labels);
[~,highSorter] = sort(mepHigh,'descend');
eegHigh = eeg(:,:,labels);
eegHigh = eegHigh(:,:,highSorter);

mepLow = mepSize(~labels);
[~,lowSorter] = sort(mepLow,'ascend');
eegLow = eeg(:,:,~labels);
eegLow = eegLow(:,:,lowSorter);

eegLabeled = cat(3, eegHigh(:,:,1:nSamples), eegLow(:,:,1:nSamples));

classId = ones(nSamples*2,1);
classId(nSamples+1:end) = 2;

%     eeg = eeg(:,:,mepSorter([1:nSamples,end-(nSamples-1):end])); % matrix of all data (chan x time x trial)
%     classId = double(labels(mepSorter([1:nSamples,end-(nSamples-1):end])));
%     classId(classId == 0) = 2;

end