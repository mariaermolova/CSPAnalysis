function [data1,data2] = prepDataforCSP(param,data,freqIdx,toiIdx,trialFlag,classFlag)
%bp filter, hilbert transform, cut to toi, separate conditions, equalize
%trial numbers

% bp filter
bpFilt = designfilt('bandpassiir','FilterOrder',param.bpfiltparam.FilterOrder, ...
    'HalfPowerFrequency1',param.freq(freqIdx),'HalfPowerFrequency2',param.freq(freqIdx)+param.freqband(freqIdx), ...
    'SampleRate',param.SampleRate,'DesignMethod',param.bpfiltparam.DesignMethod);
% fvtool(bp_filt)
% padSize = ceil(bpFilt.FilterOrder*1.5);
% padSize = round(param.SampleRate/mean([param.freq(freqIdx),param.freq(freqIdx)+param.freqband(freqIdx)]));
padSize = round(param.SampleRate);
dataPerm = permute(data,[2 1 3]); %time x chan x trial
paddedData  = padarray(dataPerm,padSize,'symmetric','both'); %pad
% filteredData = filtfilt(bpFilt,paddedData);
filteredData = filter(bpFilt,paddedData);

resampledData = [];
for trIdx = 1:size(filteredData,3)
    resampledData(:,:,trIdx) = resample(filteredData(:,:,trIdx),250,1000);
end
padSize = padSize/4;
toi1 = round(param.toi(toiIdx)/4);
toi2 = round((param.toi(toiIdx)+param.toiWindow)/4);

% hilbert transform
hilbData = hilbert(resampledData); % hilbert transform
hilbData = hilbData(padSize+1:end-padSize,:,:);
hilbData = permute(hilbData,[2 1 3]); %chan x time x trial
hilbData = hilbData(:,toi1:toi2,:); % choose time

if classFlag == 1
    data1 = hilbData(:,:,param.class.CV==1);
    data2 = hilbData(:,:,param.class.CV==2);
elseif classFlag == 2
    data1 = hilbData(:,:,param.class.Val==1);
    data2 = hilbData(:,:,param.class.Val==2);
else
    data1 = hilbData;
end

% solve inequal trial sizes
if trialFlag == 1
    nTr_1_old = size(data1,3);
    nTr_2_old = size(data2,3);
    
    if nTr_1_old < nTr_2_old
        data2 = permute(data2,[3 1 2]);
        data2 = datasample(data2,nTr_1_old,'Replace',false);
        data2 = permute(data2,[2 3 1]);
    elseif nTr_1_old > nTr_2_old
        data1 = permute(data1,[3 1 2]);
        data1 = datasample(data1,nTr_2_old,'Replace',false);
        data1 = permute(data1,[2 3 1]);
    end
end

end

