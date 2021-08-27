function [covN] = noiseCovariance(param,dataCV)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

bpFilt = designfilt('highpassiir','FilterOrder',param.bpfiltparam.FilterOrder, ...
    'HalfPowerFrequency',1, ...
    'SampleRate',param.SampleRate,'DesignMethod',param.bpfiltparam.DesignMethod);
% fvtool(bpFilt)
% bpFilt = designfilt('bandpassiir','FilterOrder',param.bpfiltparam.FilterOrder, ...
%     'HalfPowerFrequency1',1,'HalfPowerFrequency2',100, ...
%     'SampleRate',param.SampleRate,'DesignMethod',param.bpfiltparam.DesignMethod);
% bpFilt = designfilt('highpassiir','FilterOrder',param.bpfiltparam.FilterOrder, ...
%     'HalfPowerFrequency',2, ...
%     'SampleRate',param.SampleRate,'DesignMethod',param.bpfiltparam.DesignMethod);

padSize = round(param.SampleRate);
dataPerm = permute(dataCV,[2 1 3]); %time x chan x trial
paddedData  = padarray(dataPerm,padSize,'symmetric','both'); %pad
% filteredData = filtfilt(bpFilt,paddedData);
filteredData = filter(bpFilt,paddedData);

% hilbert transform
hilbData = hilbert(filteredData); % hilbert transform
hilbData = hilbData(padSize+1:end-padSize,:,:);
hilbData = permute(hilbData,[2 1 3]); %chan x time x  trial
hilbData = hilbData(:,param.toi:param.toi+param.toiWindow,:); % choose time

[nCh,~,nTr] = size(hilbData);
covNTr = zeros(nCh,nCh,nTr);
for trial = 1:nTr
    trial1 = squeeze(hilbData(:,:,trial));
    covNTr(:,:,trial) = (trial1*trial1')/trace(trial1*trial1');   % S1~[C x C]
end
covN = mean(covNTr,3);
covN = real(covN);

end

