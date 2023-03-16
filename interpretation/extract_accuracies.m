function [accs] = extract_accuracies(filename)
%Extract mean accuracies from CSP results

load(filename)

for idxSub = 1:length(allSubOut)
    for idxTime = 1:size(allSubOut{1, idxSub}.ACC,1)
        for idxFreq = 1:size(allSubOut{1, idxSub}.ACC,2)
            accs(idxSub,idxTime, idxFreq) = mean(allSubOut{1, idxSub}.ACC{idxTime, idxFreq}.mean);
        end
    end
end

accs = squeeze(accs);

end