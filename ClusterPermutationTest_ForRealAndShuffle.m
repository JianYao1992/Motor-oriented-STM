%% This code is used to perform cluster-based permutation test.

function [SigTime_Real,SigTimebelowChance_Real] = ClusterPermutationTest_ForRealAndShuffle(Real,Shuffle)

temp95Percentile= prctile(Shuffle,[2.5 97.5],1);  
for i = 1:size(Shuffle,1)  % shuffle time
    tempSigBinNum = 0; SigBinNum = []; tempSigBinNum_belowChance = []; SigBinNum_belowChance =[];
    for j = 1:size(Shuffle,2) % bin
        if Shuffle (i,j) > temp95Percentile(2,j) % higher than 97.5%
            tempSigBinNum = tempSigBinNum + 1;
            if j == size(Shuffle,2)
                SigBinNum = [SigBinNum tempSigBinNum];
            end 
        elseif Shuffle (i,j) < temp95Percentile(1,j) % lower than 2.5%
            tempSigBinNum_belowChance = tempSigBinNum_belowChance + 1;
            if j == size(Shuffle,2)
                SigBinNum_belowChance = [SigBinNum_belowChance tempSigBinNum_belowChance];
            end 
        else
            if tempSigBinNum ~= 0
                SigBinNum = [SigBinNum tempSigBinNum];
            end
            tempSigBinNum = 0;
            if tempSigBinNum_belowChance ~= 0
                SigBinNum_belowChance = [SigBinNum_belowChance tempSigBinNum_belowChance];
            end
            tempSigBinNum_belowChance = 0;
        end
    end
    if isempty(SigBinNum)
        MaxSigBinNum(i) = 0;
    else
        MaxSigBinNum(i) = max(SigBinNum);
    end
    if isempty(SigBinNum_belowChance)
        MaxSigBinNum_belowChance(i) = 0;
    else
        MaxSigBinNum_belowChance(i) = max(SigBinNum_belowChance);
    end
end 
SigBinNum = prctile(MaxSigBinNum,95);  SigBinNum_belowChance = prctile(MaxSigBinNum_belowChance,95); 
MeanRealResults = mean(Real,1); 
tempSigBinNum =0; SigBinNum_Real = [];tempSigTime = []; SigTime= [];
tempSigBinNum_belowChance =0; SigBinNumbelowChance_Real = [];tempSigTime_belowChance = []; SigTime_belowChance= [];
for i = 1:length(MeanRealResults)
    if MeanRealResults (1,i) > temp95Percentile(2,i)
        tempSigTime = [tempSigTime i];
        tempSigBinNum = tempSigBinNum + 1;
        if i == length(MeanRealResults)
            SigTime = [SigTime {tempSigTime}];
            SigBinNum_Real = [SigBinNum_Real tempSigBinNum];
        end
    elseif MeanRealResults (1,i) < temp95Percentile(1,i)
        tempSigTime_belowChance = [tempSigTime_belowChance i];
        tempSigBinNum_belowChance = tempSigBinNum_belowChance + 1;
        if i == length(MeanRealResults)
            SigTime_belowChance = [SigTime_belowChance {tempSigTime_belowChance}];
            SigBinNumbelowChance_Real = [SigBinNumbelowChance_Real tempSigBinNum_belowChance];
        end
    else
        if ~isempty(tempSigTime)
            SigTime = [SigTime {tempSigTime}];
        end
        if tempSigBinNum ~= 0
            SigBinNum_Real = [SigBinNum_Real tempSigBinNum];
        end
        tempSigBinNum = 0;
        tempSigTime = [];
        if ~isempty(tempSigTime_belowChance)
            SigTime_belowChance = [SigTime_belowChance {tempSigTime_belowChance}];
        end
        if tempSigBinNum_belowChance ~= 0
            SigBinNumbelowChance_Real = [SigBinNumbelowChance_Real tempSigBinNum_belowChance];
        end
        tempSigBinNum_belowChance = 0;
        tempSigTime_belowChance = [];
    end
end
SigTime_Real = []; %%% Get the significant time bin in real decoding or selectivity result.
for i=1:length(SigBinNum_Real)
    if SigBinNum_Real(i) > SigBinNum
        SigTime_Real = [SigTime_Real SigTime{i}];
    end
end
SigTimebelowChance_Real = []; %%% Get the significant time bin in real decoding or selectivity result.
for i=1:length(SigBinNumbelowChance_Real)
    if SigBinNumbelowChance_Real(i) > SigBinNum_belowChance
        SigTimebelowChance_Real = [SigTimebelowChance_Real SigTime_belowChance{i}];
    end
end