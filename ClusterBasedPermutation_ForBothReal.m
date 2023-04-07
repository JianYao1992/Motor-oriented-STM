%% Cluster-based permuatation test for licking rate between two groups

function [SigTime_Real,SigTimeBelowChance_Real] = ClusterBasedPermutation_ForBothReal(data1,data2)

%% Difference in recording data
RealDiffBetweenTwoGroups = mean(data1) - mean(data2);

%% Difference in shuffled data
Data = [data1;data2];
ShuffleDiffBetweenTwoGroups = [];
for itr = 1:1000
    tempOrder  = randperm(size(Data,1));
    tempdata1 =  Data(tempOrder(1:size(data1,1)),:);
    tempdata2 =  Data(tempOrder(size(data1,1)+1:end),:);
    tempShuffleDiffBetweenTwoGroups = mean(tempdata1) - mean(tempdata2);
    ShuffleDiffBetweenTwoGroups = [ShuffleDiffBetweenTwoGroups; tempShuffleDiffBetweenTwoGroups];
end

%% Number of significant bins in each shuffle /// above and below chance level ///
temp95Percentile= prctile(ShuffleDiffBetweenTwoGroups,[2.5 97.5],1);
for i = 1:size(ShuffleDiffBetweenTwoGroups,1)  % shuffle time
    tempSigBinNum = 0; SigBinNum_Shuffle = []; tempSigBinNumBelowChance = 0; SigBinNumBelowChance_Shuffle = [];
    for j = 1:size(ShuffleDiffBetweenTwoGroups,2) % bin
        if ShuffleDiffBetweenTwoGroups (i,j) > temp95Percentile(2,j) % higher than 97.5%
            tempSigBinNum = tempSigBinNum + 1;
            if j == size(ShuffleDiffBetweenTwoGroups,2)
                SigBinNum_Shuffle = [SigBinNum_Shuffle tempSigBinNum];
            end
        elseif ShuffleDiffBetweenTwoGroups (i,j) < temp95Percentile(1,j) % lower than 2.5%
            tempSigBinNumBelowChance = tempSigBinNumBelowChance + 1;
            if j == size(ShuffleDiffBetweenTwoGroups,2)
                SigBinNumBelowChance_Shuffle = [SigBinNum_Shuffle tempSigBinNumBelowChance];
            end
        else
            if tempSigBinNum ~= 0
                SigBinNum_Shuffle = [SigBinNum_Shuffle tempSigBinNum];
                tempSigBinNum = 0;
            end
            if tempSigBinNumBelowChance ~= 0
                SigBinNumBelowChance_Shuffle = [SigBinNum_Shuffle tempSigBinNumBelowChance];
                tempSigBinNumBelowChance = 0;
            end
        end
    end
    if isempty(SigBinNum_Shuffle)
        MaxSigBinNum(i) = 0;
    else
        MaxSigBinNum(i) = max(SigBinNum_Shuffle);
    end
    if isempty(SigBinNumBelowChance_Shuffle)
        MaxSigBinNum_belowChance(i) = 0;
    else
        MaxSigBinNum_belowChance(i) = max(SigBinNumBelowChance_Shuffle);
    end
end
SigBinNum_Shuffle = prctile(MaxSigBinNum,95); SigBinNumBelowChance_Shuffle = prctile(MaxSigBinNum_belowChance,95);

%% Number of significant bins in real data /// above and below chance level ///
tempSigBinNum =0; SigBinNum = [];tempSigTime = []; SigTime= [];
tempSigBinNumBelowChance = 0; SigBinNumBelowChance = [];tempSigTimeBelowChance = []; SigTimeBelowChance= [];
for i = 1:length(RealDiffBetweenTwoGroups)
    if RealDiffBetweenTwoGroups (1,i) > temp95Percentile(2,i)
        tempSigTime = [tempSigTime i];
        tempSigBinNum = tempSigBinNum + 1;
        if i == length(RealDiffBetweenTwoGroups)
            SigTime = [SigTime {tempSigTime}];
            SigBinNum = [SigBinNum tempSigBinNum];
        end
    elseif RealDiffBetweenTwoGroups (1,i) < temp95Percentile(1,i)
        tempSigTimeBelowChance = [tempSigTimeBelowChance i];
        tempSigBinNumBelowChance = tempSigBinNumBelowChance + 1;
        if i == length(RealDiffBetweenTwoGroups)
            SigTimeBelowChance = [SigTimeBelowChance {tempSigTimeBelowChance}];
            SigBinNumBelowChance = [SigBinNumBelowChance tempSigBinNumBelowChance];
        end
    else
        if ~isempty(tempSigTime)
            SigTime = [SigTime {tempSigTime}];
        end
        if tempSigBinNum ~= 0
            SigBinNum = [SigBinNum tempSigBinNum];
        end
        tempSigBinNum = 0;
        tempSigTime = [];
        if ~isempty(tempSigTimeBelowChance)
            SigTimeBelowChance = [SigTimeBelowChance {tempSigTimeBelowChance}];
        end
        if tempSigBinNumBelowChance ~= 0
            SigBinNumBelowChance = [SigBinNumBelowChance tempSigBinNumBelowChance];
        end
        tempSigBinNumBelowChance = 0;
        tempSigTimeBelowChance = [];
    end
end

%% ID of bins with significant difference
SigTime_Real = []; SigTimeBelowChance_Real = [];
for i=1:length(SigBinNum)
    if SigBinNum(i) > SigBinNum_Shuffle
        SigTime_Real = [SigTime_Real SigTime{i}];
    end
end
for i=1:length(SigBinNumBelowChance)
    if SigBinNumBelowChance(i) > SigBinNumBelowChance_Shuffle
        SigTimeBelowChance_Real = [SigTimeBelowChance_Real SigTimeBelowChance{i}];
    end
end