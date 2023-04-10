function Distance = TrajectoryDistanceCalculation(UnitBinFR_S1,UnitBinFR_S2,TrialNumForEachCondition,UnitNumInAnalysis,num_resample_runs,IsShuffle,BaseLen,SampLen,DelayLen)

Distance = []; 
TimeGain = 10;
for i = 1:num_resample_runs 
    UnitBinFR = cellfun(@(x,y) vertcat(x,y),UnitBinFR_S1,UnitBinFR_S2,'UniformOutput',0);
    
    %% Sample neurons randomly
    temp = randperm(numel(UnitBinFR));
    TargetUnitID = temp(1:UnitNumInAnalysis);
    tempUnitBinFR = UnitBinFR(:,TargetUnitID);
    
    %% FR of sampled neurons
    FR_S1 = cell(1,numel(tempUnitBinFR));
    FR_S2 = cell(1,numel(tempUnitBinFR));
    for iUnit = 1:numel(tempUnitBinFR)
        if IsShuffle == 1
            TrialIndex = 1:size(tempUnitBinFR{iUnit},1);
            temp = randperm(length(TrialIndex));
            TrialIndex1 = temp(1:size(UnitBinFR_S1{iUnit},1));
            TrialIndex2 = temp(size(UnitBinFR_S1{iUnit},1)+1:end);
        elseif IsShuffle == 0
            TrialIndex1 = 1:size(UnitBinFR_S1{iUnit},1);
            TrialIndex2 = size(UnitBinFR_S1{iUnit},1)+1:size(tempUnitBinFR{iUnit},1);
        end
        % sample S1 trials randomly
        temp1 = randperm(length(TrialIndex1));
        if ~isempty(TrialNumForEachCondition)
            NewTrialIndex1 = TrialIndex1(temp1(1:TrialNumForEachCondition));
        else
            NewTrialIndex1 = TrialIndex1;
        end
        FR_S1{1,iUnit} = tempUnitBinFR{iUnit}(NewTrialIndex1,:);
        % sample S2 trials randomly
        temp2 = randperm(length(TrialIndex2));
        if ~isempty(TrialNumForEachCondition)
            NewTrialIndex2 = TrialIndex2(temp2(1:TrialNumForEachCondition));
        else
            NewTrialIndex2 = TrialIndex2;
        end
        FR_S2{1,iUnit} = tempUnitBinFR{iUnit}(NewTrialIndex2,:);
    end
    
    %% Averaged FR over the time course of the trial
    meanFR_S1 = cellfun(@mean,FR_S1,'UniformOutput',0);
    meanFR_S2 = cellfun(@mean,FR_S2,'UniformOutput',0);
    
    %% Averaged FR over the time course of sample and delay period
    meanFR_SampDelay_S1 = cellfun(@(x) x(:,BaseLen*TimeGain+1:(BaseLen+SampLen+DelayLen)*TimeGain),meanFR_S1,'UniformOutput',0);
    meanFR_SampDelay_S2 = cellfun(@(x) x(:,BaseLen*TimeGain+1:(BaseLen+SampLen+DelayLen)*TimeGain),meanFR_S2,'UniformOutput',0);
    
    %% Eigenvector
    meanFR_S1 = (vertcat(meanFR_S1{:}))'; % whole trial
    meanFR_S2 = (vertcat(meanFR_S2{:}))';
    meanFR = [meanFR_S1; meanFR_S2];
    meanFR_SampDelay_S1 = (vertcat(meanFR_SampDelay_S1{:}))'; % sample and delay period
    meanFR_SampDelay_S2 = (vertcat(meanFR_SampDelay_S2{:}))';
    meanFR_SampDelay = [meanFR_SampDelay_S1; meanFR_SampDelay_S2];
    CenterMeanFR_SampDelay = repmat(mean(meanFR_SampDelay),size(meanFR,1),1);
    [EigenVector,score,EigenValue,~,explained] = pca(meanFR_SampDelay);
    
    %% Distance between trajectories
    FR_MeanSubtracted = meanFR - CenterMeanFR_SampDelay;
    Score_WholeTrial = FR_MeanSubtracted*EigenVector;
    TotalBinNum = size(FR_MeanSubtracted,1);
    Score_S1 = Score_WholeTrial(1:TotalBinNum/2,1:20);
    Score_S2 = Score_WholeTrial(TotalBinNum/2+1:TotalBinNum,1:20);
    Distance(i,:) = sum((Score_S1-Score_S2).^2,2).^0.5;
end

