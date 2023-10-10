function [AUC,ShuffledAUC] = CalculatePopulationFcAUC(Bin,PairID,TrialMarker,TrialSpikeTime,ShuffleTimes,WhetherShuffle)

AUC = cell(1,numel(Bin));
ShuffledAUC = cell(0);
IsSigAUC = cell(1,numel(Bin));
for bin = 1:length(Bin) 
    for iPair = 1:size(PairID{bin},1)
        fprintf('Analyze %d/%d_%dth delay bin\n',iPair,size(PairID{bin},1),bin);
        Unit1 = PairID{bin}(iPair,1); 
        Unit2 = PairID{bin}(iPair,2);
        TrialID_Go_Unit1 = find(TrialMarker{Unit1,2}(:,2)==1); 
        TrialID_NoGo_Unit1 = find(TrialMarker{Unit1,2}(:,2)==2);
        TrialID_Go_Unit2 = find(TrialMarker{Unit2,2}(:,2)==1);
        TrialID_NoGo_Unit2 = find(TrialMarker{Unit2,2}(:,2)==2);
        SpikeTime_Unit1 = TrialSpikeTime{Unit1,2};
        SpikeTime_Unit2 = TrialSpikeTime{Unit2,2};
        [tempAUC,tempShuffledAUC,tempIsSigAUC] = CalculateIndividualFcAUC(bin,SpikeTime_Unit1,SpikeTime_Unit2,TrialID_Go_Unit1,TrialID_NoGo_Unit1,TrialID_Go_Unit2,TrialID_NoGo_Unit2,ShuffleTimes,WhetherShuffle);
        AUC{1,bin}(end+1,1) = tempAUC;
        ShuffledAUC{1,bin}{iPair,1} = tempShuffledAUC;
        if WhetherShuffle == 1
            IsSigAUC{1,bin}(end+1,:) = tempIsSigAUC;
        else
            IsSigAUC{1,bin}(end+1,:) = 0;
        end
    end
end
AUC = cellfun(@(x,y,z) [x y z],PairID,AUC,IsSigAUC,'UniformOutput',0);
