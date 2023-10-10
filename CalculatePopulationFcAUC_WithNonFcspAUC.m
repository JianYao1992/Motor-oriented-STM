function AUC = CalculatePopulationFcAUC_WithNonFcspAUC(Bin,PairID,TrialMarker,TrialSpikeTime,SamplingTimes)

AUC = cell(1,numel(Bin));
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
        [tempAUC_Fcsp,tempAUC_NonFcsp] = CalculateIndividualFcAUC_WithNonFcspAUC(bin,SpikeTime_Unit1,SpikeTime_Unit2,TrialID_Go_Unit1,TrialID_NoGo_Unit1,TrialID_Go_Unit2,TrialID_NoGo_Unit2,SamplingTimes);
        AUC{1,bin}(end+1,:) = horzcat(tempAUC_Fcsp,tempAUC_NonFcsp);
    end
end
AUC = cellfun(@(x,y,z) [x y],PairID,AUC,'UniformOutput',0);
