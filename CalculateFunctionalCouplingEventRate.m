function FCRate = CalculateFunctionalCouplingEventRate(PairID,TrialMarker,TrialSpikeTime)

tempFCRate = cell(0);
for bin = 1:length(PairID)
    for iTrialType = 1:length(PairID{bin})
        switch iTrialType
            case 1
                TrialType=  'GoTrials';
            case 2
                TrialType = 'NogoTrials';
        end
        for iPair = 1:size(PairID{bin}{iTrialType},1)
            fprintf('%d/%d_%s_%dth delay bin\n',iPair,size(PairID{bin}{iTrialType},1),TrialType,bin);
            Unit1 = PairID{bin}{iTrialType}(iPair,1); Unit2 = PairID{bin}{iTrialType}(iPair,2);
            if strcmp(TrialType,'GoTrials')
                TrialID_Unit1 = find(TrialMarker{Unit1,2}(:,2)==1); 
                TrialID_Unit2 = find(TrialMarker{Unit2,2}(:,2)==1);
            else
                TrialID_Unit1 = find(TrialMarker{Unit1,2}(:,2)==2);
                TrialID_Unit2 = find(TrialMarker{Unit2,2}(:,2)==2);
            end
            SpikeTime_Unit1 = TrialSpikeTime{Unit1,2}; SpikeTime_Unit2 = TrialSpikeTime{Unit2,2};
            tempFCRate{1,bin}{iTrialType,1}(iPair,1) = CalculateSinglePairEventRateofFunctionalCoupling(bin,SpikeTime_Unit1,SpikeTime_Unit2,TrialID_Unit1,TrialID_Unit2);
        end
    end
end
FCRate = cell(0);
for bin = 1:length(PairID)
    FCRate{bin} = cellfun(@(x,y) horzcat(x,y),PairID{bin},tempFCRate{bin},'UniformOutput',0);
end

