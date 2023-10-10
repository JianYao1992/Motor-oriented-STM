function tempRate = CalculateSinglePairEventRateofFunctionalCoupling(DelayBinID,SpikeStamp1,SpikeStamp2,Marker_Unit1,Marker_Unit2)

BeforeDelay = 5;

%% Spike raster in target trials of two neurons
SpikeStamp_Unit1 = SpikeStamp1(:,Marker_Unit1); SpikeStamp_Unit2 = SpikeStamp2(:,Marker_Unit2); Num = [];

%% FCSP events number during the target bin
for iTrial = 1:length(SpikeStamp_Unit1)
    % FCSP events number of the leading neuron
    tempStamp_Unit1 = [];
    for iSpike = 1:length(SpikeStamp_Unit1{iTrial})
        if ~isempty(find(SpikeStamp_Unit2{iTrial}(:,1)-SpikeStamp_Unit1{iTrial}(iSpike) <= 0.01 & SpikeStamp_Unit2{iTrial}(:,1)-SpikeStamp_Unit1{iTrial}(iSpike) > 0.002))
            tempStamp_Unit1 = [tempStamp_Unit1 SpikeStamp_Unit1{iTrial}(iSpike)];
        end
    end
    tempStamp_Unit1(tempStamp_Unit1 < BeforeDelay+DelayBinID-1 | tempStamp_Unit1 >= BeforeDelay+DelayBinID) = [];
    % FCSP events number of the following neuron
    tempStamp_Unit2 = [];
    for iSpike = 1:length(SpikeStamp_Unit2{iTrial})
        if ~isempty(find(SpikeStamp_Unit2{iTrial}(iSpike)-SpikeStamp_Unit1{iTrial}(:,1) <= 0.01 & SpikeStamp_Unit2{iTrial}(iSpike)-SpikeStamp_Unit1{iTrial}(:,1) > 0.002))
            tempStamp_Unit2 = [tempStamp_Unit2 SpikeStamp_Unit2{iTrial}(iSpike)];
        end
    end
    tempStamp_Unit2(tempStamp_Unit2 <= BeforeDelay+DelayBinID-1 | tempStamp_Unit2 > BeforeDelay+DelayBinID) = [];
    Num(end+1,1) = min([length(tempStamp_Unit1) length(tempStamp_Unit2)]);
end
tempRate = mean(Num);


