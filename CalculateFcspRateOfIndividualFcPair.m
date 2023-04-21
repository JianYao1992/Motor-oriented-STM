function FCrate = CalculateFcspRateOfIndividualFcPair(BaseLen,SampLen,DelayBinID,TimeStamp_Unit1,TimeStamp_Unit2)

BeforeDelay = BaseLen+SampLen;
FCrate = [];

%% Num of FCSP events of delay
for iTrial = 1:numel(TimeStamp_Unit1)
    % num of FCSP events in the leading neuron
    tempStamp_Unit1 = [];
    for iSpike = 1:length(TimeStamp_Unit1{iTrial})
        if ~isempty(find(TimeStamp_Unit2{iTrial}(:,1)-TimeStamp_Unit1{iTrial}(iSpike)<=0.01 & TimeStamp_Unit2{iTrial}(:,1)-TimeStamp_Unit1{iTrial}(iSpike)>0.002))
            tempStamp_Unit1 = [tempStamp_Unit1 TimeStamp_Unit1{iTrial}(iSpike)];
        end
    end
    tempStamp_Unit1(tempStamp_Unit1<BeforeDelay+DelayBinID-1 | tempStamp_Unit1>=BeforeDelay+DelayBinID) = [];
    % num of FCSP events in the following neuron
    tempStamp_Unit2 = [];
    for iSpike = 1:length(TimeStamp_Unit2{iTrial})
        if ~isempty(find(TimeStamp_Unit2{iTrial}(iSpike)-TimeStamp_Unit1{iTrial}(:,1)<=0.01 & TimeStamp_Unit2{iTrial}(iSpike)-TimeStamp_Unit1{iTrial}(:,1)>0.002))
            tempStamp_Unit2 = [tempStamp_Unit2 TimeStamp_Unit2{iTrial}(iSpike)];
        end
    end
    tempStamp_Unit2(tempStamp_Unit2<=BeforeDelay+DelayBinID-1 | tempStamp_Unit2>BeforeDelay+DelayBinID) = [];
    FCrate(end+1,1) = min([numel(tempStamp_Unit1) numel(tempStamp_Unit2)]);
end
