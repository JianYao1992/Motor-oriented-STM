function [AUC_Fcsp,AUC_NonFcsp] = CalculateIndividualFcAUC_WithNonFcspAUC(BinID,SpikeStamp1,SpikeStamp2,MarkerS1_Unit1,MarkerS2_Unit1,MarkerS1_Unit2,MarkerS2_Unit2,SamplingTimes)

BeforeDelay = 5;

%% Spike rasters of two neurons in S1 and S2 trials 
SpikeStampS1_Unit1 = SpikeStamp1(:,MarkerS1_Unit1);
SpikeStampS2_Unit1 = SpikeStamp1(:,MarkerS2_Unit1);
SpikeStampS1_Unit2 = SpikeStamp2(:,MarkerS1_Unit2); 
SpikeStampS2_Unit2 = SpikeStamp2(:,MarkerS2_Unit2);

%% Determine the minimal number of Go trials and NoGo trials
MinS1TriNum = min(length(SpikeStampS1_Unit1),length(SpikeStampS1_Unit2));
MinS2TriNum = min(length(SpikeStampS2_Unit1),length(SpikeStampS2_Unit2));
SpikeStampS1_Unit1 = SpikeStampS1_Unit1(:,1:MinS1TriNum);
SpikeStampS2_Unit1 = SpikeStampS2_Unit1(:,1:MinS2TriNum);
SpikeStampS1_Unit2 = SpikeStampS1_Unit2(:,1:MinS1TriNum); 
SpikeStampS2_Unit2 = SpikeStampS2_Unit2(:,1:MinS2TriNum);

%% FCSP number during delay of S1 trials
SpikeNum_S1 = []; FcspNum_S1 = [];
for iTrial = 1:length(SpikeStampS1_Unit1)
    % FCSP number of leading neuron
    tempFcspStamp_Unit1 = [];
    for iSpike = 1:length(SpikeStampS1_Unit1{iTrial})
        if ~isempty(find(SpikeStampS1_Unit2{iTrial}(:,1)-SpikeStampS1_Unit1{iTrial}(iSpike)<=0.01 & SpikeStampS1_Unit2{iTrial}(:,1)-SpikeStampS1_Unit1{iTrial}(iSpike)>0.002))
            tempFcspStamp_Unit1 = [tempFcspStamp_Unit1 SpikeStampS1_Unit1{iTrial}(iSpike)];
        end
    end
    tempStamp_Unit1 = SpikeStampS1_Unit1{iTrial};
    tempStamp_Unit1(tempStamp_Unit1<=BeforeDelay+BinID-1 | tempStamp_Unit1>BeforeDelay+BinID) = [];
    tempFcspStamp_Unit1(tempFcspStamp_Unit1<BeforeDelay+BinID-1 | tempFcspStamp_Unit1>=BeforeDelay+BinID) = [];
    % FCSP number of following neuron
    tempFcspStamp_Unit2 = [];
    for iSpike = 1:length(SpikeStampS1_Unit2{iTrial})
        if ~isempty(find(SpikeStampS1_Unit2{iTrial}(iSpike)-SpikeStampS1_Unit1{iTrial}(:,1)<=0.01 & SpikeStampS1_Unit2{iTrial}(iSpike)-SpikeStampS1_Unit1{iTrial}(:,1)>0.002))
            tempFcspStamp_Unit2 = [tempFcspStamp_Unit2 SpikeStampS1_Unit2{iTrial}(iSpike)];
        end
    end
    tempStamp_Unit2 = SpikeStampS1_Unit2{iTrial};
    tempStamp_Unit2(tempStamp_Unit2<=BeforeDelay+BinID-1 | tempStamp_Unit2>BeforeDelay+BinID) = [];
    SpikeNum_S1 = [SpikeNum_S1; min(numel(tempStamp_Unit1),numel(tempStamp_Unit2))];
    tempFcspStamp_Unit2(tempFcspStamp_Unit2<=BeforeDelay+BinID-1 | tempFcspStamp_Unit2>BeforeDelay+BinID) = [];
    FcspNum_S1 = [FcspNum_S1; min(numel(tempFcspStamp_Unit1),numel(tempFcspStamp_Unit2))];
end
NonFcspNum_S1 = SpikeNum_S1 - FcspNum_S1;

%% FCSP number during delay of S2 trials
SpikeNum_S2 = []; FcspNum_S2 = [];
for iTrial = 1:length(SpikeStampS2_Unit1)
    % FCSP number of leading neuron
    tempFcspStamp_Unit1 = [];
    for iSpike = 1:length(SpikeStampS2_Unit1{iTrial})
        if ~isempty(find(SpikeStampS2_Unit2{iTrial}(:,1)-SpikeStampS2_Unit1{iTrial}(iSpike)<=0.01 & SpikeStampS2_Unit2{iTrial}(:,1)-SpikeStampS2_Unit1{iTrial}(iSpike)>0.002))
            tempFcspStamp_Unit1 = [tempFcspStamp_Unit1 SpikeStampS2_Unit1{iTrial}(iSpike)];
        end
    end
    tempStamp_Unit1 = SpikeStampS2_Unit1{iTrial};
    tempStamp_Unit1(tempStamp_Unit1<=BeforeDelay+BinID-1 | tempStamp_Unit1>BeforeDelay+BinID) = [];
    tempFcspStamp_Unit1(tempFcspStamp_Unit1<BeforeDelay+BinID-1 | tempFcspStamp_Unit1>=BeforeDelay+BinID) = [];
    % FCSP number of following neuron
    tempFcspStamp_Unit2 = [];
    for iSpike = 1:length(SpikeStampS2_Unit2{iTrial})
        if ~isempty(find(SpikeStampS2_Unit2{iTrial}(iSpike)-SpikeStampS2_Unit1{iTrial}(:,1)<=0.01 & SpikeStampS2_Unit2{iTrial}(iSpike)-SpikeStampS2_Unit1{iTrial}(:,1)>0.002))
            tempFcspStamp_Unit2 = [tempFcspStamp_Unit2 SpikeStampS2_Unit2{iTrial}(iSpike)];
        end
    end
    tempStamp_Unit2 = SpikeStampS2_Unit2{iTrial};
    tempStamp_Unit2(tempStamp_Unit2<=BeforeDelay+BinID-1 | tempStamp_Unit2>BeforeDelay+BinID) = [];
    SpikeNum_S2 = [SpikeNum_S2; min(numel(tempStamp_Unit1),numel(tempStamp_Unit2))];
    tempFcspStamp_Unit2(tempFcspStamp_Unit2<=BeforeDelay+BinID-1 | tempFcspStamp_Unit2>BeforeDelay+BinID) = [];
    FcspNum_S2 = [FcspNum_S2; min(numel(tempFcspStamp_Unit1),numel(tempFcspStamp_Unit2))];
end
NonFcspNum_S2 = SpikeNum_S2 - FcspNum_S2;

%% Decison variable
[DV_S1, DV_S2] = DecisionVariableCalculation(FcspNum_S1,FcspNum_S2);

%% TPR and FPR
[TPR, FPR] = TprFprCalculation(DV_S1,DV_S2);

%% auROC of FCSP events
AUC_Fcsp = AucAnalysis(FPR,TPR);

%% Resample non-FCSP events amd calculate AUC
AUC_NonFcsp = zeros(SamplingTimes,1);
for iRun = 1:SamplingTimes
    f(iRun) = parfeval(@SampleNonFcspToCalculateAUC,1,FcspNum_S1,FcspNum_S2,NonFcspNum_S1,NonFcspNum_S2);
end
for iRun = 1:SamplingTimes
    [~,tempAUC_NonFcsp] = fetchNext(f);
    AUC_NonFcsp(iRun,1) = tempAUC_NonFcsp;
end
AUC_NonFcsp = mean(AUC_NonFcsp);



