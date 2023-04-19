function [AUC,ShuffledAUC,IsSigAUC] = CalculateIndividualFcAUC(BinID,SpikeStamp1,SpikeStamp2,MarkerS1_Unit1,MarkerS2_Unit1,MarkerS1_Unit2,MarkerS2_Unit2,ShuffleTimes,WhetherShuffle)

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
FcspNum_S1 = [];
for iTrial = 1:length(SpikeStampS1_Unit1)
    % FCSP number of leading neuron
    tempFcspStamp_Unit1 = [];
    for iSpike = 1:length(SpikeStampS1_Unit1{iTrial})
        if ~isempty(find(SpikeStampS1_Unit2{iTrial}(:,1)-SpikeStampS1_Unit1{iTrial}(iSpike)<=0.01 & SpikeStampS1_Unit2{iTrial}(:,1)-SpikeStampS1_Unit1{iTrial}(iSpike)>0.002))
            tempFcspStamp_Unit1 = [tempFcspStamp_Unit1 SpikeStampS1_Unit1{iTrial}(iSpike)];
        end
    end
    tempFcspStamp_Unit1(tempFcspStamp_Unit1<BeforeDelay+BinID-1 | tempFcspStamp_Unit1>=BeforeDelay+BinID) = [];
    % FCSP number of following neuron
    tempFcspStamp_Unit2 = [];
    for iSpike = 1:length(SpikeStampS1_Unit2{iTrial})
        if ~isempty(find(SpikeStampS1_Unit2{iTrial}(iSpike)-SpikeStampS1_Unit1{iTrial}(:,1)<=0.01 & SpikeStampS1_Unit2{iTrial}(iSpike)-SpikeStampS1_Unit1{iTrial}(:,1)>0.002))
            tempFcspStamp_Unit2 = [tempFcspStamp_Unit2 SpikeStampS1_Unit2{iTrial}(iSpike)];
        end
    end
    tempFcspStamp_Unit2(tempFcspStamp_Unit2<=BeforeDelay+BinID-1 | tempFcspStamp_Unit2>BeforeDelay+BinID) = [];
    tempFcspNum = min(numel(tempFcspStamp_Unit1),numel(tempFcspStamp_Unit2));
    FcspNum_S1 = [FcspNum_S1; tempFcspNum];
end

%% FCSP number during delay of S2 trials
FcspNum_S2 = [];
for iTrial = 1:length(SpikeStampS2_Unit1)
    % FCSP number of leading neuron
    tempFcspStamp_Unit1 = [];
    for iSpike = 1:length(SpikeStampS2_Unit1{iTrial})
        if ~isempty(find(SpikeStampS2_Unit2{iTrial}(:,1)-SpikeStampS2_Unit1{iTrial}(iSpike)<=0.01 & SpikeStampS2_Unit2{iTrial}(:,1)-SpikeStampS2_Unit1{iTrial}(iSpike)>0.002))
            tempFcspStamp_Unit1 = [tempFcspStamp_Unit1 SpikeStampS2_Unit1{iTrial}(iSpike)];
        end
    end
    tempFcspStamp_Unit1(tempFcspStamp_Unit1<BeforeDelay+BinID-1 | tempFcspStamp_Unit1>=BeforeDelay+BinID) = [];
    % FCSP number of following neuron
    tempFcspStamp_Unit2 = [];
    for iSpike = 1:length(SpikeStampS2_Unit2{iTrial})
        if ~isempty(find(SpikeStampS2_Unit2{iTrial}(iSpike)-SpikeStampS2_Unit1{iTrial}(:,1)<=0.01 & SpikeStampS2_Unit2{iTrial}(iSpike)-SpikeStampS2_Unit1{iTrial}(:,1)>0.002))
            tempFcspStamp_Unit2 = [tempFcspStamp_Unit2 SpikeStampS2_Unit2{iTrial}(iSpike)];
        end
    end
    tempFcspStamp_Unit2(tempFcspStamp_Unit2<=BeforeDelay+BinID-1 | tempFcspStamp_Unit2>BeforeDelay+BinID) = [];
    tempFcspNum = min(numel(tempFcspStamp_Unit1),numel(tempFcspStamp_Unit2));
    FcspNum_S2 = [FcspNum_S2; tempFcspNum];
end

%% Decison variable
[DV_S1, DV_S2] = DecisionVariableCalculation(FcspNum_S1,FcspNum_S2);

%% TPR and FPR
[TPR, FPR] = TprFprCalculation(DV_S1,DV_S2);

%% auROC
AUC = AucAnalysis(FPR,TPR);

%% Shuffle trial identity, calculate AUC, and perform permutation test
if WhetherShuffle == 1
    ShuffledAUC = zeros(ShuffleTimes,1);
    TriNum_S1 = length(FcspNum_S1);
    for iShuffle = 1:ShuffleTimes
        f(iShuffle) = parfeval(@ShuffledAucCalculation,1,horzcat({FcspNum_S1},{FcspNum_S2}),TriNum_S1,1,1,1);
    end
    for iShuffle = 1:ShuffleTimes
        [~,tempShuffledAUC] = fetchNext(f);
        ShuffledAUC(iShuffle,1) = tempShuffledAUC;
    end
    IsSigAUC = SigLevelAfterPermutation(AUC,ShuffledAUC);
else
    ShuffledAUC = [];
    IsSigAUC = [];
end



