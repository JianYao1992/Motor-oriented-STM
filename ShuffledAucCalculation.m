function ShuffledAUC = ShuffledAucCalculation(FRinTwoKindsOfTrials,TrialNuminFirstKindOfTrial,SpecificWindow,TimeGain,IsDVorFrROC)

%% Step 1: Get the mean FR of all trials in two conditions after shuffling
MeanFR_S1 = mean(FRinTwoKindsOfTrials{1}(:,(1+TimeGain*(SpecificWindow-1)):TimeGain*SpecificWindow),2);
MeanFR_S2 = mean(FRinTwoKindsOfTrials{2}(:,(1+TimeGain*(SpecificWindow-1)):TimeGain*SpecificWindow),2);
AllTrialMeanFR = [MeanFR_S1; MeanFR_S2];
ShuffledOrder = randperm(size(AllTrialMeanFR,1));
ShuffledMeanFR_S1 = AllTrialMeanFR(ShuffledOrder(1:TrialNuminFirstKindOfTrial),:);
ShuffledMeanFR_S2 = AllTrialMeanFR(ShuffledOrder((1+TrialNuminFirstKindOfTrial):end),:);
if IsDVorFrROC==1
    %% Step 2: Calculate the DV between S1 and S2
    [ShuffledDV_S1,ShuffledDV_S2] = DecisionVariableCalculation(ShuffledMeanFR_S1,ShuffledMeanFR_S2);
    %% Step 3: Calculate the TPR and FPR of unit
    [TPR,FPR] = TprFprCalculation(ShuffledDV_S1,ShuffledDV_S2);
else
    %% Step 3: Calculate the TPR and FPR of unit
    [TPR,FPR] = TprFprCalculation(ShuffledMeanFR_S1',ShuffledMeanFR_S2');
end
%% Step 4: Calculate the Area under curve (AUC) of each unit
ShuffledAUC = AucAnalysis(FPR,TPR);
