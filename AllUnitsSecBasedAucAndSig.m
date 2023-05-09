%% AUC of all neurons

clear; clc; close all;
tic

%% Assignment
Group = 'CtrlGroup'; 
IsDVOrFRROC = 1; % 1: decision variable; 2: FR
IsPerformStatTest = 1;
ShuffledTimes = 1000;
WorkerNum = 10;
TimeGain = 10;
if IsPerformStatTest == 1
    poolobj = gcp('nocreate'); 
    if isempty(poolobj)
        myCluster = parcluster('local'); myCluster.NumWorkers = WorkerNum; parpool(myCluster,WorkerNum)
    end
    AllUnitsShuffledAUC = cell(0);
    IsSigAUC = [];
end

%% FR of each trial
load(['UnitsSummary_' Group '.mat']);
% ID of mPFC neurons
ID_mPFC = [];
for iUnit = 1:size(UnitsSummary.mPFC,1)
    ID_mPFC = [ID_mPFC UnitsSummary.mPFC{iUnit,1}(1)];
end
% ID of aAIC neurons
ID_aAIC = [];
for iUnit = 1:size(UnitsSummary.aAIC,1)
    ID_aAIC = [ID_aAIC UnitsSummary.aAIC{iUnit,1}(1)];
end
ID = [ID_mPFC'; ID_aAIC'];
ID = num2cell(ID);
% FR in Go and NoGo trilas of mPFC neurons
load(sprintf('mPFCFRinS1S2_%s',Group)); 
UnitsFRinS1S2_mPFC = vertcat(TargetBrainUnitsFRinS1,TargetBrainUnitsFRinS2); 
UnitsFRinS1S2_mPFC = UnitsFRinS1S2_mPFC';
% FR in Go and NoGo trilas of aAIC neurons
load(sprintf('aAICFRinS1S2_%s',Group));
UnitsFRinS1S2_aAIC = vertcat(TargetBrainUnitsFRinS1,TargetBrainUnitsFRinS2); 
UnitsFRinS1S2_aAIC = UnitsFRinS1S2_aAIC';
% sort FR with neurons ID
UnitsFRinS1S2 = vertcat(UnitsFRinS1S2_mPFC,UnitsFRinS1S2_aAIC);
UnitsFRinS1S2 = horzcat(ID,UnitsFRinS1S2);
SortedUnitsFRinS1S2 = sortrows(UnitsFRinS1S2,1); 
SortedUnitsFRinS1S2(:,1) = [];
% smooth FR
SmoothSortedUnitsFR = cell(size(SortedUnitsFRinS1S2));
windowsize = 3; % 300-msec analysis window
sliding = 1; % 100-msec sliding window
for iUnit = 1:size(SortedUnitsFRinS1S2,1)
    for iTrialType = 1:2
        for iTrial = 1:size(SortedUnitsFRinS1S2{iUnit,iTrialType},1)
            for iBin = 1:sliding:size(SortedUnitsFRinS1S2{iUnit,iTrialType},2)-(windowsize-1)
                SmoothSortedUnitsFR{iUnit,iTrialType}(iTrial,iBin) = sum(SortedUnitsFRinS1S2{iUnit,iTrialType}(iTrial,iBin:iBin+windowsize-1))/windowsize;
            end
        end
        SmoothSortedUnitsFR{iUnit,iTrialType} = SmoothSortedUnitsFR{iUnit,iTrialType}(:,11:end);
    end
end

%% auROC value of all neurons
AllUnitsAUC = [];
for iUnit = 1:size(SmoothSortedUnitsFR,1)
    disp(['Current Unit ID = ' num2str(iUnit) 'of total' num2str(size(SmoothSortedUnitsFR,1)) 'neurons']);
    % step 1: second-based FR in S1 and S2 trials
    for iSecond = 1:floor(size(SmoothSortedUnitsFR{1,1},2)/TimeGain)
        FRinBin_S1 = mean(SmoothSortedUnitsFR{iUnit,1}(:,1+(iSecond-1)*TimeGain:iSecond*TimeGain),2);
        FRinBin_S2 = mean(SmoothSortedUnitsFR{iUnit,2}(:,1+(iSecond-1)*TimeGain:iSecond*TimeGain),2);
        if IsDVOrFRROC == 1
            % step 2: DV of each unit in S1 and S2 trials
            [S1DV,S2DV] = DecisionVariableCalculation(FRinBin_S1,FRinBin_S2);
            % step 3: FPR and TPR of each unit
            [TPR,FPR] = TprFprCalculation(S1DV,S2DV);
        else
            % step 3: FPR and TPR of each unit
            [TPR,FPR] = TprFprCalculation(FRinBin_S1',FRinBin_S2');
        end
        % step 4: auROC value of each unit in each bin
        AUC = AucAnalysis(FPR,TPR);
        AllUnitsAUC(iUnit,iSecond) = AUC;
        % step 5: perform permutation test by shuffling trials
        if IsPerformStatTest == 1
            ShuffledAUC = zeros(ShuffledTimes,1);
            for iShuffle = 1:ShuffledTimes
                f(iShuffle) = parfeval(@ShuffledAucCalculation,1,SmoothSortedUnitsFR(iUnit,:),numel(FRinBin_S1),iSecond,TimeGain,IsDVOrFRROC);
            end
            for iShuffle = 1:ShuffledTimes
                [~,tempShuffledAUC] = fetchNext(f);
                ShuffledAUC(iShuffle,1) = tempShuffledAUC;
            end
            AllUnitsShuffledAUC{iUnit,iSecond} = ShuffledAUC;
            [Sig] = SigLevelAfterPermutation(AUC,ShuffledAUC);
            IsSigAUC(iUnit,iSecond) = Sig;
        end
    end
end
save(['AllNeuronsAUCValue_' Group],'AllUnitsAUC','AllUnitsShuffledAUC','IsSigAUC','ID_mPFC','ID_aAIC');
toc





