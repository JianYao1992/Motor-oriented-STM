%% Second-based FR of all neurons in Go and NoGo trials.

clear; clc; close all;

%% Assignment
Group = 'CtrlGroup'; 
TimeGain = 10;

%% FR of all neurons in each trial
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
% FR in Go and NoGo trials of mPFC neurons
load(sprintf('mPFCFRinS1S2_%s',Group)); 
UnitsFRinS1S2_mPFC = vertcat(TargetBrainUnitsFRinS1,TargetBrainUnitsFRinS2);
UnitsFRinS1S2_mPFC = UnitsFRinS1S2_mPFC';
% FR in Go and NoGo trials of aAIC neurons
load(sprintf('aAICFRinS1S2_%s',Group)); 
UnitsFRinS1S2_aAIC = vertcat(TargetBrainUnitsFRinS1,TargetBrainUnitsFRinS2);
UnitsFRinS1S2_aAIC = UnitsFRinS1S2_aAIC';
% sort FR with neurons ID
UnitsFRinS1S2 = vertcat(UnitsFRinS1S2_mPFC,UnitsFRinS1S2_aAIC); 
UnitsFRinS1S2 = horzcat(ID,UnitsFRinS1S2);
SortedUnitsFRinS1S2 = sortrows(UnitsFRinS1S2,1);
SortedUnitsFRinS1S2(:,1) = [];
% calculate FR of each unit
MeanSortedUnitsFR = cell(size(SortedUnitsFRinS1S2));
for iUnit = 1:size(SortedUnitsFRinS1S2,1)
    disp(fprintf('smoothing_%d_th Unit\n',iUnit));
    for iTrialType = 1:2
        for bin = 1:1:size(SortedUnitsFRinS1S2{1,1},2)/TimeGain
            MeanSortedUnitsFR{iUnit,iTrialType}(1,bin) = mean(mean(SortedUnitsFRinS1S2{iUnit,iTrialType}(:,1+TimeGain*(bin-1):TimeGain*bin)));
        end
        MeanSortedUnitsFR{iUnit,iTrialType} = MeanSortedUnitsFR{iUnit,iTrialType}(:,2:end);
    end
end
save(['AllNeuronsFRinGoandNogoTrials_' Group],'MeanSortedUnitsFR','-v7.3');






