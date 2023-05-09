%% Integrate sec-based FR and AUC results of all neurons.

clear; clc; close all;

%% Assignment
Group = 'CtrlGroup';

%% Integrate FR and AUC
load(sprintf('AllNeuronsFRinGoandNogoTrials_%s',Group));
load(sprintf('AllNeuronsAUCValue_%s',Group));
UnitsFR = cellfun(@(x,y) vertcat(x,y),MeanSortedUnitsFR(:,1),MeanSortedUnitsFR(:,2),'UniformOutput',0);
UnitsFRandAUC = cell(0);
for iUnit = 1:length(UnitsFR)
    UnitsFRandAUC{iUnit,1} = vertcat(UnitsFR{iUnit}(:,1:6),AllUnitsAUC(iUnit,1:6),IsSigAUC(iUnit,1:6));
end
save(sprintf('All Units FR and AUC_%s',Group),'UnitsFRandAUC','AllUnitsShuffledAUC','ID_mPFC','ID_aAIC','-v7.3');


