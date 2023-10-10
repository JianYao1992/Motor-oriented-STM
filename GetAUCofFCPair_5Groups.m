%% This code calculates ratio of pairs with selective FC in FC in non-memory, incongruent inactive,incongruent active, congruent inactive, 
%% and congruent active pairs.
clear; clc; close all;

%% Assignment
Group = 'ChR2Group'; 
Target = 'mPFC-aAIC'; 
PairType = {'congruent active','congruent inactive','incongruent active','incongruent inactive','non-memory'};

%% Load result of FCSP coding
load(sprintf('FC Rate and Coding_%s_%s',Target,Group));

%% Construct FR selectivity for each pair of neurons
load(sprintf('UnitsSummary_%s',Group));
ID_mPFC = [];
ID_aAIC = [];
for iUnit = 1:size(UnitsSummary.mPFC,1)
    ID_mPFC = [ID_mPFC UnitsSummary.mPFC{iUnit,1}(1)];
end
for iUnit = 1:size(UnitsSummary.aAIC,1)
    ID_aAIC = [ID_aAIC UnitsSummary.aAIC{iUnit,1}(1)];
end
ID = [ID_mPFC';ID_aAIC'];
ID = num2cell(ID);
AllUnitsInfo = [UnitsSummary.mPFC; UnitsSummary.aAIC];
AllUnitsInfo = [ID AllUnitsInfo]; 
SortedAllUnitsInfo = sortrows(AllUnitsInfo,1);
clear UnitsSummary AllUnitsInfo UnitTrialMark SpikeTimeIndividualTrial
AllUnitsSelectivity = cell2mat(SortedAllUnitsInfo(:,4));

%% Load information of FC pairs with functional coupling in each delay bin.
PairSelectivity = cell(0);
for bin = 1:length(UniqueFCPairAUC)
    tempID = [];
    for iPair = 1:size(UniqueFCPairAUC{bin},1)
        PairSelectivity{bin}(iPair,1:7) = AllUnitsSelectivity(UniqueFCPairAUC{bin}(iPair,1),1:7);
        PairSelectivity{bin}(iPair,8:14) = AllUnitsSelectivity(UniqueFCPairAUC{bin}(iPair,2),1:7);
    end
    % Select pairs of neurons within non-memory or memory neurons.
    tempID = find((all(PairSelectivity{bin}(:,3:6)==0,2) & all(PairSelectivity{bin}(:,10:13)==0,2))|(~((any(PairSelectivity{bin}(:,3:6)==1,2) & any(PairSelectivity{bin}(:,3:6)==2,2))|(any(PairSelectivity{bin}(:,10:13)==1,2) & any(PairSelectivity{bin}(:,10:13)==2,2))) & max(PairSelectivity{bin}(:,3:6),[],2)>0 & max(PairSelectivity{bin}(:,10:13),[],2)>0));
    UniqueFCPairAUC{bin} = UniqueFCPairAUC{bin}(tempID,:);
    PairSelectivity{bin} = PairSelectivity{bin}(tempID,:);
end

%% Classify FC pairs into different types
TargetPairAUCofFC = cell(0);
for i = 1:length(PairType)
    for bin = 1:4
        tempID = [];
        switch PairType{i}
            case 'non-memory'
                tempID = find(all(PairSelectivity{bin}(:,3:6)==0,2));
            case 'incongruent inactive'
                tempID = find(~all(PairSelectivity{bin}(:,3:6)==0,2) & PairSelectivity{bin}(:,2+bin).*PairSelectivity{bin}(:,9+bin)==0 & max(PairSelectivity{bin}(:,3:6),[],2) ~= max(PairSelectivity{bin}(:,10:13),[],2));
            case 'incongruent active'
                tempID = find(~all(PairSelectivity{bin}(:,3:6)==0,2) & PairSelectivity{bin}(:,2+bin).*PairSelectivity{bin}(:,9+bin)>0 & PairSelectivity{bin}(:,9+bin)>0 & PairSelectivity{bin}(:,2+bin)~=PairSelectivity{bin}(:,9+bin));
            case 'congruent inactive'
                tempID = find(~all(PairSelectivity{bin}(:,3:6)==0,2) & PairSelectivity{bin}(:,2+bin).*PairSelectivity{bin}(:,9+bin)==0 & max(PairSelectivity{bin}(:,3:6),[],2)==max(PairSelectivity{bin}(:,10:13),[],2));
            case 'congruent active'
                tempID = find(PairSelectivity{bin}(:,2+bin)>0 & PairSelectivity{bin}(:,9+bin)>0 & max(PairSelectivity{bin}(:,3:6),[],2)==max(PairSelectivity{bin}(:,10:13),[],2));
        end
        TargetPairAUCofFC{1,i}{1,bin} = UniqueFCPairAUC{bin}(tempID,:);
    end
end
% RatioofFCPair = [];
% for iType = 1:length(TargetPairAUCofFC)
%     temp = vertcat(TargetPairAUCofFC{iType,1}{:});
%     RatioofFCPair(1,iType) = nnz(temp(:,3)>0.5 & temp(:,4)==1);
%     RatioofFCPair(2,iType) = size(temp,1);
%     RatioofFCPair(3,iType) = nnz(temp(:,3)>0.5 & temp(:,4)==1)/size(temp,1);
% end

%% Save file
save(sprintf('Functional Coupling Coding_%dpairs_%s_%s',numel(PairType),Target,Group),'TargetPairAUCofFC','-v7.3');
