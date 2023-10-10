%% Coding ability of FCSP events for Go and NoGo trials during delay period

clear; clc; close all;

%% Assignment
Group = 'ChR2Group';
Target = 'Within aAIC';
IsConsiderNonsimulData = 0;
PairsType = {'Con Act','Con Inact','Incong Act','Incon Inact','Non-memory'};
BinSize = 1;
BinRange = 1:BinSize:4;
BinNum = numel(BinRange);
if strcmp(Target,'Within mPFC')
    tarreg1 = 1;
    tarreg2 = 1;
elseif strcmp(Target,'Within aAIC')
    tarreg1 = 2;
    tarreg2 = 2;
elseif strcmp(Target,'mPFC-aAIC')
    tarreg1 = 1;
    tarreg2 = 2;
end

%% Nonsimultaneously recorded neuronal pairs in target area
% Judge whether the real value is significantly higher than unsimultaneously recorded result.
if IsConsiderNonsimulData == 1
    load(sprintf('Unsimultaneously recorded neuron pairs_%s',Group)); % here, the pairs in UnSimulPairs have excluded the simultaneously recorded and switched neuron pairs.
    % pairs of nonsimultaneously recorded neurons in target area
    UnSimulPairs = UnSimulPairs(UnSimulPairs(:,3)==tarreg1 & UnSimulPairs(:,14)==tarreg2,:);
    TargetUnSimulPairs = cell(0);
    for i = 1:numel(PairsType)
        if strcmp(PairsType{i},'Con Act')
            for iBin = 1:BinNum
                TargetUnSimulPairs{i}{iBin} = UnSimulPairs(UnSimulPairs(:,6+iBin)==UnSimulPairs(:,17+iBin) & UnSimulPairs(:,6+iBin)>0,:);
            end
        elseif strcmp(PairsType{i},'Con Inact')
            for iBin = 1:BinNum
                TargetUnSimulPairs{i}{iBin} = UnSimulPairs(UnSimulPairs(:,6+iBin).*UnSimulPairs(:,17+iBin)==0 & max(UnSimulPairs(:,7:10),[],2)==max(UnSimulPairs(:,18:21),[],2) & max(UnSimulPairs(:,7:10),[],2)>0,:);
            end
        elseif strcmp(PairsType{i},'Incon Act')
            for iBin = 1:BinNum
                TargetUnSimulPairs{i}{iBin} = UnSimulPairs(UnSimulPairs(:,6+iBin)>0 & UnSimulPairs(:,17+iBin)>0 & UnSimulPairs(:,6+iBin)~=UnSimulPairs(:,17+iBin),:);
            end
        elseif strcmp(PairsType{i},'Incon Inact')
            for iBin = 1:BinNum
                TargetUnSimulPairs{i}{iBin} = UnSimulPairs(UnSimulPairs(:,6+iBin).*UnSimulPairs(:,17+iBin)==0 & max(UnSimulPairs(:,7:10),[],2)>0 & max(UnSimulPairs(:,18:21),[],2)>0 & max(UnSimulPairs(:,7:10),[],2)~=max(UnSimulPairs(:,18:21),[],2),:);
            end
        elseif strcmp(PairsType{i},'Non-memory')
            for iBin = 1:BinNum
                TargetUnSimulPairs{i}{iBin} = UnSimulPairs(all(UnSimulPairs(:,7:10)==0,2) & all(UnSimulPairs(:,18:21)==0,2),:);
            end
        end
    end
end

%% Parallel pool
WorkerNum = 10;
NullDistribution = 10;
SamplingRun = 10;
poolobj = gcp('nocreate');
if isempty(poolobj)
    myCluster = parcluster('local'); myCluster.NumWorkers = WorkerNum; parpool(myCluster,WorkerNum);
end

%% Spike rasters of each trial
load(['UnitsSummary_' Group '.mat']);
% ID of mPFC and aAIC neurons
ID_mPFC = [];
for iUnit = 1:size(UnitsSummary.mPFC,1)
    ID_mPFC = [ID_mPFC UnitsSummary.mPFC{iUnit,1}(1)];
end
ID_aAIC = [];
for iUnit = 1:size(UnitsSummary.aAIC,1)
    ID_aAIC = [ID_aAIC UnitsSummary.aAIC{iUnit,1}(1)];
end
% spike rasters
load(sprintf('NeuronOriginandSpikeRasterInformationfor%s.mat',Group));
ID = [ID_mPFC'; ID_aAIC'];
ID = num2cell(ID);
UnitTrialMark = horzcat(ID, [mPFCUnitTrialMark; aAICUnitTrialMark]); % sort trial mark
SortedUnitTriMark = sortrows(UnitTrialMark,1);
SpikeTimeIndivTrial = horzcat(ID, [mPFCSpikeTimeinIndividualTrial; aAICSpikeTimeinIndividualTrial]); % sort spike raster
SortedSpikeTimeIndivTrial = sortrows(SpikeTimeIndivTrial,1);
clear UnitsSummary AllUnitsInfo UnitTrialMark SpikeTimeIndivTrial

%% Load FR of all neurons
load(sprintf('All Units FR and AUC_%s',Group));

%% FC neuronal pairs of each bin of delay
AllTypesFCPairID = cell(0);
AllTypesAUC = cell(0);
AllNonsimulAUC = cell(0);
DiffinFR = cell(1,numel(PairsType));
for type = 1:numel(PairsType)
    disp(['Analyze ' PairsType{type} ' in ' Target]);
    FCPairsID = cell(1,BinNum);
    DiffinFR{type} = cell(1,numel(BinRange));
    for bin = 1:BinNum
        FcPairs_Go = [];
        FcPairs_NoGo = [];
        load(['FCofNeurons_PropertyPopulation_1msbin_' num2str(BinRange(bin)) '_' num2str(BinRange(bin)+BinSize) '_GoNoGo_' Group '.mat']);
        
        %% FC in Go trials
        PairsMarkerOfSwitUnit_Go = [];
        for iPair = 1:size(conn_chain_go,1)
            Pref_Unit1 = unique(pref_chain_go(iPair,3:6));
            Pref_Unit1(Pref_Unit1==0) = [];
            Pref_Unit2 = unique(pref_chain_go(iPair,10:13));
            Pref_Unit2(Pref_Unit2==0) = [];
            if (~isempty(Pref_Unit1) && max(Pref_Unit1) ~= min(Pref_Unit1)) || (~isempty(Pref_Unit2) && max(Pref_Unit2) ~= min(Pref_Unit2))
                PairsMarkerOfSwitUnit_Go(end+1,1) = iPair;
            end
            % left is the leading neuron, right is the following neuron
            if conn_chain_go(iPair,3) < 0
                temp = conn_chain_go(iPair,1); conn_chain_go(iPair,1) = conn_chain_go(iPair,2); conn_chain_go(iPair,2) = temp;
                conn_chain_go(iPair,3) = -1*conn_chain_go(iPair,3);
                temp = pref_chain_go(iPair,1:7); pref_chain_go(iPair,1:7) = pref_chain_go(iPair,8:14); pref_chain_go(iPair,8:14) = temp;
                temp = reg_chain_go(iPair,1); reg_chain_go(iPair,1) = reg_chain_go(iPair,2); reg_chain_go(iPair,2) = temp;
                reg_chain_go(iPair,3) = -1*reg_chain_go(iPair,3);
            end
        end
        % remove FC pairs with at least 1 switched neuron
        conn_chain_go(PairsMarkerOfSwitUnit_Go,:) = [];
        pref_chain_go(PairsMarkerOfSwitUnit_Go,:) = [];
        reg_chain_go(PairsMarkerOfSwitUnit_Go,:) = [];
        
        %% FC in NoGo trials
        PairsMarkerOfSwitUnit_NoGo = [];
        for iPair = 1:size(conn_chain_nogo,1)
            Pref_Unit1 = unique(pref_chain_nogo(iPair,3:6)); Pref_Unit1(Pref_Unit1==0) = [];
            Pref_Unit2 = unique(pref_chain_nogo(iPair,10:13)); Pref_Unit2(Pref_Unit2==0) = [];
            if (~isempty(Pref_Unit1) && max(Pref_Unit1) ~= min(Pref_Unit1)) || (~isempty(Pref_Unit2) && max(Pref_Unit2) ~= min(Pref_Unit2))
                PairsMarkerOfSwitUnit_NoGo(end+1,1) = iPair;
            end
            if conn_chain_nogo(iPair,3) < 0
                temp = conn_chain_nogo(iPair,1); conn_chain_nogo(iPair,1) = conn_chain_nogo(iPair,2); conn_chain_nogo(iPair,2) = temp;
                conn_chain_nogo(iPair,3) = -1*conn_chain_nogo(iPair,3);
                temp = pref_chain_nogo(iPair,1:7); pref_chain_nogo(iPair,1:7) = pref_chain_nogo(iPair,8:14); pref_chain_nogo(iPair,8:14) = temp;
                temp = reg_chain_nogo(iPair,1); reg_chain_nogo(iPair,1) = reg_chain_nogo(iPair,2); reg_chain_nogo(iPair,2) = temp;
                reg_chain_nogo(iPair,3) = -1*reg_chain_nogo(iPair,3);
            end
        end
        % remove FC pairs with at least 1 switched neuron
        conn_chain_nogo(PairsMarkerOfSwitUnit_NoGo,:) = [];
        pref_chain_nogo(PairsMarkerOfSwitUnit_NoGo,:) = [];
        reg_chain_nogo(PairsMarkerOfSwitUnit_NoGo,:) = [];
        
        %% Target FC pairs
        switch PairsType{type}
            case 'Non-memory'
                FcPairs_Go = conn_chain_go(ismember(reg_chain_go(:,1),tarreg1) & ismember(reg_chain_go(:,2),tarreg2) & all(pref_chain_go(:,3:6)==0,2) & all(pref_chain_go(:,10:13)==0,2),:);
                FcPairs_NoGo = conn_chain_nogo(ismember(reg_chain_nogo(:,1),tarreg1) & ismember(reg_chain_nogo(:,2),tarreg2) & all(pref_chain_nogo(:,3:6)==0,2) & all(pref_chain_nogo(:,10:13)==0,2),:);
            case 'Con Act'
                FcPairs_Go = conn_chain_go(ismember(reg_chain_go(:,1),tarreg1) & ismember(reg_chain_go(:,2),tarreg2) & pref_chain_go(:,2+bin)==pref_chain_go(:,9+bin) & pref_chain_go(:,2+bin)>0,:);
                FcPairs_NoGo = conn_chain_nogo(ismember(reg_chain_nogo(:,1),tarreg1) & ismember(reg_chain_nogo(:,2),tarreg2) & pref_chain_nogo(:,2+bin)==pref_chain_nogo(:,9+bin) & pref_chain_nogo(:,2+bin)>0,:);
            case 'Con Inact'
                FcPairs_Go = conn_chain_go(ismember(reg_chain_go(:,1),tarreg1) & ismember(reg_chain_go(:,2),tarreg2) & max(pref_chain_go(:,3:6),[],2) == max(pref_chain_go(:,10:13),[],2) & max(pref_chain_go(:,3:6),[],2)>0 & any(pref_chain_go(:,[bin+2 bin+9])<1,2),:);
                FcPairs_NoGo = conn_chain_nogo(ismember(reg_chain_nogo(:,1),tarreg1) & ismember(reg_chain_nogo(:,2),tarreg2) & max(pref_chain_nogo(:,3:6),[],2) == max(pref_chain_nogo(:,10:13),[],2) & max(pref_chain_nogo(:,3:6),[],2)>0 & any(pref_chain_nogo(:,[bin+2 bin+9])<1,2),:);
            case 'Incong Act'
                FcPairs_Go = conn_chain_go(ismember(reg_chain_go(:,1),tarreg1) & ismember(reg_chain_go(:,2),tarreg2) & pref_chain_go(:,2+bin)>0 & pref_chain_go(:,9+bin)>0 & pref_chain_go(:,2+bin)~=pref_chain_go(:,9+bin),:);
                FcPairs_NoGo = conn_chain_nogo(ismember(reg_chain_nogo(:,1),tarreg1) & ismember(reg_chain_nogo(:,2),tarreg2) & pref_chain_nogo(:,2+bin)>0 & pref_chain_nogo(:,9+bin)>0 & pref_chain_nogo(:,2+bin)~=pref_chain_nogo(:,9+bin),:);
            case 'Incon Inact'
                FcPairs_Go = conn_chain_go(ismember(reg_chain_go(:,1),tarreg1) & ismember(reg_chain_go(:,2),tarreg2) & max(pref_chain_go(:,3:6),[],2)~=max(pref_chain_go(:,10:13),[],2) & max(pref_chain_go(:,3:6),[],2)>0 & max(pref_chain_go(:,10:13),[],2)>0 & any(pref_chain_go(:,[bin+2 bin+9])<1,2),:);
                FcPairs_NoGo = conn_chain_nogo(ismember(reg_chain_nogo(:,1),tarreg1) & ismember(reg_chain_nogo(:,2),tarreg2) & max(pref_chain_nogo(:,3:6),[],2)~=max(pref_chain_nogo(:,10:13),[],2) & max(pref_chain_nogo(:,3:6),[],2)>0 & max(pref_chain_nogo(:,10:13),[],2)>0 & any(pref_chain_nogo(:,[bin+2 bin+9])<1,2),:);
        end
        % unique FCs
        if ~isempty(FcPairs_Go) && isempty(FcPairs_NoGo)
            FCPairsID{bin} = FcPairs_Go(:,1:2);
        elseif isempty(FcPairs_Go) && ~isempty(FcPairs_NoGo)
            FCPairsID{bin} = FcPairs_NoGo(:,1:2);
        elseif ~isempty(FcPairs_Go) && ~isempty(FcPairs_NoGo)
            FCPairsID{bin} = unique(vertcat(FcPairs_Go(:,1:2),FcPairs_NoGo(:,1:2)),'rows');
        end
        
        %% FR difference between Go and NoGo trials of each FC neuronal pair
        if ~isempty(FcPairs_Go) && ~isempty(FcPairs_NoGo)
            DiffinFR{type}{bin} = CalculateFcPairsFrDiff(bin,FcPairs_Go(:,1:2),FcPairs_NoGo(:,1:2),UnitsFRandAUC);
        elseif ~isempty(FcPairs_Go) && isempty(FcPairs_NoGo)
            DiffinFR{type}{bin} = CalculateFcPairsFrDiff(bin,FcPairs_Go(:,1:2),[],UnitsFRandAUC);
        elseif isempty(FcPairs_Go) && ~isempty(FcPairs_NoGo)
            DiffinFR{type}{bin} = CalculateFcPairsFrDiff(bin,[],FcPairs_NoGo(:,1:2),UnitsFRandAUC);
        end
    end
    
    %% auROC of FCSP events for Go and NoGo trials
    [AUC,~] = CalculatePopulationFcAUC(BinRange,FCPairsID,SortedUnitTriMark,SortedSpikeTimeIndivTrial,1000,0);
    AllTypesAUC{type} = AUC;
    AllTypesFCPairID{type} = FCPairsID;
    
    %% Sample pairs of unsimultaneously recorded neurons
    if IsConsiderNonsimulData == 1
        AllNonsimulAUC{type} = [];
        FCPairNum = cellfun(@(x) size(x,1),FCPairsID);
        for iNullDis = 1:NullDistribution
            f(iNullDis) = parfeval(@CalculateUnsimulPairsFunctionalCouplingAUC,1,BinRange,FCPairNum,TargetUnSimulPairs{type},SortedUnitTriMark,SortedSpikeTimeIndivTrial,SamplingRun);
        end
        for iNullDis = 1:NullDistribution
            [~, UnSimulAUC] = fetchNext(f);
            AllNonsimulAUC{type} = [AllNonsimulAUC{type}; UnSimulAUC];
        end
    end    
end

%% Save variables
if IsConsiderNonsimulData == 1
    save(sprintf('Functional Coupling Coding_1msbin_%dpairs_%s_%s',numel(PairsType),Target,Group),'AllTypesAUC','AllTypesFCPairID','AllUnSimulAUC','-v7.3');
else
    save(sprintf('Functional Coupling Coding_1msbin_%dpairs_%s_%s',numel(PairsType),Target,Group),'AllTypesAUC','AllTypesFCPairID','-v7.3');
end
save(sprintf('Functional Coupling Firing Rate Difference_1msbin_%dpairs_%s_%s',numel(PairsType),Target,Group),'DiffinFR','-v7.3');


