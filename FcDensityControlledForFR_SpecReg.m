%% Region-specific FC density of neuronal pairs controlled for FR

clear; clc; close all;

%% Assignment
Group = 'ChR2Group';
DelayBinNum = 4;
TarFR = [2 13];

%% Region-specific FC of various types of neuronal pairs, controlled for FR
IsSigFC_hit_ConAct = cell(5,numel(TarFR)); % each row represents each region, 1: all regions; 2: mPFC; 3: aAIC; 4: mPFC-aAIC; 5: aAIC-mPFC
IsSigFC_hit_ConInact = cell(5,numel(TarFR));
IsSigFC_hit_InconAct = cell(5,numel(TarFR)); 
IsSigFC_hit_InconInact = cell(5,numel(TarFR)); 
IsSigFC_hit_Nonmem = cell(5,numel(TarFR));
for i = 1:numel(TarFR)
    if i==1
        FRRange = [TarFR(i) TarFR(i)+5];
    elseif i==2
        FRRange = [TarFR(i)-6 TarFR(i)];
    end
    %% FC density of pairs of non-memory neurons
    for bin = 1:DelayBinNum % four 1-sec bins of delay
        load(sprintf('FCofNonmemoryNeurons_PropertyPopulation_1msbin_%d_%d_%s.mat',bin,bin+1,Group));
        load(sprintf('PairsFRModulation_1msbin_nonsel_conn_chain_delay4s_%d_%d_%s.mat',bin,bin+1,Group));
        % mPFC, aAIC, mPFC-aAIC, and aAIC-mPFC
        tempPairsNum_hit_Nonmem = nnz(ismember(Pair_reg(:,1),[1 2]) & ismember(Pair_reg(:,2),[1 2]) & FRPairChains_go(:,1)>=FRRange(1) & FRPairChains_go(:,1)<FRRange(2) & FRPairChains_go(:,2)>=FRRange(1) & FRPairChains_go(:,2)<FRRange(2));
        tempFcPairsNum_hit_Nonmem = nnz(ismember(reg_chain_go(:,1),[1 2]) & ismember(reg_chain_go(:,2),[1 2]) & FR_conn_Chains_go(:,1)>=FRRange(1) & FR_conn_Chains_go(:,1)<FRRange(2) & FR_conn_Chains_go(:,2)>=FRRange(1) & FR_conn_Chains_go(:,2)<FRRange(2)); 
        IsSigFC_hit_Nonmem{1,i}{1,bin} = zeros(tempPairsNum_hit_Nonmem,1);
        if tempPairsNum_hit_Nonmem > 0
            IsSigFC_hit_Nonmem{1,i}{1,bin}(1:tempFcPairsNum_hit_Nonmem) = 1;
        end
        % mPFC
        tempPairsNum_hit_mPFC_Nonmem = nnz(ismember(Pair_reg(:,1),1) & ismember(Pair_reg(:,2),1) & FRPairChains_go(:,1)>=FRRange(1) & FRPairChains_go(:,1)<FRRange(2) & FRPairChains_go(:,2)>=FRRange(1) & FRPairChains_go(:,2)<FRRange(2));
        tempFcPairsNum_hit_mPFC_Nonmem = nnz(ismember(reg_chain_go(:,1),1) & ismember(reg_chain_go(:,2),1) & FR_conn_Chains_go(:,1)>=FRRange(1) & FR_conn_Chains_go(:,1)<FRRange(2) & FR_conn_Chains_go(:,2)>=FRRange(1) & FR_conn_Chains_go(:,2)<FRRange(2));
        IsSigFC_hit_Nonmem{2,i}{1,bin} = zeros(tempPairsNum_hit_mPFC_Nonmem,1);
        if tempPairsNum_hit_mPFC_Nonmem > 0
            IsSigFC_hit_Nonmem{2,i}{1,bin}(1:tempFcPairsNum_hit_mPFC_Nonmem) = 1;
        end
        % aAIC
        tempPairsNum_hit_aAIC_Nonmem = nnz(ismember(Pair_reg(:,1),2) & ismember(Pair_reg(:,2),2) & FRPairChains_go(:,1)>=FRRange(1) & FRPairChains_go(:,1)<FRRange(2) & FRPairChains_go(:,2)>=FRRange(1) & FRPairChains_go(:,2)<FRRange(2));
        tempFcPairsNum_hit_aAIC_Nonmem = nnz(ismember(reg_chain_go(:,1),2) & ismember(reg_chain_go(:,2),2) & FR_conn_Chains_go(:,1)>=FRRange(1) & FR_conn_Chains_go(:,1)<FRRange(2) & FR_conn_Chains_go(:,2)>=FRRange(1) & FR_conn_Chains_go(:,2)<FRRange(2));
        IsSigFC_hit_Nonmem{3,i}{1,bin} = zeros(tempPairsNum_hit_aAIC_Nonmem,1);
        if tempPairsNum_hit_aAIC_Nonmem > 0
            IsSigFC_hit_Nonmem{3,i}{1,bin}(1:tempFcPairsNum_hit_aAIC_Nonmem) = 1;
        end
        % mPFC-aAIC
        tempPairs_hit_mPFCaAIC_Nonmem = Pair_chain(ismember(Pair_reg(:,1),2) & ismember(Pair_reg(:,2),1) & FRPairChains_go(:,1)>=FRRange(1) & FRPairChains_go(:,1)<FRRange(2) & FRPairChains_go(:,2)>=FRRange(1) & FRPairChains_go(:,2)<FRRange(2),:);
        tempFcPairs_hit_mPFCaAIC_Nonmem = conn_chain_go(ismember(reg_chain_go(:,1),2) & ismember(reg_chain_go(:,2),1) & reg_chain_go(:,3)<0 & FR_conn_Chains_go(:,1)>=FRRange(1) & FR_conn_Chains_go(:,1)<FRRange(2) & FR_conn_Chains_go(:,2)>=FRRange(1) & FR_conn_Chains_go(:,2)<FRRange(2),1:2);
        IsSigFC_hit_Nonmem{4,i}{1,bin} = zeros(size(tempPairs_hit_mPFCaAIC_Nonmem,1),1);
        for iFcPair = 1:size(tempFcPairs_hit_mPFCaAIC_Nonmem,1)
            tempUnit1 = tempFcPairs_hit_mPFCaAIC_Nonmem(iFcPair,1);
            tempUnit2 = tempFcPairs_hit_mPFCaAIC_Nonmem(iFcPair,2);
            tempPos = find(tempPairs_hit_mPFCaAIC_Nonmem(:,1)==tempUnit1 & tempPairs_hit_mPFCaAIC_Nonmem(:,2)==tempUnit2);
            if ~isempty(tempPos)
                IsSigFC_hit_Nonmem{4,i}{1,bin}(tempPos) = 1;
            end
        end
        % aAIC-mPFC
        tempPairs_hit_aAICmPFC_Nonmem = Pair_chain(ismember(Pair_reg(:,1),2) & ismember(Pair_reg(:,2),1) & FRPairChains_go(:,1)>=FRRange(1) & FRPairChains_go(:,1)<FRRange(2) & FRPairChains_go(:,2)>=FRRange(1) & FRPairChains_go(:,2)<FRRange(2),:);
        tempFcPairs_hit_aAICmPFC_Nonmem = conn_chain_go(ismember(reg_chain_go(:,1),2) & ismember(reg_chain_go(:,2),1) & reg_chain_go(:,3)>0 & FR_conn_Chains_go(:,1)>=FRRange(1) & FR_conn_Chains_go(:,1)<FRRange(2) & FR_conn_Chains_go(:,2)>=FRRange(1) & FR_conn_Chains_go(:,2)<FRRange(2),1:2);
        IsSigFC_hit_Nonmem{5,i}{1,bin} = zeros(size(tempPairs_hit_aAICmPFC_Nonmem,1),1);
        for iFcPair = 1:size(tempFcPairs_hit_aAICmPFC_Nonmem,1)
            tempUnit1 = tempFcPairs_hit_aAICmPFC_Nonmem(iFcPair,1);
            tempUnit2 = tempFcPairs_hit_aAICmPFC_Nonmem(iFcPair,2);
            tempPos = find(tempPairs_hit_aAICmPFC_Nonmem(:,1)==tempUnit1 & tempPairs_hit_aAICmPFC_Nonmem(:,2)==tempUnit2);
            if ~isempty(tempPos)
                IsSigFC_hit_Nonmem{5,i}{1,bin}(tempPos) = 1;
            end
        end
        clear Pair_reg reg_chain_go  FRPairChains_go  FR_conn_Chains_go 
    end
    
    %% FC density of pairs of memory neurons
    for bin = 1:DelayBinNum % four 1-second bins of delay
        load(sprintf('FCofMemoryNeurons_PropertyPopulation_1msbin_%d_%d_%s.mat',bin,bin+1,Group));
        load(sprintf('PairsFRModulation_1msbin_selec_conn_chain_delay4s_%d_%d_%s.mat',bin,bin+1,Group));
        
       %% Congruent active pairs
        % mPFC, aAIC, mPFC-aAIC, and aAIC-mPFC
        tempPairsNum_hit_ConAct = nnz(ismember(Pair_reg(:,1),[1 2]) & ismember(Pair_reg(:,2),[1 2]) & Pref_pair(:,2+bin)==Pref_pair(:,9+bin) & Pref_pair(:,2+bin)>0 & FRPairChains_go(:,1)>=FRRange(1) & FRPairChains_go(:,1)<FRRange(2) & FRPairChains_go(:,2)>=FRRange(1) & FRPairChains_go(:,2)<FRRange(2));
        tempFcPairsNum_hit_ConAct = nnz(ismember(reg_chain_go(:,1),[1 2]) & ismember(reg_chain_go(:,2),[1 2]) & pref_chain_go(:,2+bin)==pref_chain_go(:,9+bin) & pref_chain_go(:,2+bin)>0 & FR_conn_Chains_go(:,1)>=FRRange(1) & FR_conn_Chains_go(:,1)<FRRange(2) & FR_conn_Chains_go(:,2)>=FRRange(1) & FR_conn_Chains_go(:,2)<FRRange(2));
        IsSigFC_hit_ConAct{1,i}{1,bin} = zeros(tempPairsNum_hit_ConAct,1);
        if tempPairsNum_hit_ConAct > 0
            IsSigFC_hit_ConAct{1,i}{1,bin}(1:tempFcPairsNum_hit_ConAct) = 1;
        end
        % mPFC
        tempPairsNum_hit_mPFC_ConAct = nnz(ismember(Pair_reg(:,1),1) & ismember(Pair_reg(:,2),1) & Pref_pair(:,2+bin)==Pref_pair(:,9+bin) & Pref_pair(:,2+bin)>0 & FRPairChains_go(:,1)>=FRRange(1) & FRPairChains_go(:,1)<FRRange(2) & FRPairChains_go(:,2)>=FRRange(1) & FRPairChains_go(:,2)<FRRange(2));
        tempFcPairsNum_hit_mPFC_ConAct = nnz(ismember(reg_chain_go(:,1),1) & ismember(reg_chain_go(:,2),1) & pref_chain_go(:,2+bin)==pref_chain_go(:,9+bin) & pref_chain_go(:,2+bin)>0 & FR_conn_Chains_go(:,1)>=FRRange(1) & FR_conn_Chains_go(:,1)<FRRange(2) & FR_conn_Chains_go(:,2)>=FRRange(1) & FR_conn_Chains_go(:,2)<FRRange(2));
        IsSigFC_hit_ConAct{2,i}{1,bin} = zeros(tempPairsNum_hit_mPFC_ConAct,1);
        if tempPairsNum_hit_mPFC_ConAct > 0
            IsSigFC_hit_ConAct{2,i}{1,bin}(1:tempFcPairsNum_hit_mPFC_ConAct) = 1;
        end
        % aAIC
        tempParisNum_hit_aAIC_ConAct = nnz(ismember(Pair_reg(:,1),2) & ismember(Pair_reg(:,2),2) & Pref_pair(:,2+bin)==Pref_pair(:,9+bin) & Pref_pair(:,2+bin)>0 & FRPairChains_go(:,1)>=FRRange(1) & FRPairChains_go(:,1)<FRRange(2) & FRPairChains_go(:,2)>=FRRange(1) & FRPairChains_go(:,2)<FRRange(2));
        tempNum_CongActConn_aAIC_Hit = nnz(ismember(reg_chain_go(:,1),2) & ismember(reg_chain_go(:,2),2) & pref_chain_go(:,2+bin)==pref_chain_go(:,9+bin) & pref_chain_go(:,2+bin)>0 & FR_conn_Chains_go(:,1)>=FRRange(1) & FR_conn_Chains_go(:,1)<FRRange(2) & FR_conn_Chains_go(:,2)>=FRRange(1) & FR_conn_Chains_go(:,2)<FRRange(2));
        IsSigFC_hit_ConAct{3,i}{1,bin} = zeros(tempParisNum_hit_aAIC_ConAct,1);
        if tempParisNum_hit_aAIC_ConAct > 0
            IsSigFC_hit_ConAct{3,i}{1,bin}(1:tempNum_CongActConn_aAIC_Hit) = 1;
        end
        % mPFC-aAIC
        tempPairs_hit_mPFCaAIC_ConAct = Pair_chain(ismember(Pair_reg(:,1),2) & ismember(Pair_reg(:,2),1) & Pref_pair(:,2+bin)==Pref_pair(:,9+bin) & Pref_pair(:,2+bin)>0 & FRPairChains_go(:,1)>=FRRange(1) & FRPairChains_go(:,1)<FRRange(2) & FRPairChains_go(:,2)>=FRRange(1) & FRPairChains_go(:,2)<FRRange(2),:);
        tempFcPairs_hit_mPFCaAIC_ConAct = conn_chain_go(ismember(reg_chain_go(:,1),2) & ismember(reg_chain_go(:,2),1) & reg_chain_go(:,3)<0 & pref_chain_go(:,2+bin)==pref_chain_go(:,9+bin) & pref_chain_go(:,2+bin)>0 & FR_conn_Chains_go(:,1)>=FRRange(1) & FR_conn_Chains_go(:,1)<FRRange(2) & FR_conn_Chains_go(:,2)>=FRRange(1) & FR_conn_Chains_go(:,2)<FRRange(2),1:2);
        IsSigFC_hit_ConAct{4,i}{1,bin} = zeros(size(tempPairs_hit_mPFCaAIC_ConAct,1),1);
        for iFcPair = 1:size(tempFcPairs_hit_mPFCaAIC_ConAct,1)
            tempUnit1 = tempFcPairs_hit_mPFCaAIC_ConAct(iFcPair,1);
            tempUnit2 = tempFcPairs_hit_mPFCaAIC_ConAct(iFcPair,2);
            tempPos = find(tempPairs_hit_mPFCaAIC_ConAct(:,1)==tempUnit1 & tempPairs_hit_mPFCaAIC_ConAct(:,2)==tempUnit2);
            if ~isempty(tempPos)
                IsSigFC_hit_ConAct{4,i}{1,bin}(tempPos) = 1;
            end
        end
        % aAIC-mPFC
        tempPairs_hit_aAICmPFC_ConAct = Pair_chain(ismember(Pair_reg(:,1),2) & ismember(Pair_reg(:,2),1) & Pref_pair(:,2+bin)==Pref_pair(:,9+bin) & Pref_pair(:,2+bin)>0 & FRPairChains_go(:,1)>=FRRange(1) & FRPairChains_go(:,1)<FRRange(2) & FRPairChains_go(:,2)>=FRRange(1) & FRPairChains_go(:,2)<FRRange(2),:);
        tempFcPairs_hit_aAICmPFC_ConAct = conn_chain_go(ismember(reg_chain_go(:,1),2) & ismember(reg_chain_go(:,2),1) & reg_chain_go(:,3)>0 & pref_chain_go(:,2+bin)==pref_chain_go(:,9+bin) & pref_chain_go(:,2+bin)>0 & FR_conn_Chains_go(:,1)>=FRRange(1) & FR_conn_Chains_go(:,1)<FRRange(2) & FR_conn_Chains_go(:,2)>=FRRange(1) & FR_conn_Chains_go(:,2)<FRRange(2),1:2);
        IsSigFC_hit_ConAct{5,i}{1,bin} = zeros(size(tempPairs_hit_aAICmPFC_ConAct,1),1);
        for iFcPair = 1:size(tempFcPairs_hit_aAICmPFC_ConAct,1)
            tempUnit1 = tempFcPairs_hit_aAICmPFC_ConAct(iFcPair,1);
            tempUnit2 = tempFcPairs_hit_aAICmPFC_ConAct(iFcPair,2);
            tempPos = find(tempPairs_hit_aAICmPFC_ConAct(:,1)==tempUnit1 & tempPairs_hit_aAICmPFC_ConAct(:,2)==tempUnit2);
            if ~isempty(tempPos)
                IsSigFC_hit_ConAct{5,i}{1,bin}(tempPos) = 1;
            end
        end
        
       %% Congruent inactive pairs
        % mPFC, aAIC, mPFC-aAIC, and aAIC-mPFC
        tempPairsNum_hit_ConInact = nnz(ismember(Pair_reg(:,1),[1 2]) & ismember(Pair_reg(:,2),[1 2]) & max(Pref_pair(:,3:6),[],2)==max(Pref_pair(:,10:13),[],2) & any(Pref_pair(:,[bin+2 bin+9])<1,2) & FRPairChains_go(:,1)>=FRRange(1) & FRPairChains_go(:,1)<FRRange(2) & FRPairChains_go(:,2)>=FRRange(1) & FRPairChains_go(:,2)<FRRange(2));
        tempFcPairsNum_hit_ConInact = nnz(ismember(reg_chain_go(:,1),[1 2]) & ismember(reg_chain_go(:,2),[1 2]) & max(pref_chain_go(:,3:6),[],2)==max(pref_chain_go(:,10:13),[],2) & any(pref_chain_go(:,[bin+2 bin+9])<1,2) & FR_conn_Chains_go(:,1)>=FRRange(1) & FR_conn_Chains_go(:,1)<FRRange(2) & FR_conn_Chains_go(:,2)>=FRRange(1) & FR_conn_Chains_go(:,2)<FRRange(2));
        IsSigFC_hit_ConInact{1,i}{1,bin} = zeros(tempPairsNum_hit_ConInact,1);
        if tempPairsNum_hit_ConInact > 0
            IsSigFC_hit_ConInact{1,i}{1,bin}(1:tempFcPairsNum_hit_ConInact) = 1;
        end
        % mPFC
        tempPairsNum_hit_mPFC_ConInact = nnz(ismember(Pair_reg(:,1),1) & ismember(Pair_reg(:,2),1) & max(Pref_pair(:,3:6),[],2)==max(Pref_pair(:,10:13),[],2) & any(Pref_pair(:,[bin+2 bin+9])<1,2) & FRPairChains_go(:,1)>=FRRange(1) & FRPairChains_go(:,1)<FRRange(2) & FRPairChains_go(:,2)>=FRRange(1) & FRPairChains_go(:,2)<FRRange(2));
        tempFcPairsNum_hit_mPFC_ConInact = nnz(ismember(reg_chain_go(:,1),1) & ismember(reg_chain_go(:,2),1) & max(pref_chain_go(:,3:6),[],2)==max(pref_chain_go(:,10:13),[],2) & any(pref_chain_go(:,[bin+2 bin+9])<1,2) & FR_conn_Chains_go(:,1)>=FRRange(1) & FR_conn_Chains_go(:,1)<FRRange(2) & FR_conn_Chains_go(:,2)>=FRRange(1) & FR_conn_Chains_go(:,2)<FRRange(2));
        IsSigFC_hit_ConInact{2,i}{1,bin} = zeros(tempPairsNum_hit_mPFC_ConInact,1);
        if tempPairsNum_hit_mPFC_ConInact > 0
            IsSigFC_hit_ConInact{2,i}{1,bin}(1:tempFcPairsNum_hit_mPFC_ConInact) = 1;
        end
        % aAIC
        tempPairsNum_hit_aAIC_ConInact = nnz(ismember(Pair_reg(:,1),2) & ismember(Pair_reg(:,2),2) & max(Pref_pair(:,3:6),[],2)==max(Pref_pair(:,10:13),[],2) & any(Pref_pair(:,[bin+2 bin+9])<1,2) & FRPairChains_go(:,1)>=FRRange(1) & FRPairChains_go(:,1)<FRRange(2) & FRPairChains_go(:,2)>=FRRange(1) & FRPairChains_go(:,2)<FRRange(2));
        tempFcPairsNum_hit_aAIC_ConInact = nnz(ismember(reg_chain_go(:,1),2) & ismember(reg_chain_go(:,2),2) & max(pref_chain_go(:,3:6),[],2)==max(pref_chain_go(:,10:13),[],2) & any(pref_chain_go(:,[bin+2 bin+9])<1,2) & FR_conn_Chains_go(:,1)>=FRRange(1) & FR_conn_Chains_go(:,1)<FRRange(2) & FR_conn_Chains_go(:,2)>=FRRange(1) & FR_conn_Chains_go(:,2)<FRRange(2));
        IsSigFC_hit_ConInact{3,i}{1,bin} = zeros(tempPairsNum_hit_aAIC_ConInact,1);
        if tempPairsNum_hit_aAIC_ConInact > 0
            IsSigFC_hit_ConInact{3,i}{1,bin}(1:tempFcPairsNum_hit_aAIC_ConInact) = 1;
        end
        % mPFC-aAIC
        tempPairs_hit_mPFCaAIC_ConInact = Pair_chain(ismember(Pair_reg(:,1),2) & ismember(Pair_reg(:,2),1) & max(Pref_pair(:,3:6),[],2)==max(Pref_pair(:,10:13),[],2) & any(Pref_pair(:,[bin+2 bin+9])<1,2) & FRPairChains_go(:,1)>=FRRange(1) & FRPairChains_go(:,1)<FRRange(2) & FRPairChains_go(:,2)>=FRRange(1) & FRPairChains_go(:,2)<FRRange(2),:);
        tempFcPairs_hit_mPFCaAIC_ConInact = conn_chain_go(ismember(reg_chain_go(:,1),2) & ismember(reg_chain_go(:,2),1) & reg_chain_go(:,3)<0 & max(pref_chain_go(:,3:6),[],2)==max(pref_chain_go(:,10:13),[],2) & any(pref_chain_go(:,[bin+2 bin+9])<1,2) & FR_conn_Chains_go(:,1)>=FRRange(1) & FR_conn_Chains_go(:,1)<FRRange(2) & FR_conn_Chains_go(:,2)>=FRRange(1) & FR_conn_Chains_go(:,2)<FRRange(2),1:2);
        IsSigFC_hit_ConInact{4,i}{1,bin} = zeros(size(tempPairs_hit_mPFCaAIC_ConInact,1),1);
        for iFcPair = 1:size(tempFcPairs_hit_mPFCaAIC_ConInact,1)
            tempUnit1 = tempFcPairs_hit_mPFCaAIC_ConInact(iFcPair,1);
            tempUnit2 = tempFcPairs_hit_mPFCaAIC_ConInact(iFcPair,2);
            tempPos = find(tempPairs_hit_mPFCaAIC_ConInact(:,1)==tempUnit1 & tempPairs_hit_mPFCaAIC_ConInact(:,2)==tempUnit2);
            if ~isempty(tempPos)
                IsSigFC_hit_ConInact{4,i}{1,bin}(tempPos) = 1;
            end
        end
        % aAIC-mPFC
        tempPairs_hit_aAICmPFC_ConInact = Pair_chain(ismember(Pair_reg(:,1),2) & ismember(Pair_reg(:,2),1) & max(Pref_pair(:,3:6),[],2)==max(Pref_pair(:,10:13),[],2) & any(Pref_pair(:,[bin+2 bin+9])<1,2) & FRPairChains_go(:,1)>=FRRange(1) & FRPairChains_go(:,1)<FRRange(2) & FRPairChains_go(:,2)>=FRRange(1) & FRPairChains_go(:,2)<FRRange(2),:);
        tempFcPairs_hit_aAICmPFC_ConInact = conn_chain_go(ismember(reg_chain_go(:,1),2) & ismember(reg_chain_go(:,2),1) & reg_chain_go(:,3)>0 & max(pref_chain_go(:,3:6),[],2)==max(pref_chain_go(:,10:13),[],2) & any(pref_chain_go(:,[bin+2 bin+9])<1,2) & FR_conn_Chains_go(:,1)>=FRRange(1) & FR_conn_Chains_go(:,1)<FRRange(2) & FR_conn_Chains_go(:,2)>=FRRange(1) & FR_conn_Chains_go(:,2)<FRRange(2),1:2);
        IsSigFC_hit_ConInact{5,i}{1,bin} = zeros(size(tempPairs_hit_aAICmPFC_ConInact,1),1);
        for iFcPair = 1:size(tempFcPairs_hit_aAICmPFC_ConInact,1)
            tempUnit1 = tempFcPairs_hit_aAICmPFC_ConInact(iFcPair,1);
            tempUnit2 = tempFcPairs_hit_aAICmPFC_ConInact(iFcPair,2);
            tempPos = find(tempPairs_hit_aAICmPFC_ConInact(:,1)==tempUnit1 & tempPairs_hit_aAICmPFC_ConInact(:,2)==tempUnit2);
            if ~isempty(tempPos)
                IsSigFC_hit_ConInact{5,i}{1,bin}(tempPos) = 1;
            end
        end
       %% Incongruent active pairs
        % mPFC, aAIC, mPFC-aAIC, and aAIC-mPFC
        tempPairsNum_hit_InconAct = nnz(ismember(Pair_reg(:,1),[1 2]) & ismember(Pair_reg(:,2),[1 2]) & Pref_pair(:,2+bin)>0 & Pref_pair(:,9+bin)>0 & Pref_pair(:,2+bin)~=Pref_pair(:,9+bin) & FRPairChains_go(:,1)>=FRRange(1) & FRPairChains_go(:,1)<FRRange(2) & FRPairChains_go(:,2)>=FRRange(1) & FRPairChains_go(:,2)<FRRange(2));
        tempNum_IncongActConn_Hit = nnz(ismember(reg_chain_go(:,1),[1 2]) & ismember(reg_chain_go(:,2),[1 2]) & pref_chain_go(:,2+bin)>0 & pref_chain_go(:,9+bin)>0 & pref_chain_go(:,2+bin)~=pref_chain_go(:,9+bin) & FR_conn_Chains_go(:,1)>=FRRange(1) & FR_conn_Chains_go(:,1)<FRRange(2) & FR_conn_Chains_go(:,2)>=FRRange(1) & FR_conn_Chains_go(:,2)<FRRange(2));
        IsSigFC_hit_InconAct{1,i}{1,bin} = zeros(tempPairsNum_hit_InconAct,1);
        if tempPairsNum_hit_InconAct > 0
            IsSigFC_hit_InconAct{1,i}{1,bin}(1:tempNum_IncongActConn_Hit) = 1;
        end
        % mPFC
        tempPairsNum_hit_mPFC_InconAct = nnz(ismember(Pair_reg(:,1),1) & ismember(Pair_reg(:,2),1) & Pref_pair(:,2+bin)>0 & Pref_pair(:,9+bin)>0 & Pref_pair(:,2+bin)~=Pref_pair(:,9+bin) & FRPairChains_go(:,1)>=FRRange(1) & FRPairChains_go(:,1)<FRRange(2) & FRPairChains_go(:,2)>=FRRange(1) & FRPairChains_go(:,2)<FRRange(2));
        tempFcPairsNum_hit_mPFC_InconAct = nnz(ismember(reg_chain_go(:,1),1) & ismember(reg_chain_go(:,2),1) & pref_chain_go(:,2+bin)>0 & pref_chain_go(:,9+bin)>0 & pref_chain_go(:,2+bin)~=pref_chain_go(:,9+bin) & FR_conn_Chains_go(:,1)>=FRRange(1) & FR_conn_Chains_go(:,1)<FRRange(2) & FR_conn_Chains_go(:,2)>=FRRange(1) & FR_conn_Chains_go(:,2)<FRRange(2));
        IsSigFC_hit_InconAct{2,i}{1,bin} = zeros(tempPairsNum_hit_mPFC_InconAct,1);
        if tempPairsNum_hit_mPFC_InconAct > 0
            IsSigFC_hit_InconAct{2,i}{1,bin}(1:tempFcPairsNum_hit_mPFC_InconAct) = 1;
        end
        % aAIC
        tempPairsNum_hit_aAIC_InconAct = nnz(ismember(Pair_reg(:,1),2) & ismember(Pair_reg(:,2),2) & Pref_pair(:,2+bin)>0 & Pref_pair(:,9+bin)>0 & Pref_pair(:,2+bin)~=Pref_pair(:,9+bin) & FRPairChains_go(:,1)>=FRRange(1) & FRPairChains_go(:,1)<FRRange(2) & FRPairChains_go(:,2)>=FRRange(1) & FRPairChains_go(:,2)<FRRange(2));
        tempFcPairsNum_hit_aAIC_InconAct = nnz(ismember(reg_chain_go(:,1),2) & ismember(reg_chain_go(:,2),2) & pref_chain_go(:,2+bin)>0 & pref_chain_go(:,9+bin)>0 & pref_chain_go(:,2+bin)~=pref_chain_go(:,9+bin) & FR_conn_Chains_go(:,1)>=FRRange(1) & FR_conn_Chains_go(:,1)<FRRange(2) & FR_conn_Chains_go(:,2)>=FRRange(1) & FR_conn_Chains_go(:,2)<FRRange(2));
        IsSigFC_hit_InconAct{3,i}{1,bin} = zeros(tempPairsNum_hit_aAIC_InconAct,1);
        if tempPairsNum_hit_aAIC_InconAct > 0
            IsSigFC_hit_InconAct{3,i}{1,bin}(1:tempFcPairsNum_hit_aAIC_InconAct) = 1;
        end
        % mPFC-aAIC
        tempPairs_hit_mPFCaAIC_InconAct = Pair_chain(ismember(Pair_reg(:,1),2) & ismember(Pair_reg(:,2),1) & Pref_pair(:,2+bin)>0 & Pref_pair(:,9+bin)>0 & Pref_pair(:,2+bin)~=Pref_pair(:,9+bin) & FRPairChains_go(:,1)>=FRRange(1) & FRPairChains_go(:,1)<FRRange(2) & FRPairChains_go(:,2)>=FRRange(1) & FRPairChains_go(:,2)<FRRange(2),:);
        tempFcPairs_hit_mPFCaAIC_InconAct = conn_chain_go(ismember(reg_chain_go(:,1),2) & ismember(reg_chain_go(:,2),1) & reg_chain_go(:,3)<0 & pref_chain_go(:,2+bin)>0 & pref_chain_go(:,9+bin)>0 & pref_chain_go(:,2+bin)~=pref_chain_go(:,9+bin) & FR_conn_Chains_go(:,1)>=FRRange(1) & FR_conn_Chains_go(:,1)<FRRange(2) & FR_conn_Chains_go(:,2)>=FRRange(1) & FR_conn_Chains_go(:,2)<FRRange(2),1:2);
        IsSigFC_hit_InconAct{4,i}{1,bin} = zeros(size(tempPairs_hit_mPFCaAIC_InconAct,1),1);
        for iFcPair = 1:size(tempFcPairs_hit_mPFCaAIC_InconAct,1)
            tempUnit1 = tempFcPairs_hit_mPFCaAIC_InconAct(iFcPair,1);
            tempUnit2 = tempFcPairs_hit_mPFCaAIC_InconAct(iFcPair,2);
            tempPos = find(tempPairs_hit_mPFCaAIC_InconAct(:,1)==tempUnit1 & tempPairs_hit_mPFCaAIC_InconAct(:,2)==tempUnit2);
            if ~isempty(tempPos)
                IsSigFC_hit_InconAct{4,i}{1,bin}(tempPos) = 1;
            end
        end
        % aAIC-mPFC
        tempPairs_hit_aAICmPFC_InconAct = Pair_chain(ismember(Pair_reg(:,1),2) & ismember(Pair_reg(:,2),1) & Pref_pair(:,2+bin)>0 & Pref_pair(:,9+bin)>0 & Pref_pair(:,2+bin)~=Pref_pair(:,9+bin) & FRPairChains_go(:,1)>=FRRange(1) & FRPairChains_go(:,1)<FRRange(2) & FRPairChains_go(:,2)>=FRRange(1) & FRPairChains_go(:,2)<FRRange(2),:);
        tempFcPairs_hit_aAICmPFC_InconAct = conn_chain_go(ismember(reg_chain_go(:,1),2) & ismember(reg_chain_go(:,2),1) & reg_chain_go(:,3)>0 & pref_chain_go(:,2+bin)>0 & pref_chain_go(:,9+bin)>0 & pref_chain_go(:,2+bin)~=pref_chain_go(:,9+bin) & FR_conn_Chains_go(:,1)>=FRRange(1) & FR_conn_Chains_go(:,1)<FRRange(2) & FR_conn_Chains_go(:,2)>=FRRange(1) & FR_conn_Chains_go(:,2)<FRRange(2),1:2);
        IsSigFC_hit_InconAct{5,i}{1,bin} = zeros(size(tempPairs_hit_aAICmPFC_InconAct,1),1);
        for iFcPair = 1:size(tempFcPairs_hit_aAICmPFC_InconAct,1)
            tempUnit1 = tempFcPairs_hit_aAICmPFC_InconAct(iFcPair,1);
            tempUnit2 = tempFcPairs_hit_aAICmPFC_InconAct(iFcPair,2);
            tempPos = find(tempPairs_hit_aAICmPFC_InconAct(:,1)==tempUnit1 & tempPairs_hit_aAICmPFC_InconAct(:,2)==tempUnit2);
            if ~isempty(tempPos)
                IsSigFC_hit_InconAct{5,i}{1,bin}(tempPos) = 1;
            end
        end
        
       %% Incongruent inactive pairs
        % mPFC, aAIC, mPFC-aAIC, and aAIC-mPFC
        tempPairsNum_hit_InconInact = nnz(ismember(Pair_reg(:,1),[1 2]) & ismember(Pair_reg(:,2),[1 2]) & max(Pref_pair(:,3:6),[],2)~=max(Pref_pair(:,10:13),[],2) & max(Pref_pair(:,3:6),[],2)>0 & max(Pref_pair(:,10:13),[],2)>0 & any(Pref_pair(:,[bin+2 bin+9])<1,2) & FRPairChains_go(:,1)>=FRRange(1) & FRPairChains_go(:,1)<FRRange(2) & FRPairChains_go(:,2)>=FRRange(1) & FRPairChains_go(:,2)<FRRange(2));
        tempFcPairsNum_hit_InconInact = nnz(ismember(reg_chain_go(:,1),[1 2]) & ismember(reg_chain_go(:,2),[1 2]) & max(pref_chain_go(:,3:6),[],2)~=max(pref_chain_go(:,10:13),[],2) & max(pref_chain_go(:,3:6),[],2)>0 & max(pref_chain_go(:,10:13),[],2)>0 & any(pref_chain_go(:,[bin+2 bin+9])<1,2) & FR_conn_Chains_go(:,1)>=FRRange(1) & FR_conn_Chains_go(:,1)<FRRange(2) & FR_conn_Chains_go(:,2)>=FRRange(1) & FR_conn_Chains_go(:,2)<FRRange(2));
        IsSigFC_hit_InconInact{1,i}{1,bin} = zeros(tempPairsNum_hit_InconInact,1);
        if tempPairsNum_hit_InconInact > 0
            IsSigFC_hit_InconInact{1,i}{1,bin}(1:tempFcPairsNum_hit_InconInact) = 1;
        end
        % mPFC
        tempPairsNum_hit_mPFC_InconInact = nnz(ismember(Pair_reg(:,1),1) & ismember(Pair_reg(:,2),1) & max(Pref_pair(:,3:6),[],2)~=max(Pref_pair(:,10:13),[],2) & max(Pref_pair(:,3:6),[],2)>0 & max(Pref_pair(:,10:13),[],2)>0 & any(Pref_pair(:,[bin+2 bin+9])<1,2) & FRPairChains_go(:,1)>=FRRange(1) & FRPairChains_go(:,1)<FRRange(2) & FRPairChains_go(:,2)>=FRRange(1) & FRPairChains_go(:,2)<FRRange(2));
        tempFcPairsNum_hit_mPFC_InconInact = nnz(ismember(reg_chain_go(:,1),1) & ismember(reg_chain_go(:,2),1) & max(pref_chain_go(:,3:6),[],2)~=max(pref_chain_go(:,10:13),[],2) & max(pref_chain_go(:,3:6),[],2)>0 & max(pref_chain_go(:,10:13),[],2)>0 & any(pref_chain_go(:,[bin+2 bin+9])<1,2) & FR_conn_Chains_go(:,1)>=FRRange(1) & FR_conn_Chains_go(:,1)<FRRange(2) & FR_conn_Chains_go(:,2)>=FRRange(1) & FR_conn_Chains_go(:,2)<FRRange(2));
        IsSigFC_hit_InconInact{2,i}{1,bin} = zeros(tempPairsNum_hit_mPFC_InconInact,1);
        if tempPairsNum_hit_mPFC_InconInact > 0
            IsSigFC_hit_InconInact{2,i}{1,bin}(1:tempFcPairsNum_hit_mPFC_InconInact) = 1;
        end
        % aAIC
        tempPairsNum_hit_aAIC_InconInact = nnz(ismember(Pair_reg(:,1),2) & ismember(Pair_reg(:,2),2) & max(Pref_pair(:,3:6),[],2)~=max(Pref_pair(:,10:13),[],2) & max(Pref_pair(:,3:6),[],2)>0 & max(Pref_pair(:,10:13),[],2)>0 & any(Pref_pair(:,[bin+2 bin+9])<1,2) & FRPairChains_go(:,1)>=FRRange(1) & FRPairChains_go(:,1)<FRRange(2) & FRPairChains_go(:,2)>=FRRange(1) & FRPairChains_go(:,2)<FRRange(2));
        tempFcPairsNum_hit_aAIC_InconInact = nnz(ismember(reg_chain_go(:,1),2) & ismember(reg_chain_go(:,2),2) & max(pref_chain_go(:,3:6),[],2)~=max(pref_chain_go(:,10:13),[],2) & max(pref_chain_go(:,3:6),[],2)>0 & max(pref_chain_go(:,10:13),[],2)>0 & any(pref_chain_go(:,[bin+2 bin+9])<1,2) & FR_conn_Chains_go(:,1)>=FRRange(1) & FR_conn_Chains_go(:,1)<FRRange(2) & FR_conn_Chains_go(:,2)>=FRRange(1) & FR_conn_Chains_go(:,2)<FRRange(2));
        IsSigFC_hit_InconInact{3,i}{1,bin} = zeros(tempPairsNum_hit_aAIC_InconInact,1);
        if tempPairsNum_hit_aAIC_InconInact > 0
            IsSigFC_hit_InconInact{3,i}{1,bin}(1:tempFcPairsNum_hit_aAIC_InconInact) = 1;
        end
        % mPFC-aAIC
        tempPairs_hit_mPFCaAIC_InconInact = Pair_chain(ismember(Pair_reg(:,1),2) & ismember(Pair_reg(:,2),1) & max(Pref_pair(:,3:6),[],2)~=max(Pref_pair(:,10:13),[],2) & max(Pref_pair(:,3:6),[],2)>0 & max(Pref_pair(:,10:13),[],2)>0 & any(Pref_pair(:,[bin+2 bin+9])<1,2) & FRPairChains_go(:,1)>=FRRange(1) & FRPairChains_go(:,1)<FRRange(2) & FRPairChains_go(:,2)>=FRRange(1) & FRPairChains_go(:,2)<FRRange(2),:);
        tempFcPairs_hit_mPFCaAIC_InconInact = conn_chain_go(ismember(reg_chain_go(:,1),2) & ismember(reg_chain_go(:,2),1) & reg_chain_go(:,3)<0 & max(pref_chain_go(:,3:6),[],2)~=max(pref_chain_go(:,10:13),[],2) & max(pref_chain_go(:,3:6),[],2)>0 & max(pref_chain_go(:,10:13),[],2)>0 & any(pref_chain_go(:,[bin+2 bin+9])<1,2) & FR_conn_Chains_go(:,1)>=FRRange(1) & FR_conn_Chains_go(:,1)<FRRange(2) & FR_conn_Chains_go(:,2)>=FRRange(1) & FR_conn_Chains_go(:,2)<FRRange(2),1:2);
        IsSigFC_hit_InconInact{4,i}{1,bin} = zeros(size(tempPairs_hit_mPFCaAIC_InconInact,1),1);
        for iFcPair = 1:size(tempFcPairs_hit_mPFCaAIC_InconInact,1)
            tempUnit1 = tempFcPairs_hit_mPFCaAIC_InconInact(iFcPair,1);
            tempUnit2 = tempFcPairs_hit_mPFCaAIC_InconInact(iFcPair,2);
            tempPos = find(tempPairs_hit_mPFCaAIC_InconInact(:,1)==tempUnit1 & tempPairs_hit_mPFCaAIC_InconInact(:,2)==tempUnit2);
            if ~isempty(tempPos)
                IsSigFC_hit_InconInact{4,i}{1,bin}(tempPos) = 1;
            end
        end
        % aAIC-mPFC
        tempPairs_hit_aAICmPFC_InconInact = Pair_chain(ismember(Pair_reg(:,1),2) & ismember(Pair_reg(:,2),1) & max(Pref_pair(:,3:6),[],2)~=max(Pref_pair(:,10:13),[],2) & max(Pref_pair(:,3:6),[],2)>0 & max(Pref_pair(:,10:13),[],2)>0 & any(Pref_pair(:,[bin+2 bin+9])<1,2) & FRPairChains_go(:,1)>=FRRange(1) & FRPairChains_go(:,1)<FRRange(2) & FRPairChains_go(:,2)>=FRRange(1) & FRPairChains_go(:,2)<FRRange(2),:);
        tempFcPairs_hit_aAICmPFC_InconInact = conn_chain_go(ismember(reg_chain_go(:,1),2) & ismember(reg_chain_go(:,2),1) & reg_chain_go(:,3)>0 & max(pref_chain_go(:,3:6),[],2)~=max(pref_chain_go(:,10:13),[],2) & max(pref_chain_go(:,3:6),[],2)>0 & max(pref_chain_go(:,10:13),[],2)>0 & any(pref_chain_go(:,[bin+2 bin+9])<1,2) & FR_conn_Chains_go(:,1)>=FRRange(1) & FR_conn_Chains_go(:,1)<FRRange(2) & FR_conn_Chains_go(:,2)>=FRRange(1) & FR_conn_Chains_go(:,2)<FRRange(2),1:2);
        IsSigFC_hit_InconInact{5,i}{1,bin} = zeros(size(tempPairs_hit_aAICmPFC_InconInact,1),1);
        for iFcPair = 1:size(tempFcPairs_hit_aAICmPFC_InconInact,1)
            tempUnit1 = tempFcPairs_hit_aAICmPFC_InconInact(iFcPair,1);
            tempUnit2 = tempFcPairs_hit_aAICmPFC_InconInact(iFcPair,2);
            tempPos = find(tempPairs_hit_aAICmPFC_InconInact(:,1)==tempUnit1 & tempPairs_hit_aAICmPFC_InconInact(:,2)==tempUnit2);
            if ~isempty(tempPos)
                IsSigFC_hit_InconInact{5,i}{1,bin}(tempPos) = 1;
            end
        end
    end
end

%% Save variables
save(['Region-specific FC density within FR ' num2str(TarFR(1)) ' to ' num2str(TarFR(1)+5) ' to ' num2str(TarFR(end)) '-' Group '_1msbin'],'IsSigFC_hit_ConAct','IsSigFC_hit_ConInact','IsSigFC_hit_InconAct','IsSigFC_hit_InconInact','IsSigFC_hit_Nonmem','-v7.3');

