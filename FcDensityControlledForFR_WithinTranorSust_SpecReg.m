%% Region-specific FC density of neuronal pairs controlled for FR
% //////within transient, within sustained neurons//////
clear; clc; close all;

%% Assignment
Group = 'ChR2Group';
DelayBinNum = 4;
TarFR = [2 13];

%% Region-specific FC of various types of neuronal pairs, controlled for FR
% //////each row represents each region, 1: all regions; 2: mPFC; 3: aAIC; 4: mPFC-aAIC; 5: aAIC-mPFC//////
IsSigFC_hit_WithinTran_ConAct = cell(4,numel(TarFR)); % congruent active
IsSigFC_hit_WithinSust_ConAct = cell(4,numel(TarFR)); 
IsSigFC_hit_WithinTran_ConInact = cell(4,numel(TarFR)); % congruent inactive
IsSigFC_hit_WithinTran_InconAct = cell(4,numel(TarFR)); % incongruent active
IsSigFC_hit_WithinSust_InconAct = cell(4,numel(TarFR));
IsSigFC_hit_WithinTran_InconInact = cell(4,numel(TarFR)); % incongruent inactive
IsSigFC_hit_Nonmem = cell(4,numel(TarFR));
for i = 1:numel(TarFR)
    if i==1
        FRRange = [TarFR(i) TarFR(i)+5];
    elseif i==2
        FRRange = [TarFR(i)-6 TarFR(i)];
    end
    
    %% FC density of pairs of memory neurons
    for bin = 1:DelayBinNum % four 1-second bins of delay
        load(sprintf('FCofMemoryNeurons_PropertyPopulation_%d_%d_%s.mat',bin,bin+1,Group));
        load(sprintf('PairsFRModulation_selec_conn_chain_delay4s_%d_%d_%s.mat',bin,bin+1,Group));
        
       %% Congruent active pairs
        % mPFC, aAIC, mPFC-aAIC, and aAIC-mPFC
        tempPairsNum_hit_WithinTran_ConAct = nnz(ismember(Pair_reg(:,1),[1 2]) & ismember(Pair_reg(:,2),[1 2]) & Pref_pair(:,2+bin)==Pref_pair(:,9+bin) & Pref_pair(:,2+bin)>0 & ~all(Pref_pair(:,3:6),2) & ~all(Pref_pair(:,10:13),2) & FRPairChains_go(:,1)>=FRRange(1) & FRPairChains_go(:,1)<FRRange(2) & FRPairChains_go(:,2)>=FRRange(1) & FRPairChains_go(:,2)<FRRange(2));
        tempFcPairsNum_hit_WithinTran_ConAct = nnz(ismember(reg_chain_go(:,1),[1 2]) & ismember(reg_chain_go(:,2),[1 2]) & pref_chain_go(:,2+bin)==pref_chain_go(:,9+bin) & pref_chain_go(:,2+bin)>0 & ~all(pref_chain_go(:,3:6),2) & ~all(pref_chain_go(:,10:13),2) & FR_conn_Chains_go(:,1)>=FRRange(1) & FR_conn_Chains_go(:,1)<FRRange(2) & FR_conn_Chains_go(:,2)>=FRRange(1) & FR_conn_Chains_go(:,2)<FRRange(2));
        IsSigFC_hit_WithinTran_ConAct{1,i}{1,bin} = zeros(tempPairsNum_hit_WithinTran_ConAct,1);
        if tempPairsNum_hit_WithinTran_ConAct > 0
            IsSigFC_hit_WithinTran_ConAct{1,i}{1,bin}(1:tempFcPairsNum_hit_WithinTran_ConAct) = 1; % within transient neurons
        end
        tempPairsNum_hit_WithinSust_ConAct = nnz(ismember(Pair_reg(:,1),[1 2]) & ismember(Pair_reg(:,2),[1 2]) & Pref_pair(:,2+bin)==Pref_pair(:,9+bin) & Pref_pair(:,2+bin)>0 & all(Pref_pair(:,3:6),2) & all(Pref_pair(:,10:13),2) & FRPairChains_go(:,1)>=FRRange(1) & FRPairChains_go(:,1)<FRRange(2) & FRPairChains_go(:,2)>=FRRange(1) & FRPairChains_go(:,2)<FRRange(2));
        tempFcPairsNum_hit_WithinSust_ConAct = nnz(ismember(reg_chain_go(:,1),[1 2]) & ismember(reg_chain_go(:,2),[1 2]) & pref_chain_go(:,2+bin)==pref_chain_go(:,9+bin) & pref_chain_go(:,2+bin)>0 & all(pref_chain_go(:,3:6),2) & all(pref_chain_go(:,10:13),2) & FR_conn_Chains_go(:,1)>=FRRange(1) & FR_conn_Chains_go(:,1)<FRRange(2) & FR_conn_Chains_go(:,2)>=FRRange(1) & FR_conn_Chains_go(:,2)<FRRange(2));
        IsSigFC_hit_WithinSust_ConAct{1,i}{1,bin} = zeros(tempPairsNum_hit_WithinSust_ConAct,1);
        if tempPairsNum_hit_WithinSust_ConAct > 0
            IsSigFC_hit_WithinSust_ConAct{1,i}{1,bin}(1:tempFcPairsNum_hit_WithinSust_ConAct) = 1; % within sustained neurons
        end
        % mPFC
        tempPairsNum_hit_WithinTran_mPFC_ConAct = nnz(ismember(Pair_reg(:,1),1) & ismember(Pair_reg(:,2),1) & Pref_pair(:,2+bin)==Pref_pair(:,9+bin) & Pref_pair(:,2+bin)>0 & ~all(Pref_pair(:,3:6),2) & ~all(Pref_pair(:,10:13),2) & FRPairChains_go(:,1)>=FRRange(1) & FRPairChains_go(:,1)<FRRange(2) & FRPairChains_go(:,2)>=FRRange(1) & FRPairChains_go(:,2)<FRRange(2));
        tempFcPairsNum_hit_WithinTran_mPFC_ConAct = nnz(ismember(reg_chain_go(:,1),1) & ismember(reg_chain_go(:,2),1) & pref_chain_go(:,2+bin)==pref_chain_go(:,9+bin) & pref_chain_go(:,2+bin)>0 & ~all(pref_chain_go(:,3:6),2) & ~all(pref_chain_go(:,10:13),2) & FR_conn_Chains_go(:,1)>=FRRange(1) & FR_conn_Chains_go(:,1)<FRRange(2) & FR_conn_Chains_go(:,2)>=FRRange(1) & FR_conn_Chains_go(:,2)<FRRange(2));
        IsSigFC_hit_WithinTran_ConAct{2,i}{1,bin} = zeros(tempPairsNum_hit_WithinTran_mPFC_ConAct,1);
        if tempPairsNum_hit_WithinTran_mPFC_ConAct > 0
            IsSigFC_hit_WithinTran_ConAct{2,i}{1,bin}(1:tempFcPairsNum_hit_WithinTran_mPFC_ConAct) = 1; % within transient neurons
        end
        tempPairsNum_hit_WithinSust_mPFC_ConAct = nnz(ismember(Pair_reg(:,1),1) & ismember(Pair_reg(:,2),1) & Pref_pair(:,2+bin)==Pref_pair(:,9+bin) & Pref_pair(:,2+bin)>0 & all(Pref_pair(:,3:6),2) & all(Pref_pair(:,10:13),2) & FRPairChains_go(:,1)>=FRRange(1) & FRPairChains_go(:,1)<FRRange(2) & FRPairChains_go(:,2)>=FRRange(1) & FRPairChains_go(:,2)<FRRange(2));
        tempFcPairsNum_hit_WithinSust_mPFC_ConAct = nnz(ismember(reg_chain_go(:,1),1) & ismember(reg_chain_go(:,2),1) & pref_chain_go(:,2+bin)==pref_chain_go(:,9+bin) & pref_chain_go(:,2+bin)>0 & all(pref_chain_go(:,3:6),2) & all(pref_chain_go(:,10:13),2) & FR_conn_Chains_go(:,1)>=FRRange(1) & FR_conn_Chains_go(:,1)<FRRange(2) & FR_conn_Chains_go(:,2)>=FRRange(1) & FR_conn_Chains_go(:,2)<FRRange(2));
        IsSigFC_hit_WithinSust_ConAct{2,i}{1,bin} = zeros(tempPairsNum_hit_WithinSust_mPFC_ConAct,1);
        if tempPairsNum_hit_WithinSust_mPFC_ConAct > 0
            IsSigFC_hit_WithinSust_ConAct{2,i}{1,bin}(1:tempFcPairsNum_hit_WithinSust_mPFC_ConAct) = 1; % within sustained neurons
        end
        % aAIC
        tempParisNum_hit_WithinTran_aAIC_ConAct = nnz(ismember(Pair_reg(:,1),2) & ismember(Pair_reg(:,2),2) & Pref_pair(:,2+bin)==Pref_pair(:,9+bin) & Pref_pair(:,2+bin)>0 & ~all(Pref_pair(:,3:6),2) & ~all(Pref_pair(:,10:13),2) & FRPairChains_go(:,1)>=FRRange(1) & FRPairChains_go(:,1)<FRRange(2) & FRPairChains_go(:,2)>=FRRange(1) & FRPairChains_go(:,2)<FRRange(2));
        tempFcPairsNum_hit_WithinTran_aAIC_ConAct = nnz(ismember(reg_chain_go(:,1),2) & ismember(reg_chain_go(:,2),2) & pref_chain_go(:,2+bin)==pref_chain_go(:,9+bin) & pref_chain_go(:,2+bin)>0 & ~all(pref_chain_go(:,3:6),2) & ~all(pref_chain_go(:,10:13),2) & FR_conn_Chains_go(:,1)>=FRRange(1) & FR_conn_Chains_go(:,1)<FRRange(2) & FR_conn_Chains_go(:,2)>=FRRange(1) & FR_conn_Chains_go(:,2)<FRRange(2));
        IsSigFC_hit_WithinTran_ConAct{3,i}{1,bin} = zeros(tempParisNum_hit_WithinTran_aAIC_ConAct,1);
        if tempParisNum_hit_WithinTran_aAIC_ConAct > 0
            IsSigFC_hit_WithinTran_ConAct{3,i}{1,bin}(1:tempFcPairsNum_hit_WithinTran_aAIC_ConAct) = 1; % within transient neurons
        end
        tempParisNum_hit_WithinSust_aAIC_ConAct = nnz(ismember(Pair_reg(:,1),2) & ismember(Pair_reg(:,2),2) & Pref_pair(:,2+bin)==Pref_pair(:,9+bin) & Pref_pair(:,2+bin)>0 & all(Pref_pair(:,3:6),2) & all(Pref_pair(:,10:13),2) & FRPairChains_go(:,1)>=FRRange(1) & FRPairChains_go(:,1)<FRRange(2) & FRPairChains_go(:,2)>=FRRange(1) & FRPairChains_go(:,2)<FRRange(2));
        tempFcPairsNum_hit_WithinSust_aAIC_ConAct = nnz(ismember(reg_chain_go(:,1),2) & ismember(reg_chain_go(:,2),2) & pref_chain_go(:,2+bin)==pref_chain_go(:,9+bin) & pref_chain_go(:,2+bin)>0 & all(pref_chain_go(:,3:6),2) & all(pref_chain_go(:,10:13),2) & FR_conn_Chains_go(:,1)>=FRRange(1) & FR_conn_Chains_go(:,1)<FRRange(2) & FR_conn_Chains_go(:,2)>=FRRange(1) & FR_conn_Chains_go(:,2)<FRRange(2));
        IsSigFC_hit_WithinSust_ConAct{3,i}{1,bin} = zeros(tempParisNum_hit_WithinSust_aAIC_ConAct,1);
        if tempParisNum_hit_WithinSust_aAIC_ConAct > 0
            IsSigFC_hit_WithinSust_ConAct{3,i}{1,bin}(1:tempFcPairsNum_hit_WithinSust_aAIC_ConAct) = 1; % within sustained neurons
        end
        % mPFC-aAIC
        tempPairsNum_hit_WithinTran_mPFCaAIC_ConAct = nnz(ismember(Pair_reg(:,1),2) & ismember(Pair_reg(:,2),1) & Pref_pair(:,2+bin)==Pref_pair(:,9+bin) & Pref_pair(:,2+bin)>0 & ~all(Pref_pair(:,3:6),2) & ~all(Pref_pair(:,10:13),2) & FRPairChains_go(:,1)>=FRRange(1) & FRPairChains_go(:,1)<FRRange(2) & FRPairChains_go(:,2)>=FRRange(1) & FRPairChains_go(:,2)<FRRange(2));
        tempFcPairsNum_hit_WithinTran_mPFCaAIC_ConAct = nnz(ismember(reg_chain_go(:,1),2) & ismember(reg_chain_go(:,2),1) & reg_chain_go(:,3)<0 & pref_chain_go(:,2+bin)==pref_chain_go(:,9+bin) & pref_chain_go(:,2+bin)>0 & ~all(pref_chain_go(:,3:6),2) & ~all(pref_chain_go(:,10:13),2) & FR_conn_Chains_go(:,1)>=FRRange(1) & FR_conn_Chains_go(:,1)<FRRange(2) & FR_conn_Chains_go(:,2)>=FRRange(1) & FR_conn_Chains_go(:,2)<FRRange(2));
        IsSigFC_hit_WithinTran_ConAct{4,i}{1,bin} = zeros(tempPairsNum_hit_WithinTran_mPFCaAIC_ConAct,1);
        if tempPairsNum_hit_WithinTran_mPFCaAIC_ConAct > 0
            IsSigFC_hit_WithinTran_ConAct{4,i}{1,bin}(1:tempFcPairsNum_hit_WithinTran_mPFCaAIC_ConAct) = 1; % within transient neurons
        end
        tempPairsNum_hit_WithinSust_mPFCaAIC_ConAct = nnz(ismember(Pair_reg(:,1),2) & ismember(Pair_reg(:,2),1) & Pref_pair(:,2+bin)==Pref_pair(:,9+bin) & Pref_pair(:,2+bin)>0 & all(Pref_pair(:,3:6),2) & all(Pref_pair(:,10:13),2) & FRPairChains_go(:,1)>=FRRange(1) & FRPairChains_go(:,1)<FRRange(2) & FRPairChains_go(:,2)>=FRRange(1) & FRPairChains_go(:,2)<FRRange(2));
        tempFcPairsNum_hit_WithinSust_mPFCaAIC_ConAct = nnz(ismember(reg_chain_go(:,1),2) & ismember(reg_chain_go(:,2),1) & reg_chain_go(:,3)<0 & pref_chain_go(:,2+bin)==pref_chain_go(:,9+bin) & pref_chain_go(:,2+bin)>0 & all(pref_chain_go(:,3:6),2) & all(pref_chain_go(:,10:13),2) & FR_conn_Chains_go(:,1)>=FRRange(1) & FR_conn_Chains_go(:,1)<FRRange(2) & FR_conn_Chains_go(:,2)>=FRRange(1) & FR_conn_Chains_go(:,2)<FRRange(2));
        IsSigFC_hit_WithinSust_ConAct{4,i}{1,bin} = zeros(tempPairsNum_hit_WithinSust_mPFCaAIC_ConAct,1);
        if tempPairsNum_hit_WithinSust_mPFCaAIC_ConAct > 0
            IsSigFC_hit_WithinSust_ConAct{4,i}{1,bin}(1:tempFcPairsNum_hit_WithinSust_mPFCaAIC_ConAct) = 1; % within sustained neurons
        end
        % aAIC-mPFC
        tempPairsNum_hit_WithinTran_aAICmPFC_ConAct = nnz(ismember(Pair_reg(:,1),2) & ismember(Pair_reg(:,2),1) & Pref_pair(:,2+bin)==Pref_pair(:,9+bin) & Pref_pair(:,2+bin)>0 & ~all(Pref_pair(:,3:6),2) & ~all(Pref_pair(:,10:13),2) & FRPairChains_go(:,1)>=FRRange(1) & FRPairChains_go(:,1)<FRRange(2) & FRPairChains_go(:,2)>=FRRange(1) & FRPairChains_go(:,2)<FRRange(2));
        tempFcPairsNum_hit_WithinTran_aAICmPFC_ConAct = nnz(ismember(reg_chain_go(:,1),2) & ismember(reg_chain_go(:,2),1) & reg_chain_go(:,3)>0 & pref_chain_go(:,2+bin)==pref_chain_go(:,9+bin) & pref_chain_go(:,2+bin)>0 & ~all(pref_chain_go(:,3:6),2) & ~all(pref_chain_go(:,10:13),2) & FR_conn_Chains_go(:,1)>=FRRange(1) & FR_conn_Chains_go(:,1)<FRRange(2) & FR_conn_Chains_go(:,2)>=FRRange(1) & FR_conn_Chains_go(:,2)<FRRange(2));
        IsSigFC_hit_WithinTran_ConAct{5,i}{1,bin} = zeros(tempPairsNum_hit_WithinTran_aAICmPFC_ConAct,1);
        if tempPairsNum_hit_WithinTran_aAICmPFC_ConAct > 0
            IsSigFC_hit_WithinTran_ConAct{5,i}{1,bin}(1:tempFcPairsNum_hit_WithinTran_aAICmPFC_ConAct) = 1; % within transient neurons
        end
        tempPairsNum_hit_WithinSust_aAICmPFC_ConAct = nnz(ismember(Pair_reg(:,1),2) & ismember(Pair_reg(:,2),1) & Pref_pair(:,2+bin)==Pref_pair(:,9+bin) & Pref_pair(:,2+bin)>0 & all(Pref_pair(:,3:6),2) & all(Pref_pair(:,10:13),2) & FRPairChains_go(:,1)>=FRRange(1) & FRPairChains_go(:,1)<FRRange(2) & FRPairChains_go(:,2)>=FRRange(1) & FRPairChains_go(:,2)<FRRange(2));
        tempFcPairsNum_hit_WithinSust_aAICmPFC_ConAct = nnz(ismember(reg_chain_go(:,1),2) & ismember(reg_chain_go(:,2),1) & reg_chain_go(:,3)>0 & pref_chain_go(:,2+bin)==pref_chain_go(:,9+bin) & pref_chain_go(:,2+bin)>0 & all(pref_chain_go(:,3:6),2) & all(pref_chain_go(:,10:13),2) & FR_conn_Chains_go(:,1)>=FRRange(1) & FR_conn_Chains_go(:,1)<FRRange(2) & FR_conn_Chains_go(:,2)>=FRRange(1) & FR_conn_Chains_go(:,2)<FRRange(2));
        IsSigFC_hit_WithinSust_ConAct{5,i}{1,bin} = zeros(tempPairsNum_hit_WithinSust_aAICmPFC_ConAct,1);
        if tempPairsNum_hit_WithinSust_aAICmPFC_ConAct > 0
            IsSigFC_hit_WithinSust_ConAct{5,i}{1,bin}(1:tempFcPairsNum_hit_WithinSust_aAICmPFC_ConAct) = 1; % within sustained neurons
        end
       %% Congruent inactive pairs //////only within transient neurons//////
        % mPFC, aAIC, mPFC-aAIC, and aAIC-mPFC
        tempPairsNum_hit_WithinTran_ConInact = nnz(ismember(Pair_reg(:,1),[1 2]) & ismember(Pair_reg(:,2),[1 2]) & max(Pref_pair(:,3:6),[],2)==max(Pref_pair(:,10:13),[],2) & any(Pref_pair(:,[bin+2 bin+9])<1,2) & ~all(Pref_pair(:,3:6),2) & ~all(Pref_pair(:,10:13),2) & FRPairChains_go(:,1)>=FRRange(1) & FRPairChains_go(:,1)<FRRange(2) & FRPairChains_go(:,2)>=FRRange(1) & FRPairChains_go(:,2)<FRRange(2));
        tempFcPairsNum_hit_WithinTran_ConInact = nnz(ismember(reg_chain_go(:,1),[1 2]) & ismember(reg_chain_go(:,2),[1 2]) & max(pref_chain_go(:,3:6),[],2)==max(pref_chain_go(:,10:13),[],2) & any(pref_chain_go(:,[bin+2 bin+9])<1,2) & ~all(pref_chain_go(:,3:6),2) & ~all(pref_chain_go(:,10:13),2) & FR_conn_Chains_go(:,1)>=FRRange(1) & FR_conn_Chains_go(:,1)<FRRange(2) & FR_conn_Chains_go(:,2)>=FRRange(1) & FR_conn_Chains_go(:,2)<FRRange(2));
        IsSigFC_hit_WithinTran_ConInact{1,i}{1,bin} = zeros(tempPairsNum_hit_WithinTran_ConInact,1);
        if tempPairsNum_hit_WithinTran_ConInact > 0
            IsSigFC_hit_WithinTran_ConInact{1,i}{1,bin}(1:tempFcPairsNum_hit_WithinTran_ConInact) = 1;
        end
        % mPFC
        tempPairsNum_hit_WithinTran_mPFC_ConInact = nnz(ismember(Pair_reg(:,1),1) & ismember(Pair_reg(:,2),1) & max(Pref_pair(:,3:6),[],2)==max(Pref_pair(:,10:13),[],2) & any(Pref_pair(:,[bin+2 bin+9])<1,2) & ~all(Pref_pair(:,3:6),2) & ~all(Pref_pair(:,10:13),2) & FRPairChains_go(:,1)>=FRRange(1) & FRPairChains_go(:,1)<FRRange(2) & FRPairChains_go(:,2)>=FRRange(1) & FRPairChains_go(:,2)<FRRange(2));
        tempFcPairsNum_hit_WithinTran_mPFC_ConInact = nnz(ismember(reg_chain_go(:,1),1) & ismember(reg_chain_go(:,2),1) & max(pref_chain_go(:,3:6),[],2)==max(pref_chain_go(:,10:13),[],2) & any(pref_chain_go(:,[bin+2 bin+9])<1,2) & ~all(pref_chain_go(:,3:6),2) & ~all(pref_chain_go(:,10:13),2) & FR_conn_Chains_go(:,1)>=FRRange(1) & FR_conn_Chains_go(:,1)<FRRange(2) & FR_conn_Chains_go(:,2)>=FRRange(1) & FR_conn_Chains_go(:,2)<FRRange(2));
        IsSigFC_hit_WithinTran_ConInact{2,i}{1,bin} = zeros(tempPairsNum_hit_WithinTran_mPFC_ConInact,1);
        if tempPairsNum_hit_WithinTran_mPFC_ConInact > 0
            IsSigFC_hit_WithinTran_ConInact{2,i}{1,bin}(1:tempFcPairsNum_hit_WithinTran_mPFC_ConInact) = 1;
        end
        % aAIC
        tempPairsNum_hit_WithinTran_aAIC_ConInact = nnz(ismember(Pair_reg(:,1),2) & ismember(Pair_reg(:,2),2) & max(Pref_pair(:,3:6),[],2)==max(Pref_pair(:,10:13),[],2) & any(Pref_pair(:,[bin+2 bin+9])<1,2) & ~all(Pref_pair(:,3:6),2) & ~all(Pref_pair(:,10:13),2) & FRPairChains_go(:,1)>=FRRange(1) & FRPairChains_go(:,1)<FRRange(2) & FRPairChains_go(:,2)>=FRRange(1) & FRPairChains_go(:,2)<FRRange(2));
        tempFcPairsNum_hit_WithinTran_aAIC_ConInact = nnz(ismember(reg_chain_go(:,1),2) & ismember(reg_chain_go(:,2),2) & max(pref_chain_go(:,3:6),[],2)==max(pref_chain_go(:,10:13),[],2) & any(pref_chain_go(:,[bin+2 bin+9])<1,2) & ~all(pref_chain_go(:,3:6),2) & ~all(pref_chain_go(:,10:13),2) & FR_conn_Chains_go(:,1)>=FRRange(1) & FR_conn_Chains_go(:,1)<FRRange(2) & FR_conn_Chains_go(:,2)>=FRRange(1) & FR_conn_Chains_go(:,2)<FRRange(2));
        IsSigFC_hit_WithinTran_ConInact{3,i}{1,bin} = zeros(tempPairsNum_hit_WithinTran_aAIC_ConInact,1);
        if tempPairsNum_hit_WithinTran_aAIC_ConInact > 0
            IsSigFC_hit_WithinTran_ConInact{3,i}{1,bin}(1:tempFcPairsNum_hit_WithinTran_aAIC_ConInact) = 1;
        end
        % mPFC-aAIC
        tempPairsNum_hit_WithinTran_mPFCaAIC_ConInact = nnz(ismember(Pair_reg(:,1),2) & ismember(Pair_reg(:,2),1) & max(Pref_pair(:,3:6),[],2)==max(Pref_pair(:,10:13),[],2) & any(Pref_pair(:,[bin+2 bin+9])<1,2) & ~all(Pref_pair(:,3:6),2) & ~all(Pref_pair(:,10:13),2) & FRPairChains_go(:,1)>=FRRange(1) & FRPairChains_go(:,1)<FRRange(2) & FRPairChains_go(:,2)>=FRRange(1) & FRPairChains_go(:,2)<FRRange(2));
        tempFcPairsNum_hit_WithinTran_mPFCaAIC_ConInact = nnz(ismember(reg_chain_go(:,1),2) & ismember(reg_chain_go(:,2),1) & reg_chain_go(:,3)<0 & max(pref_chain_go(:,3:6),[],2)==max(pref_chain_go(:,10:13),[],2) & any(pref_chain_go(:,[bin+2 bin+9])<1,2) & ~all(pref_chain_go(:,3:6),2) & ~all(pref_chain_go(:,10:13),2) & FR_conn_Chains_go(:,1)>=FRRange(1) & FR_conn_Chains_go(:,1)<FRRange(2) & FR_conn_Chains_go(:,2)>=FRRange(1) & FR_conn_Chains_go(:,2)<FRRange(2));
        IsSigFC_hit_WithinTran_ConInact{4,i}{1,bin} = zeros(tempPairsNum_hit_WithinTran_mPFCaAIC_ConInact,1);
        if tempPairsNum_hit_WithinTran_mPFCaAIC_ConInact > 0
            IsSigFC_hit_WithinTran_ConInact{4,i}{1,bin}(1:tempFcPairsNum_hit_WithinTran_mPFCaAIC_ConInact) = 1;
        end
        % aAIC-mPFC
        tempPairsNum_hit_WithinTran_aAICmPFC_ConInact = nnz(ismember(Pair_reg(:,1),2) & ismember(Pair_reg(:,2),1) & max(Pref_pair(:,3:6),[],2)==max(Pref_pair(:,10:13),[],2) & any(Pref_pair(:,[bin+2 bin+9])<1,2) & ~all(Pref_pair(:,3:6),2) & ~all(Pref_pair(:,10:13),2) & FRPairChains_go(:,1)>=FRRange(1) & FRPairChains_go(:,1)<FRRange(2) & FRPairChains_go(:,2)>=FRRange(1) & FRPairChains_go(:,2)<FRRange(2));
        tempFcPairsNum_hit_WithinTran_aAICmPFC_ConInact = nnz(ismember(reg_chain_go(:,1),2) & ismember(reg_chain_go(:,2),1) & reg_chain_go(:,3)>0 & max(pref_chain_go(:,3:6),[],2)==max(pref_chain_go(:,10:13),[],2) & any(pref_chain_go(:,[bin+2 bin+9])<1,2) & ~all(pref_chain_go(:,3:6),2) & ~all(pref_chain_go(:,10:13),2) & FR_conn_Chains_go(:,1)>=FRRange(1) & FR_conn_Chains_go(:,1)<FRRange(2) & FR_conn_Chains_go(:,2)>=FRRange(1) & FR_conn_Chains_go(:,2)<FRRange(2));
        IsSigFC_hit_WithinTran_ConInact{5,i}{1,bin} = zeros(tempPairsNum_hit_WithinTran_aAICmPFC_ConInact,1);
        if tempPairsNum_hit_WithinTran_aAICmPFC_ConInact > 0
            IsSigFC_hit_WithinTran_ConInact{5,i}{1,bin}(1:tempFcPairsNum_hit_WithinTran_aAICmPFC_ConInact) = 1;
        end
       %% Incongruent active pairs
        % mPFC, aAIC, mPFC-aAIC, and aAIC-mPFC
        tempPairsNum_hit_WithinTran_InconAct = nnz(ismember(Pair_reg(:,1),[1 2]) & ismember(Pair_reg(:,2),[1 2]) & Pref_pair(:,2+bin)>0 & Pref_pair(:,9+bin)>0 & Pref_pair(:,2+bin)~=Pref_pair(:,9+bin) & ~all(Pref_pair(:,3:6),2) & ~all(Pref_pair(:,10:13),2) & FRPairChains_go(:,1)>=FRRange(1) & FRPairChains_go(:,1)<FRRange(2) & FRPairChains_go(:,2)>=FRRange(1) & FRPairChains_go(:,2)<FRRange(2));
        tempFcPairsNum_hit_WithinTran_InconAct = nnz(ismember(reg_chain_go(:,1),[1 2]) & ismember(reg_chain_go(:,2),[1 2]) & pref_chain_go(:,2+bin)>0 & pref_chain_go(:,9+bin)>0 & pref_chain_go(:,2+bin)~=pref_chain_go(:,9+bin) & ~all(pref_chain_go(:,3:6),2) & ~all(pref_chain_go(:,10:13),2) & FR_conn_Chains_go(:,1)>=FRRange(1) & FR_conn_Chains_go(:,1)<FRRange(2) & FR_conn_Chains_go(:,2)>=FRRange(1) & FR_conn_Chains_go(:,2)<FRRange(2));
        IsSigFC_hit_WithinTran_InconAct{1,i}{1,bin} = zeros(tempPairsNum_hit_WithinTran_InconAct,1);
        if tempPairsNum_hit_WithinTran_InconAct > 0
            IsSigFC_hit_WithinTran_InconAct{1,i}{1,bin}(1:tempFcPairsNum_hit_WithinTran_InconAct) = 1; % within transient neurons
        end
        tempPairsNum_hit_WithinSust_InconAct = nnz(ismember(Pair_reg(:,1),[1 2]) & ismember(Pair_reg(:,2),[1 2]) & Pref_pair(:,2+bin)>0 & Pref_pair(:,9+bin)>0 & Pref_pair(:,2+bin)~=Pref_pair(:,9+bin) & all(Pref_pair(:,3:6),2) & all(Pref_pair(:,10:13),2) & FRPairChains_go(:,1)>=FRRange(1) & FRPairChains_go(:,1)<FRRange(2) & FRPairChains_go(:,2)>=FRRange(1) & FRPairChains_go(:,2)<FRRange(2));
        tempFcPairsNum_hit_WithinSust_InconAct = nnz(ismember(reg_chain_go(:,1),[1 2]) & ismember(reg_chain_go(:,2),[1 2]) & pref_chain_go(:,2+bin)>0 & pref_chain_go(:,9+bin)>0 & pref_chain_go(:,2+bin)~=pref_chain_go(:,9+bin) & all(pref_chain_go(:,3:6),2) & all(pref_chain_go(:,10:13),2) & FR_conn_Chains_go(:,1)>=FRRange(1) & FR_conn_Chains_go(:,1)<FRRange(2) & FR_conn_Chains_go(:,2)>=FRRange(1) & FR_conn_Chains_go(:,2)<FRRange(2));
        IsSigFC_hit_WithinSust_InconAct{1,i}{1,bin} = zeros(tempPairsNum_hit_WithinSust_InconAct,1);
        if tempPairsNum_hit_WithinSust_InconAct > 0
            IsSigFC_hit_WithinSust_InconAct{1,i}{1,bin}(1:tempFcPairsNum_hit_WithinSust_InconAct) = 1; % within sustained neurons
        end
        % mPFC
        tempPairsNum_hit_WithinTran_mPFC_InconAct = nnz(ismember(Pair_reg(:,1),1) & ismember(Pair_reg(:,2),1) & Pref_pair(:,2+bin)>0 & Pref_pair(:,9+bin)>0 & Pref_pair(:,2+bin)~=Pref_pair(:,9+bin) & ~all(Pref_pair(:,3:6),2) & ~all(Pref_pair(:,10:13),2) & FRPairChains_go(:,1)>=FRRange(1) & FRPairChains_go(:,1)<FRRange(2) & FRPairChains_go(:,2)>=FRRange(1) & FRPairChains_go(:,2)<FRRange(2));
        tempFcPairsNum_hit_WithinTran_mPFC_InconAct = nnz(ismember(reg_chain_go(:,1),1) & ismember(reg_chain_go(:,2),1) & pref_chain_go(:,2+bin)>0 & pref_chain_go(:,9+bin)>0 & pref_chain_go(:,2+bin)~=pref_chain_go(:,9+bin) & ~all(pref_chain_go(:,3:6),2) & ~all(pref_chain_go(:,10:13),2) & FR_conn_Chains_go(:,1)>=FRRange(1) & FR_conn_Chains_go(:,1)<FRRange(2) & FR_conn_Chains_go(:,2)>=FRRange(1) & FR_conn_Chains_go(:,2)<FRRange(2));
        IsSigFC_hit_WithinTran_InconAct{2,i}{1,bin} = zeros(tempPairsNum_hit_WithinTran_mPFC_InconAct,1);
        if tempPairsNum_hit_WithinTran_mPFC_InconAct > 0
            IsSigFC_hit_WithinTran_InconAct{2,i}{1,bin}(1:tempFcPairsNum_hit_WithinTran_mPFC_InconAct) = 1; % within transient neurons
        end
        tempPairsNum_hit_WithinSust_mPFC_InconAct = nnz(ismember(Pair_reg(:,1),1) & ismember(Pair_reg(:,2),1) & Pref_pair(:,2+bin)>0 & Pref_pair(:,9+bin)>0 & Pref_pair(:,2+bin)~=Pref_pair(:,9+bin) & all(Pref_pair(:,3:6),2) & all(Pref_pair(:,10:13),2) & FRPairChains_go(:,1)>=FRRange(1) & FRPairChains_go(:,1)<FRRange(2) & FRPairChains_go(:,2)>=FRRange(1) & FRPairChains_go(:,2)<FRRange(2));
        tempFcPairsNum_hit_WithinSust_mPFC_InconAct = nnz(ismember(reg_chain_go(:,1),1) & ismember(reg_chain_go(:,2),1) & pref_chain_go(:,2+bin)>0 & pref_chain_go(:,9+bin)>0 & pref_chain_go(:,2+bin)~=pref_chain_go(:,9+bin) & all(pref_chain_go(:,3:6),2) & all(pref_chain_go(:,10:13),2) & FR_conn_Chains_go(:,1)>=FRRange(1) & FR_conn_Chains_go(:,1)<FRRange(2) & FR_conn_Chains_go(:,2)>=FRRange(1) & FR_conn_Chains_go(:,2)<FRRange(2));
        IsSigFC_hit_WithinSust_InconAct{2,i}{1,bin} = zeros(tempPairsNum_hit_WithinSust_mPFC_InconAct,1);
        if tempPairsNum_hit_WithinSust_mPFC_InconAct > 0
            IsSigFC_hit_WithinSust_InconAct{2,i}{1,bin}(1:tempFcPairsNum_hit_WithinSust_mPFC_InconAct) = 1; % within sustained neurons
        end
        % aAIC
        tempPairsNum_hit_WithinTran_aAIC_InconAct = nnz(ismember(Pair_reg(:,1),2) & ismember(Pair_reg(:,2),2) & Pref_pair(:,2+bin)>0 & Pref_pair(:,9+bin)>0 & Pref_pair(:,2+bin)~=Pref_pair(:,9+bin) & ~all(Pref_pair(:,3:6),2) & ~all(Pref_pair(:,10:13),2) & FRPairChains_go(:,1)>=FRRange(1) & FRPairChains_go(:,1)<FRRange(2) & FRPairChains_go(:,2)>=FRRange(1) & FRPairChains_go(:,2)<FRRange(2));
        tempFcPairsNum_hit_WithinTran_aAIC_InconAct = nnz(ismember(reg_chain_go(:,1),2) & ismember(reg_chain_go(:,2),2) & pref_chain_go(:,2+bin)>0 & pref_chain_go(:,9+bin)>0 & pref_chain_go(:,2+bin)~=pref_chain_go(:,9+bin) & ~all(pref_chain_go(:,3:6),2) & ~all(pref_chain_go(:,10:13),2) & FR_conn_Chains_go(:,1)>=FRRange(1) & FR_conn_Chains_go(:,1)<FRRange(2) & FR_conn_Chains_go(:,2)>=FRRange(1) & FR_conn_Chains_go(:,2)<FRRange(2));
        IsSigFC_hit_WithinTran_InconAct{3,i}{1,bin} = zeros(tempPairsNum_hit_WithinTran_aAIC_InconAct,1);
        if tempPairsNum_hit_WithinTran_aAIC_InconAct > 0
            IsSigFC_hit_WithinTran_InconAct{3,i}{1,bin}(1:tempFcPairsNum_hit_WithinTran_aAIC_InconAct) = 1; % within transient neurons
        end
        tempPairsNum_hit_WithinSust_aAIC_InconAct = nnz(ismember(Pair_reg(:,1),2) & ismember(Pair_reg(:,2),2) & Pref_pair(:,2+bin)>0 & Pref_pair(:,9+bin)>0 & Pref_pair(:,2+bin)~=Pref_pair(:,9+bin) & all(Pref_pair(:,3:6),2) & all(Pref_pair(:,10:13),2) & FRPairChains_go(:,1)>=FRRange(1) & FRPairChains_go(:,1)<FRRange(2) & FRPairChains_go(:,2)>=FRRange(1) & FRPairChains_go(:,2)<FRRange(2));
        tempFcPairsNum_hit_WithinSust_aAIC_InconAct = nnz(ismember(reg_chain_go(:,1),2) & ismember(reg_chain_go(:,2),2) & pref_chain_go(:,2+bin)>0 & pref_chain_go(:,9+bin)>0 & pref_chain_go(:,2+bin)~=pref_chain_go(:,9+bin) & all(pref_chain_go(:,3:6),2) & all(pref_chain_go(:,10:13),2) & FR_conn_Chains_go(:,1)>=FRRange(1) & FR_conn_Chains_go(:,1)<FRRange(2) & FR_conn_Chains_go(:,2)>=FRRange(1) & FR_conn_Chains_go(:,2)<FRRange(2));
        IsSigFC_hit_WithinSust_InconAct{3,i}{1,bin} = zeros(tempPairsNum_hit_WithinSust_aAIC_InconAct,1);
        if tempPairsNum_hit_WithinSust_aAIC_InconAct > 0
            IsSigFC_hit_WithinSust_InconAct{3,i}{1,bin}(1:tempFcPairsNum_hit_WithinSust_aAIC_InconAct) = 1; % within sustained neurons
        end
        % mPFC-aAIC
        tempPairsNum_hit_WithinTran_mPFCaAIC_InconAct = nnz(ismember(Pair_reg(:,1),2) & ismember(Pair_reg(:,2),1) & Pref_pair(:,2+bin)>0 & Pref_pair(:,9+bin)>0 & Pref_pair(:,2+bin)~=Pref_pair(:,9+bin) & ~all(Pref_pair(:,3:6),2) & ~all(Pref_pair(:,10:13),2) & FRPairChains_go(:,1)>=FRRange(1) & FRPairChains_go(:,1)<FRRange(2) & FRPairChains_go(:,2)>=FRRange(1) & FRPairChains_go(:,2)<FRRange(2));
        tempFcPairsNum_hit_WithinTran_mPFCaAIC_InconAct = nnz(ismember(reg_chain_go(:,1),2) & ismember(reg_chain_go(:,2),1) & reg_chain_go(:,3)<0 & pref_chain_go(:,2+bin)>0 & pref_chain_go(:,9+bin)>0 & pref_chain_go(:,2+bin)~=pref_chain_go(:,9+bin) & ~all(pref_chain_go(:,3:6),2) & ~all(pref_chain_go(:,10:13),2) & FR_conn_Chains_go(:,1)>=FRRange(1) & FR_conn_Chains_go(:,1)<FRRange(2) & FR_conn_Chains_go(:,2)>=FRRange(1) & FR_conn_Chains_go(:,2)<FRRange(2));
        IsSigFC_hit_WithinTran_InconAct{4,i}{1,bin} = zeros(tempPairsNum_hit_WithinTran_mPFCaAIC_InconAct,1);
        if tempPairsNum_hit_WithinTran_mPFCaAIC_InconAct > 0
            IsSigFC_hit_WithinTran_InconAct{4,i}{1,bin}(1:tempFcPairsNum_hit_WithinTran_mPFCaAIC_InconAct) = 1; % within transient neurons
        end
        tempPairsNum_hit_WithinSust_mPFCaAIC_InconAct = nnz(ismember(Pair_reg(:,1),2) & ismember(Pair_reg(:,2),1) & Pref_pair(:,2+bin)>0 & Pref_pair(:,9+bin)>0 & Pref_pair(:,2+bin)~=Pref_pair(:,9+bin) & all(Pref_pair(:,3:6),2) & all(Pref_pair(:,10:13),2) & FRPairChains_go(:,1)>=FRRange(1) & FRPairChains_go(:,1)<FRRange(2) & FRPairChains_go(:,2)>=FRRange(1) & FRPairChains_go(:,2)<FRRange(2));
        tempFcPairsNum_hit_WithinSust_mPFCaAIC_InconAct = nnz(ismember(reg_chain_go(:,1),2) & ismember(reg_chain_go(:,2),1) & reg_chain_go(:,3)<0 & pref_chain_go(:,2+bin)>0 & pref_chain_go(:,9+bin)>0 & pref_chain_go(:,2+bin)~=pref_chain_go(:,9+bin) & all(pref_chain_go(:,3:6),2) & all(pref_chain_go(:,10:13),2) & FR_conn_Chains_go(:,1)>=FRRange(1) & FR_conn_Chains_go(:,1)<FRRange(2) & FR_conn_Chains_go(:,2)>=FRRange(1) & FR_conn_Chains_go(:,2)<FRRange(2));
        IsSigFC_hit_WithinSust_InconAct{4,i}{1,bin} = zeros(tempPairsNum_hit_WithinSust_mPFCaAIC_InconAct,1);
        if tempPairsNum_hit_WithinSust_mPFCaAIC_InconAct > 0
            IsSigFC_hit_WithinSust_InconAct{4,i}{1,bin}(1:tempFcPairsNum_hit_WithinSust_mPFCaAIC_InconAct) = 1; % within sustained neurons
        end
        % aAIC-mPFC
        tempPairsNum_hit_WithinTran_aAICmPFC_InconAct = nnz(ismember(Pair_reg(:,1),2) & ismember(Pair_reg(:,2),1) & Pref_pair(:,2+bin)>0 & Pref_pair(:,9+bin)>0 & Pref_pair(:,2+bin)~=Pref_pair(:,9+bin) & ~all(Pref_pair(:,3:6),2) & ~all(Pref_pair(:,10:13),2) & FRPairChains_go(:,1)>=FRRange(1) & FRPairChains_go(:,1)<FRRange(2) & FRPairChains_go(:,2)>=FRRange(1) & FRPairChains_go(:,2)<FRRange(2));
        tempFcPairsNum_hit_WithinTran_aAICmPFC_InconAct = nnz(ismember(reg_chain_go(:,1),2) & ismember(reg_chain_go(:,2),1) & reg_chain_go(:,3)>0 & pref_chain_go(:,2+bin)>0 & pref_chain_go(:,9+bin)>0 & pref_chain_go(:,2+bin)~=pref_chain_go(:,9+bin) & ~all(pref_chain_go(:,3:6),2) & ~all(pref_chain_go(:,10:13),2) & FR_conn_Chains_go(:,1)>=FRRange(1) & FR_conn_Chains_go(:,1)<FRRange(2) & FR_conn_Chains_go(:,2)>=FRRange(1) & FR_conn_Chains_go(:,2)<FRRange(2));
        IsSigFC_hit_WithinTran_InconAct{5,i}{1,bin} = zeros(tempPairsNum_hit_WithinTran_aAICmPFC_InconAct,1);
        if tempPairsNum_hit_WithinTran_aAICmPFC_InconAct > 0
            IsSigFC_hit_WithinTran_InconAct{5,i}{1,bin}(1:tempFcPairsNum_hit_WithinTran_aAICmPFC_InconAct) = 1; % within transient neurons
        end
        tempPairsNum_hit_WithinSust_aAICmPFC_InconAct = nnz(ismember(Pair_reg(:,1),2) & ismember(Pair_reg(:,2),1) & Pref_pair(:,2+bin)>0 & Pref_pair(:,9+bin)>0 & Pref_pair(:,2+bin)~=Pref_pair(:,9+bin) & all(Pref_pair(:,3:6),2) & all(Pref_pair(:,10:13),2) & FRPairChains_go(:,1)>=FRRange(1) & FRPairChains_go(:,1)<FRRange(2) & FRPairChains_go(:,2)>=FRRange(1) & FRPairChains_go(:,2)<FRRange(2));
        tempFcPairsNum_hit_WithinSust_aAICmPFC_InconAct = nnz(ismember(reg_chain_go(:,1),2) & ismember(reg_chain_go(:,2),1) & reg_chain_go(:,3)>0 & pref_chain_go(:,2+bin)>0 & pref_chain_go(:,9+bin)>0 & pref_chain_go(:,2+bin)~=pref_chain_go(:,9+bin) & all(pref_chain_go(:,3:6),2) & all(pref_chain_go(:,10:13),2) & FR_conn_Chains_go(:,1)>=FRRange(1) & FR_conn_Chains_go(:,1)<FRRange(2) & FR_conn_Chains_go(:,2)>=FRRange(1) & FR_conn_Chains_go(:,2)<FRRange(2));
        IsSigFC_hit_WithinSust_InconAct{5,i}{1,bin} = zeros(tempPairsNum_hit_WithinSust_aAICmPFC_InconAct,1);
        if tempPairsNum_hit_WithinSust_aAICmPFC_InconAct > 0
            IsSigFC_hit_WithinSust_InconAct{5,i}{1,bin}(1:tempFcPairsNum_hit_WithinSust_aAICmPFC_InconAct) = 1; % within sustained neurons
        end
        
       %% Incongruent inactive pairs
        % mPFC, aAIC, mPFC-aAIC, and aAIC-mPFC
        tempPairsNum_hit_WithinTran_InconInact = nnz(ismember(Pair_reg(:,1),[1 2]) & ismember(Pair_reg(:,2),[1 2]) & max(Pref_pair(:,3:6),[],2)~=max(Pref_pair(:,10:13),[],2) & max(Pref_pair(:,3:6),[],2)>0 & max(Pref_pair(:,10:13),[],2)>0 & any(Pref_pair(:,[bin+2 bin+9])<1,2) & ~all(Pref_pair(:,3:6),2) & ~all(Pref_pair(:,10:13),2) & FRPairChains_go(:,1)>=FRRange(1) & FRPairChains_go(:,1)<FRRange(2) & FRPairChains_go(:,2)>=FRRange(1) & FRPairChains_go(:,2)<FRRange(2));
        tempFcPairsNum_hit_WithinTran_InconInact = nnz(ismember(reg_chain_go(:,1),[1 2]) & ismember(reg_chain_go(:,2),[1 2]) & max(pref_chain_go(:,3:6),[],2)~=max(pref_chain_go(:,10:13),[],2) & max(pref_chain_go(:,3:6),[],2)>0 & max(pref_chain_go(:,10:13),[],2)>0 & any(pref_chain_go(:,[bin+2 bin+9])<1,2) & ~all(pref_chain_go(:,3:6),2) & ~all(pref_chain_go(:,10:13),2) & FR_conn_Chains_go(:,1)>=FRRange(1) & FR_conn_Chains_go(:,1)<FRRange(2) & FR_conn_Chains_go(:,2)>=FRRange(1) & FR_conn_Chains_go(:,2)<FRRange(2));
        IsSigFC_hit_WithinTran_InconInact{1,i}{1,bin} = zeros(tempPairsNum_hit_WithinTran_InconInact,1);
        if tempPairsNum_hit_WithinTran_InconInact > 0
            IsSigFC_hit_WithinTran_InconInact{1,i}{1,bin}(1:tempFcPairsNum_hit_WithinTran_InconInact) = 1;
        end
        % mPFC
        tempPairsNum_hit_WithinTran_mPFC_InconInact = nnz(ismember(Pair_reg(:,1),1) & ismember(Pair_reg(:,2),1) & max(Pref_pair(:,3:6),[],2)~=max(Pref_pair(:,10:13),[],2) & max(Pref_pair(:,3:6),[],2)>0 & max(Pref_pair(:,10:13),[],2)>0 & any(Pref_pair(:,[bin+2 bin+9])<1,2) & ~all(Pref_pair(:,3:6),2) & ~all(Pref_pair(:,10:13),2) & FRPairChains_go(:,1)>=FRRange(1) & FRPairChains_go(:,1)<FRRange(2) & FRPairChains_go(:,2)>=FRRange(1) & FRPairChains_go(:,2)<FRRange(2));
        tempFcPairsNum_hit_WithinTran_mPFC_InconInact = nnz(ismember(reg_chain_go(:,1),1) & ismember(reg_chain_go(:,2),1) & max(pref_chain_go(:,3:6),[],2)~=max(pref_chain_go(:,10:13),[],2) & max(pref_chain_go(:,3:6),[],2)>0 & max(pref_chain_go(:,10:13),[],2)>0 & any(pref_chain_go(:,[bin+2 bin+9])<1,2) & ~all(pref_chain_go(:,3:6),2) & ~all(pref_chain_go(:,10:13),2) & FR_conn_Chains_go(:,1)>=FRRange(1) & FR_conn_Chains_go(:,1)<FRRange(2) & FR_conn_Chains_go(:,2)>=FRRange(1) & FR_conn_Chains_go(:,2)<FRRange(2));
        IsSigFC_hit_WithinTran_InconInact{2,i}{1,bin} = zeros(tempPairsNum_hit_WithinTran_mPFC_InconInact,1);
        if tempPairsNum_hit_WithinTran_mPFC_InconInact > 0
            IsSigFC_hit_WithinTran_InconInact{2,i}{1,bin}(1:tempFcPairsNum_hit_WithinTran_mPFC_InconInact) = 1;
        end
        % aAIC
        tempPairsNum_hit_WithinTran_aAIC_InconInact = nnz(ismember(Pair_reg(:,1),2) & ismember(Pair_reg(:,2),2) & max(Pref_pair(:,3:6),[],2)~=max(Pref_pair(:,10:13),[],2) & max(Pref_pair(:,3:6),[],2)>0 & max(Pref_pair(:,10:13),[],2)>0 & any(Pref_pair(:,[bin+2 bin+9])<1,2) & ~all(Pref_pair(:,3:6),2) & ~all(Pref_pair(:,10:13),2) & FRPairChains_go(:,1)>=FRRange(1) & FRPairChains_go(:,1)<FRRange(2) & FRPairChains_go(:,2)>=FRRange(1) & FRPairChains_go(:,2)<FRRange(2));
        tempFcPairsNum_hit_WithinTran_aAIC_InconInact = nnz(ismember(reg_chain_go(:,1),2) & ismember(reg_chain_go(:,2),2) & max(pref_chain_go(:,3:6),[],2)~=max(pref_chain_go(:,10:13),[],2) & max(pref_chain_go(:,3:6),[],2)>0 & max(pref_chain_go(:,10:13),[],2)>0 & any(pref_chain_go(:,[bin+2 bin+9])<1,2) & ~all(pref_chain_go(:,3:6),2) & ~all(pref_chain_go(:,10:13),2) & FR_conn_Chains_go(:,1)>=FRRange(1) & FR_conn_Chains_go(:,1)<FRRange(2) & FR_conn_Chains_go(:,2)>=FRRange(1) & FR_conn_Chains_go(:,2)<FRRange(2));
        IsSigFC_hit_WithinTran_InconInact{3,i}{1,bin} = zeros(tempPairsNum_hit_WithinTran_aAIC_InconInact,1);
        if tempPairsNum_hit_WithinTran_aAIC_InconInact > 0
            IsSigFC_hit_WithinTran_InconInact{3,i}{1,bin}(1:tempFcPairsNum_hit_WithinTran_aAIC_InconInact) = 1;
        end
        % mPFC-aAIC
        tempPairsNum_hit_WithinTran_mPFCaAIC_InconInact = nnz(ismember(Pair_reg(:,1),2) & ismember(Pair_reg(:,2),1) & max(Pref_pair(:,3:6),[],2)~=max(Pref_pair(:,10:13),[],2) & max(Pref_pair(:,3:6),[],2)>0 & max(Pref_pair(:,10:13),[],2)>0 & any(Pref_pair(:,[bin+2 bin+9])<1,2) & ~all(Pref_pair(:,3:6),2) & ~all(Pref_pair(:,10:13),2) & FRPairChains_go(:,1)>=FRRange(1) & FRPairChains_go(:,1)<FRRange(2) & FRPairChains_go(:,2)>=FRRange(1) & FRPairChains_go(:,2)<FRRange(2));
        tempFcPairsNum_hit_WithinTran_mPFCaAIC_InconInact = nnz(ismember(reg_chain_go(:,1),2) & ismember(reg_chain_go(:,2),1) & reg_chain_go(:,3)<0 & max(pref_chain_go(:,3:6),[],2)~=max(pref_chain_go(:,10:13),[],2) & max(pref_chain_go(:,3:6),[],2)>0 & max(pref_chain_go(:,10:13),[],2)>0 & any(pref_chain_go(:,[bin+2 bin+9])<1,2) & ~all(pref_chain_go(:,3:6),2) & ~all(pref_chain_go(:,10:13),2) & FR_conn_Chains_go(:,1)>=FRRange(1) & FR_conn_Chains_go(:,1)<FRRange(2) & FR_conn_Chains_go(:,2)>=FRRange(1) & FR_conn_Chains_go(:,2)<FRRange(2));
        IsSigFC_hit_WithinTran_InconInact{4,i}{1,bin} = zeros(tempPairsNum_hit_WithinTran_mPFCaAIC_InconInact,1);
        if tempPairsNum_hit_WithinTran_mPFCaAIC_InconInact > 0
            IsSigFC_hit_WithinTran_InconInact{4,i}{1,bin}(1:tempFcPairsNum_hit_WithinTran_mPFCaAIC_InconInact) = 1;
        end
        % aAIC-mPFC
        tempPairsNum_hit_WithinTran_aAICmPFC_InconInact = nnz(ismember(Pair_reg(:,1),2) & ismember(Pair_reg(:,2),1) & max(Pref_pair(:,3:6),[],2)~=max(Pref_pair(:,10:13),[],2) & max(Pref_pair(:,3:6),[],2)>0 & max(Pref_pair(:,10:13),[],2)>0 & any(Pref_pair(:,[bin+2 bin+9])<1,2) & ~all(Pref_pair(:,3:6),2) & ~all(Pref_pair(:,10:13),2) & FRPairChains_go(:,1)>=FRRange(1) & FRPairChains_go(:,1)<FRRange(2) & FRPairChains_go(:,2)>=FRRange(1) & FRPairChains_go(:,2)<FRRange(2));
        tempFcPairsNum_hit_WithinTran_aAICmPFC_InconInact = nnz(ismember(reg_chain_go(:,1),2) & ismember(reg_chain_go(:,2),1) & reg_chain_go(:,3)>0 & max(pref_chain_go(:,3:6),[],2)~=max(pref_chain_go(:,10:13),[],2) & max(pref_chain_go(:,3:6),[],2)>0 & max(pref_chain_go(:,10:13),[],2)>0 & any(pref_chain_go(:,[bin+2 bin+9])<1,2) & ~all(pref_chain_go(:,3:6),2) & ~all(pref_chain_go(:,10:13),2) & FR_conn_Chains_go(:,1)>=FRRange(1) & FR_conn_Chains_go(:,1)<FRRange(2) & FR_conn_Chains_go(:,2)>=FRRange(1) & FR_conn_Chains_go(:,2)<FRRange(2));
        IsSigFC_hit_WithinTran_InconInact{5,i}{1,bin} = zeros(tempPairsNum_hit_WithinTran_aAICmPFC_InconInact,1);
        if tempPairsNum_hit_WithinTran_aAICmPFC_InconInact > 0
            IsSigFC_hit_WithinTran_InconInact{5,i}{1,bin}(1:tempFcPairsNum_hit_WithinTran_aAICmPFC_InconInact) = 1;
        end
    end
end

%% Save variables
save(['Region-specific FC density_within transient or sustained neurons_within FR ' num2str(TarFR(1)) ' to ' num2str(TarFR(1)+5) ' to ' num2str(TarFR(end)) '-' Group],'IsSigFC_hit_WithinTran_ConAct','IsSigFC_hit_WithinSust_ConAct','IsSigFC_hit_WithinTran_ConInact','IsSigFC_hit_WithinTran_InconAct','IsSigFC_hit_WithinSust_InconAct','IsSigFC_hit_WithinTran_InconInact','-v7.3');

