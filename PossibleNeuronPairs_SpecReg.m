%% Region-specific possible neuronal pairs

clear; clc; close all;

%% Assignment
Group = 'CtrlGroup';
DelayBinNum = 4;

%% Various types of neuronal pairs
Pairs_ConAct = cell(4,DelayBinNum); % each row represents each region, 1: all regions; 2: mPFC; 3: aAIC; 4: between mPFC and aAIC
Pairs_ConInact = cell(4,DelayBinNum);
Pairs_InconAct = cell(4,DelayBinNum);
Pairs_InconInact = cell(4,DelayBinNum);
Pairs_Nonmem = cell(4,DelayBinNum);

%% Pairs of non-memory neurons
for bin = 1:DelayBinNum % four 1-sec bins of delay
    load(sprintf('FCofNonmemoryNeurons_PropertyPopulation_%d_%d_%s.mat',bin,bin+1,Group));
    % mPFC, aAIC, and between mPFC and aAIC
    Pairs_Nonmem{1,bin} = Pair_chain(ismember(Pair_reg(:,1),[1 2]) & ismember(Pair_reg(:,2),[1 2]),:);
    % mPFC
    Pairs_Nonmem{2,bin} = Pair_chain(ismember(Pair_reg(:,1),1) & ismember(Pair_reg(:,2),1),:);
    % aAIC
    Pairs_Nonmem{3,bin} = Pair_chain(ismember(Pair_reg(:,1),2) & ismember(Pair_reg(:,2),2),:);
    % mPFC-aAIC
    Pairs_Nonmem{4,bin} = Pair_chain(ismember(Pair_reg(:,1),2) & ismember(Pair_reg(:,2),1),:);
    clear Pair_reg reg_chain_go FRPairChains_go FR_conn_Chains_go
end

%% Pairs of memory neurons
for bin = 1:DelayBinNum % four 1-second bins of delay
    load(sprintf('FCofMemoryNeurons_PropertyPopulation_%d_%d_%s.mat',bin,bin+1,Group)); 
    %% Congruent active pairs
    % mPFC, aAIC, and between mPFC and aAIC
    Pairs_ConAct{1,bin}{1,1} = Pair_chain(ismember(Pair_reg(:,1),[1 2]) & ismember(Pair_reg(:,2),[1 2]) & Pref_pair(:,2+bin)==Pref_pair(:,9+bin) & Pref_pair(:,2+bin)==1,:);
    Pairs_ConAct{1,bin}{2,1} = Pair_chain(ismember(Pair_reg(:,1),[1 2]) & ismember(Pair_reg(:,2),[1 2]) & Pref_pair(:,2+bin)==Pref_pair(:,9+bin) & Pref_pair(:,2+bin)==2,:);
    % mPFC
    Pairs_ConAct{2,bin}{1,1} = Pair_chain(ismember(Pair_reg(:,1),1) & ismember(Pair_reg(:,2),1) & Pref_pair(:,2+bin)==Pref_pair(:,9+bin) & Pref_pair(:,2+bin)==1,:);
    Pairs_ConAct{2,bin}{2,1} = Pair_chain(ismember(Pair_reg(:,1),1) & ismember(Pair_reg(:,2),1) & Pref_pair(:,2+bin)==Pref_pair(:,9+bin) & Pref_pair(:,2+bin)==2,:);
    % aAIC
    Pairs_ConAct{3,bin}{1,1} = Pair_chain(ismember(Pair_reg(:,1),2) & ismember(Pair_reg(:,2),2) & Pref_pair(:,2+bin)==Pref_pair(:,9+bin) & Pref_pair(:,2+bin)==1,:);
    Pairs_ConAct{3,bin}{2,1} = Pair_chain(ismember(Pair_reg(:,1),2) & ismember(Pair_reg(:,2),2) & Pref_pair(:,2+bin)==Pref_pair(:,9+bin) & Pref_pair(:,2+bin)==2,:);
    % between mPFC and aAIC
    Pairs_ConAct{4,bin}{1,1} = Pair_chain(ismember(Pair_reg(:,1),2) & ismember(Pair_reg(:,2),1) & Pref_pair(:,2+bin)==Pref_pair(:,9+bin) & Pref_pair(:,2+bin)==1,:);
    Pairs_ConAct{4,bin}{2,1} = Pair_chain(ismember(Pair_reg(:,1),2) & ismember(Pair_reg(:,2),1) & Pref_pair(:,2+bin)==Pref_pair(:,9+bin) & Pref_pair(:,2+bin)==2,:);
    
    %% Congruent inactive pairs
    % mPFC, aAIC, and between mPFC and aAIC
    Pairs_ConInact{1,bin}{1,1} = Pair_chain(ismember(Pair_reg(:,1),[1 2]) & ismember(Pair_reg(:,2),[1 2]) & max(Pref_pair(:,3:6),[],2)==max(Pref_pair(:,10:13),[],2) & max(Pref_pair(:,3:6),[],2)==1  & any(Pref_pair(:,[bin+2 bin+9])<1,2),:);
    Pairs_ConInact{1,bin}{2,1} = Pair_chain(ismember(Pair_reg(:,1),[1 2]) & ismember(Pair_reg(:,2),[1 2]) & max(Pref_pair(:,3:6),[],2)==max(Pref_pair(:,10:13),[],2) & max(Pref_pair(:,3:6),[],2)==2  & any(Pref_pair(:,[bin+2 bin+9])<1,2),:);
    % mPFC
    Pairs_ConInact{2,bin}{1,1} = Pair_chain(ismember(Pair_reg(:,1),1) & ismember(Pair_reg(:,2),1) & max(Pref_pair(:,3:6),[],2)==max(Pref_pair(:,10:13),[],2) & max(Pref_pair(:,3:6),[],2)==1  & any(Pref_pair(:,[bin+2 bin+9])<1,2),:);
    Pairs_ConInact{2,bin}{2,1} = Pair_chain(ismember(Pair_reg(:,1),1) & ismember(Pair_reg(:,2),1) & max(Pref_pair(:,3:6),[],2)==max(Pref_pair(:,10:13),[],2) & max(Pref_pair(:,3:6),[],2)==2  & any(Pref_pair(:,[bin+2 bin+9])<1,2),:);
    % aAIC
    Pairs_ConInact{3,bin}{1,1} = Pair_chain(ismember(Pair_reg(:,1),2) & ismember(Pair_reg(:,2),2) & max(Pref_pair(:,3:6),[],2)==max(Pref_pair(:,10:13),[],2) & max(Pref_pair(:,3:6),[],2)==1  & any(Pref_pair(:,[bin+2 bin+9])<1,2),:);
    Pairs_ConInact{3,bin}{2,1} = Pair_chain(ismember(Pair_reg(:,1),2) & ismember(Pair_reg(:,2),2) & max(Pref_pair(:,3:6),[],2)==max(Pref_pair(:,10:13),[],2) & max(Pref_pair(:,3:6),[],2)==2  & any(Pref_pair(:,[bin+2 bin+9])<1,2),:);
    % between mPFC and aAIC
    Pairs_ConInact{4,bin}{1,1} = Pair_chain(ismember(Pair_reg(:,1),2) & ismember(Pair_reg(:,2),1) & max(Pref_pair(:,3:6),[],2)==max(Pref_pair(:,10:13),[],2) & max(Pref_pair(:,3:6),[],2)==1  & any(Pref_pair(:,[bin+2 bin+9])<1,2),:);
    Pairs_ConInact{4,bin}{2,1} = Pair_chain(ismember(Pair_reg(:,1),2) & ismember(Pair_reg(:,2),1) & max(Pref_pair(:,3:6),[],2)==max(Pref_pair(:,10:13),[],2) & max(Pref_pair(:,3:6),[],2)==2  & any(Pref_pair(:,[bin+2 bin+9])<1,2),:);
    
    %% Incongruent active pairs
    % mPFC, aAIC, and between mPFC and aAIC
    Pairs_InconAct{1,bin} = Pair_chain(ismember(Pair_reg(:,1),[1 2]) & ismember(Pair_reg(:,2),[1 2]) & Pref_pair(:,2+bin)>0 & Pref_pair(:,9+bin)>0 & Pref_pair(:,2+bin)~=Pref_pair(:,9+bin),:);
    % mPFC
    Pairs_InconAct{2,bin} = Pair_chain(ismember(Pair_reg(:,1),1) & ismember(Pair_reg(:,2),1) & Pref_pair(:,2+bin)>0 & Pref_pair(:,9+bin)>0 & Pref_pair(:,2+bin)~=Pref_pair(:,9+bin),:);
    % aAIC
    Pairs_InconAct{3,bin} = Pair_chain(ismember(Pair_reg(:,1),2) & ismember(Pair_reg(:,2),2) & Pref_pair(:,2+bin)>0 & Pref_pair(:,9+bin)>0 & Pref_pair(:,2+bin)~=Pref_pair(:,9+bin),:);
    % between mPFC and aAIC
    Pairs_InconAct{4,bin} = Pair_chain(ismember(Pair_reg(:,1),2) & ismember(Pair_reg(:,2),1) & Pref_pair(:,2+bin)>0 & Pref_pair(:,9+bin)>0 & Pref_pair(:,2+bin)~=Pref_pair(:,9+bin),:);

    %% Incongruent inactive pairs
    % mPFC, aAIC, and between mPFC and aAIC
    Pairs_InconInact{1,bin} = Pair_chain(ismember(Pair_reg(:,1),[1 2]) & ismember(Pair_reg(:,2),[1 2]) & max(Pref_pair(:,3:6),[],2)~=max(Pref_pair(:,10:13),[],2) & max(Pref_pair(:,3:6),[],2)>0 & max(Pref_pair(:,10:13),[],2)>0 & any(Pref_pair(:,[bin+2 bin+9])<1,2),:);
    % mPFC
    Pairs_InconInact{2,bin} = Pair_chain(ismember(Pair_reg(:,1),1) & ismember(Pair_reg(:,2),1) & max(Pref_pair(:,3:6),[],2)~=max(Pref_pair(:,10:13),[],2) & max(Pref_pair(:,3:6),[],2)>0 & max(Pref_pair(:,10:13),[],2)>0 & any(Pref_pair(:,[bin+2 bin+9])<1,2),:);
    % aAIC
    Pairs_InconInact{3,bin} = Pair_chain(ismember(Pair_reg(:,1),2) & ismember(Pair_reg(:,2),2) & max(Pref_pair(:,3:6),[],2)~=max(Pref_pair(:,10:13),[],2) & max(Pref_pair(:,3:6),[],2)>0 & max(Pref_pair(:,10:13),[],2)>0 & any(Pref_pair(:,[bin+2 bin+9])<1,2),:);
    % between mPFC and aAIC
    Pairs_InconInact{4,bin} = Pair_chain(ismember(Pair_reg(:,1),2) & ismember(Pair_reg(:,2),1) & max(Pref_pair(:,3:6),[],2)~=max(Pref_pair(:,10:13),[],2) & max(Pref_pair(:,3:6),[],2)>0 & max(Pref_pair(:,10:13),[],2)>0 & any(Pref_pair(:,[bin+2 bin+9])<1,2),:);
end

%% Save variables
save(['Region-specific possible pairs-' Group],'Pairs_ConAct','Pairs_ConInact','Pairs_InconAct','Pairs_InconInact','Pairs_Nonmem','-v7.3');

