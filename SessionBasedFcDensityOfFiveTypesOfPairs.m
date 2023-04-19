%% FC within congruent active, congruent inactive, incongruent active, incongruent inactive, non-memory (both neurons are non-memory in all delay bins)
% /// session based ///

clear; clc; close all;

%% Assignment
reg = 'Within mPFC';
if strcmp(reg,'Within mPFC')
    TarReg1 = 1; TarReg2 = 1;
elseif strcmp(reg,'Within aAIC')
    TarReg1 = 2; TarReg2 = 2;
elseif strcmp(reg,'mPFC-aAIC')
    TarReg1 = 2; TarReg2 = 1;
end
pairrange = 'all';
PairsNumthres = 10;
BinNumCriteria = 2;

%% FC density of pairs of non-memory neurons
PairNum_Nonmem = [];
FcPairNum_hit_Nonmem = [];
for bin = 1:4 % four 1-second bin of delay
    load(sprintf('FCofNonmemoryNeurons_PropertyPopulation_%d_%d_%s.mat',bin,bin+1,Group));
    load(sprintf('PairsFRModulation_nonsel_conn_chain_delay4s_%d_%d.mat',bin,bin+1));
    for iFile = 1:48 % session ID
        PairNum_Nonmem(iFile,bin) = nnz(Pair_mouse==iFile & Pair_reg(:,1)==TarReg1 & Pair_reg(:,2)==TarReg2);
        FcPairNum_hit_Nonmem(iFile,bin) = nnz(mouse_chain_go==iFile & reg_chain_go(:,1)==TarReg1 & reg_chain_go(:,2)==TarReg2);
    end
    clear Pair_reg Pair_chain Pref_pair reg_chain_go conn_chain_go pref_chain_go mouse_chain_go
end
FcDensity_hit_Nonmem = CalculateSessionBasedFcDensity(PairNum_Nonmem,FcPairNum_hit_Nonmem,PairsNumthres,BinNumCriteria);

%% FC density of pairs of memory neurons
% congruent active
PairNum_ConAct = [];
TrialNum_hit_ConAct = [];
FcPairNum_hit_ConAct = [];
% congruent inactive
PairNum_ConInact = [];
TrialNum_hit_ConInact = [];
FcPairNum_hit_ConInact = [];
% incongruent active
PairNum_InconAct = [];
TrialNum_hit_InconAct = [];
FcPairNum_hit_InconAct = [];
% incongruent inactive
PairNum_InconInact = [];
TrialNum_hit_IncongInact = [];
FcPairNum_hit_InconInact = [];
for bin = 1:4 % four 1-second bin of delay
    load(sprintf('FCofMemoryNeurons_PropertyPopulation_%d_%d_%s.mat',bin,bin+1,Group));
    % pair
    reg_pairs{bin} = Pair_reg;
    pref_pairs{bin} = Pref_pair;
    pair_mice{bin} = Pair_mouse;
    % FC pair
    reg_chains_go{bin} = reg_chain_go;
    pref_chains_go{bin} = pref_chain_go;
    mouse_chains_go{bin} = mouse_chain_go;
    clear Pair_reg Pref_pair Pair_mouse conn_chain_go pref_chain_go mouse_chain_go
end
for bin = 1:4 % four 1-second bin of delay
    load(sprintf('PairsFRModulation_selec_conn_chain_delay4s_%d_%d.mat',bin,bin+1));
    for iFile = 1:48
        %% Number of pairs
        PairNum_ConAct(iFile,bin) = nnz(reg_pairs{bin}(:,1)==TarReg1 & reg_pairs{bin}(:,2)==TarReg2 & pref_pairs{bin}(:,2+bin)==pref_pairs{bin}(:,9+bin) & pref_pairs{bin}(:,2+bin)>0 & pair_mice{bin}==iFile);
        PairNum_ConInact(iFile,bin) = nnz(reg_pairs{bin}(:,1)==TarReg1 & reg_pairs{bin}(:,2)==TarReg2 & max(pref_pairs{bin}(:,3:6),[],2)==max(pref_pairs{bin}(:,10:13),[],2) & any(pref_pairs{bin}(:,[bin+2 bin+9])<1,2) & pair_mice{bin}==iFile);
        PairNum_InconAct(iFile,bin) = nnz(reg_pairs{bin}(:,1)==TarReg1 & reg_pairs{bin}(:,2)==TarReg2 & pref_pairs{bin}(:,2+bin)>0 & pref_pairs{bin}(:,9+bin)>0 & pref_pairs{bin}(:,2+bin)~=pref_pairs{bin}(:,9+bin) & pair_mice{bin}==iFile);
        PairNum_InconInact(iFile,bin) = nnz(reg_pairs{bin}(:,1)==TarReg1 & reg_pairs{bin}(:,2)==TarReg2 & max(pref_pairs{bin}(:,3:6),[],2)~=max(pref_pairs{bin}(:,10:13),[],2) & max(pref_pairs{bin}(:,3:6),[],2)>0 & max(pref_pairs{bin}(:,10:13),[],2)>0 & any(pref_pairs{bin}(:,[bin+2 bin+9])<1,2) & pair_mice{bin}==iFile);
        
        %% Number of FC pairs in hit trials
        FcPairNum_hit_ConAct(iFile,bin) = nnz(reg_chains_go{bin}(:,1)==TarReg1 & reg_chains_go{bin}(:,2)==TarReg2 & pref_chains_go{bin}(:,2+bin)==pref_chains_go{bin}(:,9+bin) & pref_chains_go{bin}(:,2+bin)>0 & mouse_chains_go{bin}==iFile);
        FcPairNum_hit_ConInact(iFile,bin) = nnz(reg_chains_go{bin}(:,1)==TarReg1 & reg_chains_go{bin}(:,2)==TarReg2 & max(pref_chains_go{bin}(:,3:6),[],2) == max(pref_chains_go{bin}(:,10:13),[],2) & any(pref_chains_go{bin}(:,[bin+2 bin+9])<1,2) & mouse_chains_go{bin}==iFile);
        FcPairNum_hit_InconAct(iFile,bin) = nnz(reg_chains_go{bin}(:,1)==TarReg1 & reg_chains_go{bin}(:,2)==TarReg2 & pref_chains_go{bin}(:,2+bin)>0 & pref_chains_go{bin}(:,9+bin)>0 & pref_chains_go{bin}(:,2+bin)~=pref_chains_go{bin}(:,9+bin) & mouse_chains_go{bin}==iFile);
        FcPairNum_hit_InconInact(iFile,bin) = nnz(reg_chains_go{bin}(:,1)==TarReg1 & reg_chains_go{bin}(:,2)==TarReg2 & max(pref_chains_go{bin}(:,3:6),[],2)~=max(pref_chains_go{bin}(:,10:13),[],2) & max(pref_chains_go{bin}(:,3:6),[],2)>0 & max(pref_chains_go{bin}(:,10:13),[],2)>0 & any(pref_chains_go{bin}(:,[bin+2 bin+9])<1,2) & mouse_chains_go{bin}==iFile);
    end
end
FcDensity_hit_ConAct = CalculateSessionBasedFcDensity(PairNum_ConAct,FcPairNum_hit_ConAct,PairsNumthres,BinNumCriteria);
FcDensity_hit_ConInact = CalculateSessionBasedFcDensity(PairNum_ConInact,FcPairNum_hit_ConInact,PairsNumthres,BinNumCriteria);
FcDensity_hit_InconAct = CalculateSessionBasedFcDensity(PairNum_InconAct,FcPairNum_hit_InconAct,PairsNumthres,BinNumCriteria);
FcDensity_hit_InconInact = CalculateSessionBasedFcDensity(PairNum_InconInact,FcPairNum_hit_InconInact,PairsNumthres,BinNumCriteria);

%% Save result
FcDensity_hit = [{FcDensity_hit_Nonmem} {FcDensity_hit_InconInact} {FcDensity_hit_InconAct} {FcDensity_hit_ConInact} {FcDensity_hit_ConAct}];
save(sprintf('PairsNumThres %d BinNumThres %d SessionBased FC density_5types_hit_%s',PairsNumthres,BinNumCriteria,reg),'FcDensity_hit','-v7.3');

