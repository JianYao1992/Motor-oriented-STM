%% FC within transient neurons (4 types of pairs) or sustained neurons (2 types)
% //////session based//////

clear; clc; close all;

%% Assignment
Group = 'CtrlGroup';
reg = 'aAIC-mPFC';
if strcmp(reg,'Within mPFC')
    TarReg1 = 1; TarReg2 = 1;
elseif strcmp(reg,'Within aAIC')
    TarReg1 = 2; TarReg2 = 2;
elseif strcmp(reg,'mPFC-aAIC') || strcmp(reg,'aAIC-mPFC')
    TarReg1 = 2; TarReg2 = 1;
end
PairsNumthres = 12;
BinNumCriteria = 1;

%% Summarize FC across delay bins
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

%% FC density of pairs of transient neurons
% congruent active
PairNum_ConAct = [];
FcPairNum_hit_ConAct = [];
% congruent inactive
PairNum_ConInact = [];
FcPairNum_hit_ConInact = [];
% incongruent active
PairNum_InconAct = [];
FcPairNum_hit_InconAct = [];
% incongruent inactive
PairNum_InconInact = [];
FcPairNum_hit_InconInact = [];
for bin = 1:4 % four 1-sec bin of delay
    for iFile = 1:48
        %% Number of pairs
        PairNum_ConAct(iFile,bin) = nnz(reg_pairs{bin}(:,1)==TarReg1 & reg_pairs{bin}(:,2)==TarReg2 & pref_pairs{bin}(:,2+bin)==pref_pairs{bin}(:,9+bin) & pref_pairs{bin}(:,2+bin)>0 & pair_mice{bin}==iFile & ~all(pref_pairs{bin}(:,3:6),2) & ~all(pref_pairs{bin}(:,10:13),2));
        PairNum_ConInact(iFile,bin) = nnz(reg_pairs{bin}(:,1)==TarReg1 & reg_pairs{bin}(:,2)==TarReg2 & max(pref_pairs{bin}(:,3:6),[],2)==max(pref_pairs{bin}(:,10:13),[],2) & any(pref_pairs{bin}(:,[bin+2 bin+9])<1,2) & pair_mice{bin}==iFile & ~all(pref_pairs{bin}(:,3:6),2) & ~all(pref_pairs{bin}(:,10:13),2));
        PairNum_InconAct(iFile,bin) = nnz(reg_pairs{bin}(:,1)==TarReg1 & reg_pairs{bin}(:,2)==TarReg2 & pref_pairs{bin}(:,2+bin)>0 & pref_pairs{bin}(:,9+bin)>0 & pref_pairs{bin}(:,2+bin)~=pref_pairs{bin}(:,9+bin) & pair_mice{bin}==iFile & ~all(pref_pairs{bin}(:,3:6),2) & ~all(pref_pairs{bin}(:,10:13),2));
        PairNum_InconInact(iFile,bin) = nnz(reg_pairs{bin}(:,1)==TarReg1 & reg_pairs{bin}(:,2)==TarReg2 & max(pref_pairs{bin}(:,3:6),[],2)~=max(pref_pairs{bin}(:,10:13),[],2) & max(pref_pairs{bin}(:,3:6),[],2)>0 & max(pref_pairs{bin}(:,10:13),[],2)>0 & any(pref_pairs{bin}(:,[bin+2 bin+9])<1,2) & pair_mice{bin}==iFile & ~all(pref_pairs{bin}(:,3:6),2) & ~all(pref_pairs{bin}(:,10:13),2));
        
        %% Number of FC pairs in hit trials
        if strcmp(reg,'mPFC-aAIC')
            FcPairNum_hit_ConAct(iFile,bin) = nnz(reg_chains_go{bin}(:,1)==TarReg1 & reg_chains_go{bin}(:,2)==TarReg2 & reg_chains_go{bin}(:,3)<0 & pref_chains_go{bin}(:,2+bin)==pref_chains_go{bin}(:,9+bin) & pref_chains_go{bin}(:,2+bin)>0 & mouse_chains_go{bin}==iFile & ~all(pref_chains_go{bin}(:,3:6),2) & ~all(pref_chains_go{bin}(:,10:13),2));
            FcPairNum_hit_ConInact(iFile,bin) = nnz(reg_chains_go{bin}(:,1)==TarReg1 & reg_chains_go{bin}(:,2)==TarReg2 & reg_chains_go{bin}(:,3)<0 & max(pref_chains_go{bin}(:,3:6),[],2) == max(pref_chains_go{bin}(:,10:13),[],2) & any(pref_chains_go{bin}(:,[bin+2 bin+9])<1,2) & mouse_chains_go{bin}==iFile & ~all(pref_chains_go{bin}(:,3:6),2) & ~all(pref_chains_go{bin}(:,10:13),2));
            FcPairNum_hit_InconAct(iFile,bin) = nnz(reg_chains_go{bin}(:,1)==TarReg1 & reg_chains_go{bin}(:,2)==TarReg2 & reg_chains_go{bin}(:,3)<0 & pref_chains_go{bin}(:,2+bin)>0 & pref_chains_go{bin}(:,9+bin)>0 & pref_chains_go{bin}(:,2+bin)~=pref_chains_go{bin}(:,9+bin) & mouse_chains_go{bin}==iFile & ~all(pref_chains_go{bin}(:,3:6),2) & ~all(pref_chains_go{bin}(:,10:13),2));
            FcPairNum_hit_InconInact(iFile,bin) = nnz(reg_chains_go{bin}(:,1)==TarReg1 & reg_chains_go{bin}(:,2)==TarReg2 & reg_chains_go{bin}(:,3)<0 & max(pref_chains_go{bin}(:,3:6),[],2)~=max(pref_chains_go{bin}(:,10:13),[],2) & max(pref_chains_go{bin}(:,3:6),[],2)>0 & max(pref_chains_go{bin}(:,10:13),[],2)>0 & any(pref_chains_go{bin}(:,[bin+2 bin+9])<1,2) & mouse_chains_go{bin}==iFile & ~all(pref_chains_go{bin}(:,3:6),2) & ~all(pref_chains_go{bin}(:,10:13),2));
        elseif strcmp(reg,'aAIC-mPFC')
             FcPairNum_hit_ConAct(iFile,bin) = nnz(reg_chains_go{bin}(:,1)==TarReg1 & reg_chains_go{bin}(:,2)==TarReg2 & reg_chains_go{bin}(:,3)>0 & pref_chains_go{bin}(:,2+bin)==pref_chains_go{bin}(:,9+bin) & pref_chains_go{bin}(:,2+bin)>0 & mouse_chains_go{bin}==iFile & ~all(pref_chains_go{bin}(:,3:6),2) & ~all(pref_chains_go{bin}(:,10:13),2));
            FcPairNum_hit_ConInact(iFile,bin) = nnz(reg_chains_go{bin}(:,1)==TarReg1 & reg_chains_go{bin}(:,2)==TarReg2 & reg_chains_go{bin}(:,3)>0 & max(pref_chains_go{bin}(:,3:6),[],2) == max(pref_chains_go{bin}(:,10:13),[],2) & any(pref_chains_go{bin}(:,[bin+2 bin+9])<1,2) & mouse_chains_go{bin}==iFile & ~all(pref_chains_go{bin}(:,3:6),2) & ~all(pref_chains_go{bin}(:,10:13),2));
            FcPairNum_hit_InconAct(iFile,bin) = nnz(reg_chains_go{bin}(:,1)==TarReg1 & reg_chains_go{bin}(:,2)==TarReg2 & reg_chains_go{bin}(:,3)>0 & pref_chains_go{bin}(:,2+bin)>0 & pref_chains_go{bin}(:,9+bin)>0 & pref_chains_go{bin}(:,2+bin)~=pref_chains_go{bin}(:,9+bin) & mouse_chains_go{bin}==iFile & ~all(pref_chains_go{bin}(:,3:6),2) & ~all(pref_chains_go{bin}(:,10:13),2));
            FcPairNum_hit_InconInact(iFile,bin) = nnz(reg_chains_go{bin}(:,1)==TarReg1 & reg_chains_go{bin}(:,2)==TarReg2 & reg_chains_go{bin}(:,3)>0 & max(pref_chains_go{bin}(:,3:6),[],2)~=max(pref_chains_go{bin}(:,10:13),[],2) & max(pref_chains_go{bin}(:,3:6),[],2)>0 & max(pref_chains_go{bin}(:,10:13),[],2)>0 & any(pref_chains_go{bin}(:,[bin+2 bin+9])<1,2) & mouse_chains_go{bin}==iFile & ~all(pref_chains_go{bin}(:,3:6),2) & ~all(pref_chains_go{bin}(:,10:13),2));
        else
            FcPairNum_hit_ConAct(iFile,bin) = nnz(reg_chains_go{bin}(:,1)==TarReg1 & reg_chains_go{bin}(:,2)==TarReg2 & pref_chains_go{bin}(:,2+bin)==pref_chains_go{bin}(:,9+bin) & pref_chains_go{bin}(:,2+bin)>0 & mouse_chains_go{bin}==iFile & ~all(pref_chains_go{bin}(:,3:6),2) & ~all(pref_chains_go{bin}(:,10:13),2));
            FcPairNum_hit_ConInact(iFile,bin) = nnz(reg_chains_go{bin}(:,1)==TarReg1 & reg_chains_go{bin}(:,2)==TarReg2 & max(pref_chains_go{bin}(:,3:6),[],2) == max(pref_chains_go{bin}(:,10:13),[],2) & any(pref_chains_go{bin}(:,[bin+2 bin+9])<1,2) & mouse_chains_go{bin}==iFile & ~all(pref_chains_go{bin}(:,3:6),2) & ~all(pref_chains_go{bin}(:,10:13),2));
            FcPairNum_hit_InconAct(iFile,bin) = nnz(reg_chains_go{bin}(:,1)==TarReg1 & reg_chains_go{bin}(:,2)==TarReg2 & pref_chains_go{bin}(:,2+bin)>0 & pref_chains_go{bin}(:,9+bin)>0 & pref_chains_go{bin}(:,2+bin)~=pref_chains_go{bin}(:,9+bin) & mouse_chains_go{bin}==iFile & ~all(pref_chains_go{bin}(:,3:6),2) & ~all(pref_chains_go{bin}(:,10:13),2));
            FcPairNum_hit_InconInact(iFile,bin) = nnz(reg_chains_go{bin}(:,1)==TarReg1 & reg_chains_go{bin}(:,2)==TarReg2 & max(pref_chains_go{bin}(:,3:6),[],2)~=max(pref_chains_go{bin}(:,10:13),[],2) & max(pref_chains_go{bin}(:,3:6),[],2)>0 & max(pref_chains_go{bin}(:,10:13),[],2)>0 & any(pref_chains_go{bin}(:,[bin+2 bin+9])<1,2) & mouse_chains_go{bin}==iFile & ~all(pref_chains_go{bin}(:,3:6),2) & ~all(pref_chains_go{bin}(:,10:13),2));
        end
    end
end
FcDensity_Tran_hit_ConAct = CalculateSessionBasedFcDensity(PairNum_ConAct,FcPairNum_hit_ConAct,PairsNumthres,BinNumCriteria);
FcDensity_Tran_hit_ConInact = CalculateSessionBasedFcDensity(PairNum_ConInact,FcPairNum_hit_ConInact,PairsNumthres,BinNumCriteria);
FcDensity_Tran_hit_InconAct = CalculateSessionBasedFcDensity(PairNum_InconAct,FcPairNum_hit_InconAct,PairsNumthres,BinNumCriteria);
FcDensity_Tran_hit_InconInact = CalculateSessionBasedFcDensity(PairNum_InconInact,FcPairNum_hit_InconInact,PairsNumthres,BinNumCriteria);
FcDensity_Tran_hit = [{FcDensity_Tran_hit_InconInact} {FcDensity_Tran_hit_InconAct} {FcDensity_Tran_hit_ConInact} {FcDensity_Tran_hit_ConAct}];

%% FC density of pairs of sustained neurons
% congruent active
PairNum_ConAct = [];
FcPairNum_hit_ConAct = [];
% incongruent active
PairNum_InconAct = [];
FcPairNum_hit_InconAct = [];
for bin = 1:4 % four 1-sec bin of delay
    for iFile = 1:48
        %% Number of pairs
        PairNum_ConAct(iFile,bin) = nnz(reg_pairs{bin}(:,1)==TarReg1 & reg_pairs{bin}(:,2)==TarReg2 & pref_pairs{bin}(:,2+bin)==pref_pairs{bin}(:,9+bin) & pref_pairs{bin}(:,2+bin)>0 & pair_mice{bin}==iFile & all(pref_pairs{bin}(:,3:6),2) & all(pref_pairs{bin}(:,10:13),2));
        PairNum_InconAct(iFile,bin) = nnz(reg_pairs{bin}(:,1)==TarReg1 & reg_pairs{bin}(:,2)==TarReg2 & pref_pairs{bin}(:,2+bin)>0 & pref_pairs{bin}(:,9+bin)>0 & pref_pairs{bin}(:,2+bin)~=pref_pairs{bin}(:,9+bin) & pair_mice{bin}==iFile & all(pref_pairs{bin}(:,3:6),2) & all(pref_pairs{bin}(:,10:13),2));
        
        %% Number of FC pairs in hit trials
        if strcmp(reg,'mPFC-aAIC')
            FcPairNum_hit_ConAct(iFile,bin) = nnz(reg_chains_go{bin}(:,1)==TarReg1 & reg_chains_go{bin}(:,2)==TarReg2 & reg_chains_go{bin}(:,3)<0 & pref_chains_go{bin}(:,2+bin)==pref_chains_go{bin}(:,9+bin) & pref_chains_go{bin}(:,2+bin)>0 & mouse_chains_go{bin}==iFile & all(pref_chains_go{bin}(:,3:6),2) & all(pref_chains_go{bin}(:,10:13),2));
            FcPairNum_hit_InconAct(iFile,bin) = nnz(reg_chains_go{bin}(:,1)==TarReg1 & reg_chains_go{bin}(:,2)==TarReg2 & reg_chains_go{bin}(:,3)<0 & pref_chains_go{bin}(:,2+bin)>0 & pref_chains_go{bin}(:,9+bin)>0 & pref_chains_go{bin}(:,2+bin)~=pref_chains_go{bin}(:,9+bin) & mouse_chains_go{bin}==iFile & all(pref_chains_go{bin}(:,3:6),2) & all(pref_chains_go{bin}(:,10:13),2));
        elseif strcmp(reg,'aAIC-mPFC')
             FcPairNum_hit_ConAct(iFile,bin) = nnz(reg_chains_go{bin}(:,1)==TarReg1 & reg_chains_go{bin}(:,2)==TarReg2 & reg_chains_go{bin}(:,3)>0 & pref_chains_go{bin}(:,2+bin)==pref_chains_go{bin}(:,9+bin) & pref_chains_go{bin}(:,2+bin)>0 & mouse_chains_go{bin}==iFile & all(pref_chains_go{bin}(:,3:6),2) & all(pref_chains_go{bin}(:,10:13),2));
            FcPairNum_hit_InconAct(iFile,bin) = nnz(reg_chains_go{bin}(:,1)==TarReg1 & reg_chains_go{bin}(:,2)==TarReg2 & reg_chains_go{bin}(:,3)>0 & pref_chains_go{bin}(:,2+bin)>0 & pref_chains_go{bin}(:,9+bin)>0 & pref_chains_go{bin}(:,2+bin)~=pref_chains_go{bin}(:,9+bin) & mouse_chains_go{bin}==iFile & all(pref_chains_go{bin}(:,3:6),2) & all(pref_chains_go{bin}(:,10:13),2));
        else
            FcPairNum_hit_ConAct(iFile,bin) = nnz(reg_chains_go{bin}(:,1)==TarReg1 & reg_chains_go{bin}(:,2)==TarReg2 & pref_chains_go{bin}(:,2+bin)==pref_chains_go{bin}(:,9+bin) & pref_chains_go{bin}(:,2+bin)>0 & mouse_chains_go{bin}==iFile & all(pref_chains_go{bin}(:,3:6),2) & all(pref_chains_go{bin}(:,10:13),2));
            FcPairNum_hit_InconAct(iFile,bin) = nnz(reg_chains_go{bin}(:,1)==TarReg1 & reg_chains_go{bin}(:,2)==TarReg2 & pref_chains_go{bin}(:,2+bin)>0 & pref_chains_go{bin}(:,9+bin)>0 & pref_chains_go{bin}(:,2+bin)~=pref_chains_go{bin}(:,9+bin) & mouse_chains_go{bin}==iFile & all(pref_chains_go{bin}(:,3:6),2) & all(pref_chains_go{bin}(:,10:13),2));
        end
    end
end
FcDensity_Sust_hit_ConAct = CalculateSessionBasedFcDensity(PairNum_ConAct,FcPairNum_hit_ConAct,2,BinNumCriteria);
FcDensity_Sust_hit_InconAct = CalculateSessionBasedFcDensity(PairNum_InconAct,FcPairNum_hit_InconAct,2,BinNumCriteria);
FcDensity_Sust_hit = [{FcDensity_Sust_hit_InconAct} {FcDensity_Sust_hit_ConAct}];

%% Save result
save(sprintf('SessionBased FC density within transient or sustained neurons_hit_%s',reg),'FcDensity_Tran_hit','FcDensity_Sust_hit','-v7.3');

