%% FC within congruent active, congruent inactive, incongruent active, incongruent inactive, non-memory (both neurons are non-memory in all delay bins)
% //////over delay bins//////
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
PairsNumthres = 10;

%% FC density of pairs of non-memory neurons
FcDensity_hit_Nonmem = cell(0);
for bin = 1:4 % four 1-sec bins of delay
    load(sprintf('FCofNonmemoryNeurons_PropertyPopulation_%d_%d_%s.mat',bin,bin+1,Group));
    tempFcPairNum_hit_Nonmem = [];
    tempPairNum_Nonmem = [];
    for iFile = 1:48 % session ID
        tempPairNum_Nonmem(iFile,1) = nnz(Pair_mouse==iFile & Pair_reg(:,1)==TarReg1 & Pair_reg(:,2)==TarReg2);
        if strcmp(reg,'mPFC-aAIC')
            tempFcPairNum_hit_Nonmem(iFile,1) = nnz(mouse_chain_go==iFile& reg_chain_go(:,1)==TarReg1 & reg_chain_go(:,2)==TarReg2 & reg_chain_go(:,3)<0);
        elseif strcmp(reg,'aAIC-mPFC')
            tempFcPairNum_hit_Nonmem(iFile,1) = nnz(mouse_chain_go==iFile & reg_chain_go(:,1)==TarReg1 & reg_chain_go(:,2)==TarReg2 & reg_chain_go(:,3)>0);
        else
            tempFcPairNum_hit_Nonmem(iFile,1) = nnz(mouse_chain_go==iFile & reg_chain_go(:,1)==TarReg1 & reg_chain_go(:,2)==TarReg2);
        end
    end
    AboveThres = find(tempPairNum_Nonmem>=PairsNumthres);
    tempFcPairNum_hit_Nonmem = tempFcPairNum_hit_Nonmem(AboveThres);
    tempPairNum_Nonmem = tempPairNum_Nonmem(AboveThres);
    FcDensity_hit_Nonmem{1,end+1} = tempFcPairNum_hit_Nonmem./tempPairNum_Nonmem;
    clear Pair_reg Pair_chain Pref_pair reg_chain_go conn_chain_go pref_chain_go mouse_chain_go
end

%% FC density of pairs of memory neurons
for bin = 1:4 % four 1-sec bin of delay
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
FcDensity_hit_ConAct = cell(0);
FcDensity_hit_ConInact = cell(0);
FcDensity_hit_InconAct = cell(0);
FcDensity_hit_InconInact = cell(0);
for bin = 1:4 % four 1-second bin of delay
    tempPairNum_ConAct = [];
    tempPairNum_ConInact = [];
    tempPairNum_InconAct = [];
    tempPairNum_InconInact = [];
    tempFcPairNum_hit_ConAct = [];
    tempFcPairNum_hit_ConInact = [];
    tempFcPairNum_hit_InconAct = [];
    tempFcPairNum_hit_InconInact = [];
    for iFile = 1:48
        %% Number of pairs
        tempPairNum_ConAct(iFile,1) = nnz(reg_pairs{bin}(:,1)==TarReg1 & reg_pairs{bin}(:,2)==TarReg2 & pref_pairs{bin}(:,2+bin)==pref_pairs{bin}(:,9+bin) & pref_pairs{bin}(:,2+bin)>0 & pair_mice{bin}==iFile);
        tempPairNum_ConInact(iFile,1) = nnz(reg_pairs{bin}(:,1)==TarReg1 & reg_pairs{bin}(:,2)==TarReg2 & max(pref_pairs{bin}(:,3:6),[],2)==max(pref_pairs{bin}(:,10:13),[],2) & any(pref_pairs{bin}(:,[bin+2 bin+9])<1,2) & pair_mice{bin}==iFile);
        tempPairNum_InconAct(iFile,1) = nnz(reg_pairs{bin}(:,1)==TarReg1 & reg_pairs{bin}(:,2)==TarReg2 & pref_pairs{bin}(:,2+bin)>0 & pref_pairs{bin}(:,9+bin)>0 & pref_pairs{bin}(:,2+bin)~=pref_pairs{bin}(:,9+bin) & pair_mice{bin}==iFile);
        tempPairNum_InconInact(iFile,1) = nnz(reg_pairs{bin}(:,1)==TarReg1 & reg_pairs{bin}(:,2)==TarReg2 & max(pref_pairs{bin}(:,3:6),[],2)~=max(pref_pairs{bin}(:,10:13),[],2) & max(pref_pairs{bin}(:,3:6),[],2)>0 & max(pref_pairs{bin}(:,10:13),[],2)>0 & any(pref_pairs{bin}(:,[bin+2 bin+9])<1,2) & pair_mice{bin}==iFile);
        
        %% Number of FC pairs in hit trials
        if strcmp(reg,'mPFC-aAIC')
            tempFcPairNum_hit_ConAct(iFile,1) = nnz(reg_chains_go{bin}(:,1)==TarReg1 & reg_chains_go{bin}(:,2)==TarReg2 & reg_chains_go{bin}(:,3)<0 & pref_chains_go{bin}(:,2+bin)==pref_chains_go{bin}(:,9+bin) & pref_chains_go{bin}(:,2+bin)>0 & mouse_chains_go{bin}==iFile);
            tempFcPairNum_hit_ConInact(iFile,1) = nnz(reg_chains_go{bin}(:,1)==TarReg1 & reg_chains_go{bin}(:,2)==TarReg2 & reg_chains_go{bin}(:,3)<0 & max(pref_chains_go{bin}(:,3:6),[],2) == max(pref_chains_go{bin}(:,10:13),[],2) & any(pref_chains_go{bin}(:,[bin+2 bin+9])<1,2) & mouse_chains_go{bin}==iFile);
            tempFcPairNum_hit_InconAct(iFile,1) = nnz(reg_chains_go{bin}(:,1)==TarReg1 & reg_chains_go{bin}(:,2)==TarReg2 & reg_chains_go{bin}(:,3)<0 & pref_chains_go{bin}(:,2+bin)>0 & pref_chains_go{bin}(:,9+bin)>0 & pref_chains_go{bin}(:,2+bin)~=pref_chains_go{bin}(:,9+bin) & mouse_chains_go{bin}==iFile);
            tempFcPairNum_hit_InconInact(iFile,1) = nnz(reg_chains_go{bin}(:,1)==TarReg1 & reg_chains_go{bin}(:,2)==TarReg2 & reg_chains_go{bin}(:,3)<0 & max(pref_chains_go{bin}(:,3:6),[],2)~=max(pref_chains_go{bin}(:,10:13),[],2) & max(pref_chains_go{bin}(:,3:6),[],2)>0 & max(pref_chains_go{bin}(:,10:13),[],2)>0 & any(pref_chains_go{bin}(:,[bin+2 bin+9])<1,2) & mouse_chains_go{bin}==iFile);
        elseif strcmp(reg,'aAIC-mPFC')
            tempFcPairNum_hit_ConAct(iFile,1) = nnz(reg_chains_go{bin}(:,1)==TarReg1 & reg_chains_go{bin}(:,2)==TarReg2 & reg_chains_go{bin}(:,3)>0 & pref_chains_go{bin}(:,2+bin)==pref_chains_go{bin}(:,9+bin) & pref_chains_go{bin}(:,2+bin)>0 & mouse_chains_go{bin}==iFile);
            tempFcPairNum_hit_ConInact(iFile,1) = nnz(reg_chains_go{bin}(:,1)==TarReg1 & reg_chains_go{bin}(:,2)==TarReg2 & reg_chains_go{bin}(:,3)>0 & max(pref_chains_go{bin}(:,3:6),[],2) == max(pref_chains_go{bin}(:,10:13),[],2) & any(pref_chains_go{bin}(:,[bin+2 bin+9])<1,2) & mouse_chains_go{bin}==iFile);
            tempFcPairNum_hit_InconAct(iFile,1) = nnz(reg_chains_go{bin}(:,1)==TarReg1 & reg_chains_go{bin}(:,2)==TarReg2 & reg_chains_go{bin}(:,3)>0 & pref_chains_go{bin}(:,2+bin)>0 & pref_chains_go{bin}(:,9+bin)>0 & pref_chains_go{bin}(:,2+bin)~=pref_chains_go{bin}(:,9+bin) & mouse_chains_go{bin}==iFile);
            tempFcPairNum_hit_InconInact(iFile,1) = nnz(reg_chains_go{bin}(:,1)==TarReg1 & reg_chains_go{bin}(:,2)==TarReg2 & reg_chains_go{bin}(:,3)>0 & max(pref_chains_go{bin}(:,3:6),[],2)~=max(pref_chains_go{bin}(:,10:13),[],2) & max(pref_chains_go{bin}(:,3:6),[],2)>0 & max(pref_chains_go{bin}(:,10:13),[],2)>0 & any(pref_chains_go{bin}(:,[bin+2 bin+9])<1,2) & mouse_chains_go{bin}==iFile);
        else
            tempFcPairNum_hit_ConAct(iFile,1) = nnz(reg_chains_go{bin}(:,1)==TarReg1 & reg_chains_go{bin}(:,2)==TarReg2 & pref_chains_go{bin}(:,2+bin)==pref_chains_go{bin}(:,9+bin) & pref_chains_go{bin}(:,2+bin)>0 & mouse_chains_go{bin}==iFile);
            tempFcPairNum_hit_ConInact(iFile,1) = nnz(reg_chains_go{bin}(:,1)==TarReg1 & reg_chains_go{bin}(:,2)==TarReg2 & max(pref_chains_go{bin}(:,3:6),[],2) == max(pref_chains_go{bin}(:,10:13),[],2) & any(pref_chains_go{bin}(:,[bin+2 bin+9])<1,2) & mouse_chains_go{bin}==iFile);
            tempFcPairNum_hit_InconAct(iFile,1) = nnz(reg_chains_go{bin}(:,1)==TarReg1 & reg_chains_go{bin}(:,2)==TarReg2 & pref_chains_go{bin}(:,2+bin)>0 & pref_chains_go{bin}(:,9+bin)>0 & pref_chains_go{bin}(:,2+bin)~=pref_chains_go{bin}(:,9+bin) & mouse_chains_go{bin}==iFile);
            tempFcPairNum_hit_InconInact(iFile,1) = nnz(reg_chains_go{bin}(:,1)==TarReg1 & reg_chains_go{bin}(:,2)==TarReg2 & max(pref_chains_go{bin}(:,3:6),[],2)~=max(pref_chains_go{bin}(:,10:13),[],2) & max(pref_chains_go{bin}(:,3:6),[],2)>0 & max(pref_chains_go{bin}(:,10:13),[],2)>0 & any(pref_chains_go{bin}(:,[bin+2 bin+9])<1,2) & mouse_chains_go{bin}==iFile);
        end
    end
    % congruent active pairs
    AboveThres = find(tempPairNum_ConAct>=PairsNumthres);
    tempFcPairNum_hit_ConAct = tempFcPairNum_hit_ConAct(AboveThres);
    tempPairNum_ConAct = tempPairNum_ConAct(AboveThres);
    FcDensity_hit_ConAct{1,end+1} = tempFcPairNum_hit_ConAct./tempPairNum_ConAct;
    % congruent inactive pairs
    AboveThres = find(tempPairNum_ConInact>=PairsNumthres);
    tempFcPairNum_hit_ConInact = tempFcPairNum_hit_ConInact(AboveThres);
    tempPairNum_ConInact = tempPairNum_ConInact(AboveThres);
    FcDensity_hit_ConInact{1,end+1} = tempFcPairNum_hit_ConInact./tempPairNum_ConInact;
    % incongruent active pairs
    AboveThres = find(tempPairNum_InconAct>=PairsNumthres);
    tempFcPairNum_hit_InconAct = tempFcPairNum_hit_InconAct(AboveThres);
    tempPairNum_InconAct = tempPairNum_InconAct(AboveThres);
    FcDensity_hit_InconAct{1,end+1} = tempFcPairNum_hit_InconAct./tempPairNum_InconAct;
    % incongruent inactive pairs
    AboveThres = find(tempPairNum_InconInact>=PairsNumthres);
    tempFcPairNum_hit_InconInact = tempFcPairNum_hit_InconInact(AboveThres);
    tempPairNum_InconInact = tempPairNum_InconInact(AboveThres);
    FcDensity_hit_InconInact{1,end+1} = tempFcPairNum_hit_InconInact./tempPairNum_InconInact;
end

%% Save result
FcDensity_hit = [{FcDensity_hit_Nonmem} {FcDensity_hit_InconInact} {FcDensity_hit_InconAct} {FcDensity_hit_ConInact} {FcDensity_hit_ConAct}];
save(sprintf('SessionBased FC density_OverDelayBins_5types_hit_%s',reg),'FcDensity_hit','-v7.3');

