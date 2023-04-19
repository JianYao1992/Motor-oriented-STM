%% FC density of congruent active, congruent inactive, incongruent active, incongruent inactive, non-memory (both neurons are non-coding in all delay bins)
% /// second based ///
% /// FR modulation ///

clear; clc; close all;

%% Assignment
Group = 'CtrlGroup';

%% FC density of pairs of non-memory neurons /// II, DD, and ID ///
PairsNum_hit_II_Nonmem = []; FcPairsNum_hit_II_Nonmem = []; % I: increased FR; D: decreased FR
PairsNum_hit_DD_Nonmem = []; FcPairsNum_hit_DD_Nonmem = [];
PairsNum_hit_ID_Nonmem = []; FcPairsNum_hit_ID_Nonmem = [];

for bin = 1:4 % four 1-second bin of delay
    load(sprintf('FCofNonmemoryNeurons_PropertyPopulation_%d_%d_%s.mat',bin,bin+1,Group));
    load(sprintf('PairsFRModulation_nonsel_conn_chain_delay4s_%d_%d_CtrlGroup.mat',bin,bin+1));
    
    %% Number of pairs
    PairsNum_hit_II_Nonmem(bin) = nnz(all(FRModulPairChains_go>0,2));
    PairsNum_hit_DD_Nonmem(bin) = nnz(all(FRModulPairChains_go<0,2));
    PairsNum_hit_ID_Nonmem(bin) = nnz(FRModulPairChains_go(:,1).*FRModulPairChains_go(:,2)==-1);
    
    %% Number of FC pairs in hit trials
    FcPairsNum_hit_II_Nonmem(bin) = nnz(all(FRModul_conn_chain_go>0,2));
    FcPairsNum_hit_DD_Nonmem(bin) = nnz(all(FRModul_conn_chain_go<0,2));
    FcPairsNum_hit_ID_Nonmem(bin) = nnz(FRModul_conn_chain_go(:,1).*FRModul_conn_chain_go(:,2)==-1);
end

%% FC density of pairs of memory neurons /// II, DD, and ID ///
for bin = 1:4 % four 1-second bin of delay
    load(sprintf('FCofMemoryNeurons_PropertyPopulation_%d_%d_%s.mat',bin,bin+1,Group));
    % pair
    pref_pairs{bin} = Pref_pair;
    % FC pair
    pref_chains_go{bin} = pref_chain_go;
    clear Pref_pair pref_chain_go
end
for bin = 1:4 % four 1-second bin of delay
    load(sprintf('PairsFRModulation_selec_conn_chain_delay4s_%d_%d_CtrlGroup.mat',bin,bin+1));
    
    %% Number of pairs
    % congruent active pairs
    PairsNum_II_ConAct(1,bin) = nnz(pref_pairs{bin}(:,2+bin)==pref_pairs{bin}(:,9+bin) & pref_pairs{bin}(:,2+bin)>0 & all(FRModulPairChains_go>0,2));
    PairsNum_DD_ConAct(1,bin) = nnz(pref_pairs{bin}(:,2+bin)==pref_pairs{bin}(:,9+bin) & pref_pairs{bin}(:,2+bin)>0 & all(FRModulPairChains_go<0,2));
    PairsNum_ID_ConAct(1,bin) = nnz(pref_pairs{bin}(:,2+bin)==pref_pairs{bin}(:,9+bin) & pref_pairs{bin}(:,2+bin)>0 & FRModulPairChains_go(:,1).*FRModulPairChains_go(:,2)==-1);
    % congruent inactive pairs
    PairsNum_II_ConInact(1,bin) = nnz(max(pref_pairs{bin}(:,3:6),[],2) == max(pref_pairs{bin}(:,10:13),[],2) & any(pref_pairs{bin}(:,[bin+2 bin+9])<1,2) & all(FRModulPairChains_go>0,2));
    PairsNum_DD_ConInact(1,bin) = nnz(max(pref_pairs{bin}(:,3:6),[],2) == max(pref_pairs{bin}(:,10:13),[],2) & any(pref_pairs{bin}(:,[bin+2 bin+9])<1,2) & all(FRModulPairChains_go<0,2));
    PairsNum_ID_ConInact(1,bin) = nnz(max(pref_pairs{bin}(:,3:6),[],2) == max(pref_pairs{bin}(:,10:13),[],2) & any(pref_pairs{bin}(:,[bin+2 bin+9])<1,2) & FRModulPairChains_go(:,1).*FRModulPairChains_go(:,2)==-1);
    % incongruent active pairs
    PairsNum_II_InconAct(1,bin) = nnz(pref_pairs{bin}(:,2+bin)>0 & pref_pairs{bin}(:,9+bin)>0 & pref_pairs{bin}(:,2+bin)~=pref_pairs{bin}(:,9+bin) & all(FRModulPairChains_go>0,2));
    PairsNum_DD_InconAct(1,bin) = nnz(pref_pairs{bin}(:,2+bin)>0 & pref_pairs{bin}(:,9+bin)>0 & pref_pairs{bin}(:,2+bin)~=pref_pairs{bin}(:,9+bin) & all(FRModulPairChains_go<0,2));
    PairsNum_ID_InconAct(1,bin) = nnz(pref_pairs{bin}(:,2+bin)>0 & pref_pairs{bin}(:,9+bin)>0 & pref_pairs{bin}(:,2+bin)~=pref_pairs{bin}(:,9+bin) & FRModulPairChains_go(:,1).*FRModulPairChains_go(:,2)==-1);
    % incongruent inactive pairs
    PairsNum_II_InconInact(1,bin) = nnz(max(pref_pairs{bin}(:,3:6),[],2)~=max(pref_pairs{bin}(:,10:13),[],2) & max(pref_pairs{bin}(:,3:6),[],2)>0 & max(pref_pairs{bin}(:,10:13),[],2)>0 & any(pref_pairs{bin}(:,[bin+2 bin+9])<1,2) & all(FRModulPairChains_go>0,2));
    PairsNum_DD_InconInact(1,bin) = nnz(max(pref_pairs{bin}(:,3:6),[],2)~=max(pref_pairs{bin}(:,10:13),[],2) & max(pref_pairs{bin}(:,3:6),[],2)>0 & max(pref_pairs{bin}(:,10:13),[],2)>0 & any(pref_pairs{bin}(:,[bin+2 bin+9])<1,2) & all(FRModulPairChains_go<0,2));
    PairsNum_ID_InconInact(1,bin) = nnz(max(pref_pairs{bin}(:,3:6),[],2)~=max(pref_pairs{bin}(:,10:13),[],2) & max(pref_pairs{bin}(:,3:6),[],2)>0 & max(pref_pairs{bin}(:,10:13),[],2)>0 & any(pref_pairs{bin}(:,[bin+2 bin+9])<1,2) & FRModulPairChains_go(:,1).*FRModulPairChains_go(:,2)==-1);
    
    %% Number of FC pairs in hit trials
    % congruent active pairs
    FcPairsNum_hit_II_ConAct(1,bin) = nnz(pref_chains_go{bin}(:,2+bin)==pref_chains_go{bin}(:,9+bin) & pref_chains_go{bin}(:,2+bin)>0 & all(FRModul_conn_chain_go>0,2));
    FcPairsNum_hit_DD_ConAct(1,bin) = nnz(pref_chains_go{bin}(:,2+bin)==pref_chains_go{bin}(:,9+bin) & pref_chains_go{bin}(:,2+bin)>0 & all(FRModul_conn_chain_go<0,2));
    FcPairsNum_hit_ID_ConAct(1,bin) = nnz(pref_chains_go{bin}(:,2+bin)==pref_chains_go{bin}(:,9+bin) & pref_chains_go{bin}(:,2+bin)>0 & FRModul_conn_chain_go(:,1).*FRModul_conn_chain_go(:,2)==-1);
    % congruent inactive pairs
    FcPairsNum_hit_II_ConInact(1,bin) = nnz(max(pref_chains_go{bin}(:,3:6),[],2)==max(pref_chains_go{bin}(:,10:13),[],2) & any(pref_chains_go{bin}(:,[bin+2 bin+9])<1,2) & all(FRModul_conn_chain_go>0,2));
    FcPairsNum_hit_DD_ConInact(1,bin) = nnz(max(pref_chains_go{bin}(:,3:6),[],2)==max(pref_chains_go{bin}(:,10:13),[],2) & any(pref_chains_go{bin}(:,[bin+2 bin+9])<1,2) & all(FRModul_conn_chain_go<0,2));
    FcPairsNum_hit_ID_ConInact(1,bin) = nnz(max(pref_chains_go{bin}(:,3:6),[],2)==max(pref_chains_go{bin}(:,10:13),[],2) & any(pref_chains_go{bin}(:,[bin+2 bin+9])<1,2) & FRModul_conn_chain_go(:,1).*FRModul_conn_chain_go(:,2)==-1);
    % incongruent active pairs
    FcPairsNum_hit_II_InconAct(1,bin) = nnz(pref_chains_go{bin}(:,2+bin)>0 & pref_chains_go{bin}(:,9+bin)>0 & pref_chains_go{bin}(:,2+bin)~=pref_chains_go{bin}(:,9+bin) & all(FRModul_conn_chain_go>0,2));
    FcPairsNum_hit_DD_InconAct(1,bin) = nnz(pref_chains_go{bin}(:,2+bin)>0 & pref_chains_go{bin}(:,9+bin)>0 & pref_chains_go{bin}(:,2+bin)~=pref_chains_go{bin}(:,9+bin) & all(FRModul_conn_chain_go<0,2));
    FcPairsNum_hit_ID_InconAct(1,bin) = nnz(pref_chains_go{bin}(:,2+bin)>0 & pref_chains_go{bin}(:,9+bin)>0 & pref_chains_go{bin}(:,2+bin)~=pref_chains_go{bin}(:,9+bin) & FRModul_conn_chain_go(:,1).*FRModul_conn_chain_go(:,2)==-1);
    % incongruent inactive pairs
    FcPairsNum_hit_II_InconInact(1,bin) = nnz(max(pref_chains_go{bin}(:,3:6),[],2)~=max(pref_chains_go{bin}(:,10:13),[],2) & max(pref_chains_go{bin}(:,3:6),[],2)>0 & max(pref_chains_go{bin}(:,10:13),[],2)>0 & any(pref_chains_go{bin}(:,[bin+2 bin+9])<1,2) & all(FRModul_conn_chain_go>0,2));
    FcPairsNum_hit_DD_InconInact(1,bin) = nnz(max(pref_chains_go{bin}(:,3:6),[],2)~=max(pref_chains_go{bin}(:,10:13),[],2) & max(pref_chains_go{bin}(:,3:6),[],2)>0 & max(pref_chains_go{bin}(:,10:13),[],2)>0 & any(pref_chains_go{bin}(:,[bin+2 bin+9])<1,2) & all(FRModul_conn_chain_go<0,2));
    FcPairsNum_hit_ID_InconInact(1,bin) = nnz(max(pref_chains_go{bin}(:,3:6),[],2)~=max(pref_chains_go{bin}(:,10:13),[],2) & max(pref_chains_go{bin}(:,3:6),[],2)>0 & max(pref_chains_go{bin}(:,10:13),[],2)>0 & any(pref_chains_go{bin}(:,[bin+2 bin+9])<1,2) & FRModul_conn_chain_go(:,1).*FRModul_conn_chain_go(:,2)==-1);
end

%% Compare FC density of congruent active, congruent inactive, incongruent ative, incongruent inactive, and non-memory pairs /// II, DD, and ID ///
PairsNumThres = 5;
% non-memory pairs
FcDensity_hit_DD_Nonmem = (FcPairsNum_hit_DD_Nonmem./PairsNum_hit_DD_Nonmem)';
FcDensity_hit_DD_Nonmem = FcDensity_hit_DD_Nonmem(PairsNum_hit_DD_Nonmem>=PairsNumThres);
FcDensity_hit_ID_Nonmem = (FcPairsNum_hit_ID_Nonmem./PairsNum_hit_ID_Nonmem)';
FcDensity_hit_ID_Nonmem = FcDensity_hit_ID_Nonmem(PairsNum_hit_ID_Nonmem>=PairsNumThres);
FcDensity_hit_II_Nonmem = (FcPairsNum_hit_II_Nonmem./PairsNum_hit_II_Nonmem)';
FcDensity_hit_II_Nonmem = FcDensity_hit_II_Nonmem(PairsNum_hit_II_Nonmem>=PairsNumThres);
% incongruent inactive pairs
FcDensity_hit_DD_InconInact = (FcPairsNum_hit_DD_InconInact./PairsNum_DD_InconInact)';
FcDensity_hit_DD_InconInact = FcDensity_hit_DD_InconInact(PairsNum_DD_InconInact>=PairsNumThres);
FcDensity_hit_ID_InconInact = (FcPairsNum_hit_ID_InconInact./PairsNum_ID_InconInact)';
FcDensity_hit_ID_InconInact = FcDensity_hit_ID_InconInact(PairsNum_ID_InconInact>=PairsNumThres);
FcDensity_hit_II_InconInact = (FcPairsNum_hit_II_InconInact./PairsNum_II_InconInact)';
FcDensity_hit_II_InconInact = FcDensity_hit_II_InconInact(PairsNum_II_InconInact>=PairsNumThres);
% incongruent active pairs
FcDensity_hit_DD_InconAct = (FcPairsNum_hit_DD_InconAct./PairsNum_DD_InconAct)';
FcDensity_hit_DD_InconAct = FcDensity_hit_DD_InconAct(PairsNum_DD_InconAct>=PairsNumThres);
FcDensity_hit_ID_InconAct = (FcPairsNum_hit_ID_InconAct./PairsNum_ID_InconAct)';
FcDensity_hit_ID_InconAct = FcDensity_hit_ID_InconAct(PairsNum_ID_InconAct>=PairsNumThres);
FcDensity_hit_II_InconAct = (FcPairsNum_hit_II_InconAct./PairsNum_II_InconAct)';
FcDensity_hit_II_InconAct = FcDensity_hit_II_InconAct(PairsNum_II_InconAct>=PairsNumThres);
% congruent inactive pairs
FcDensity_hit_DD_ConInact = (FcPairsNum_hit_DD_ConInact./PairsNum_DD_ConInact)';
FcDensity_hit_DD_ConInact = FcDensity_hit_DD_ConInact(PairsNum_DD_ConInact>=PairsNumThres);
FcDensity_hit_ID_ConInact = (FcPairsNum_hit_ID_ConInact./PairsNum_ID_ConInact)';
FcDensity_hit_ID_ConInact = FcDensity_hit_ID_ConInact(PairsNum_ID_ConInact>=PairsNumThres);
FcDensity_hit_II_ConInact = (FcPairsNum_hit_II_ConInact./PairsNum_II_ConInact)';
FcDensity_hit_II_ConInact = FcDensity_hit_II_ConInact(PairsNum_II_ConInact>=PairsNumThres);
% congruent active pairs
FcDensity_hit_DD_ConAct = (FcPairsNum_hit_DD_ConAct./PairsNum_DD_ConAct)';
FcDensity_hit_DD_ConAct = FcDensity_hit_DD_ConAct(PairsNum_DD_ConAct>=PairsNumThres);
FcDensity_hit_ID_ConAct = (FcPairsNum_hit_ID_ConAct./PairsNum_ID_ConAct)';
FcDensity_hit_ID_ConAct = FcDensity_hit_ID_ConAct(PairsNum_ID_ConAct>=PairsNumThres);
FcDensity_hit_II_ConAct = (FcPairsNum_hit_II_ConAct./PairsNum_II_ConAct)';
FcDensity_hit_II_ConAct = FcDensity_hit_II_ConAct(PairsNum_II_ConAct>=PairsNumThres);
% Two-way ANOVA
[p,tb1,stats] = unbalanced_anova_test(horzcat({FcDensity_hit_DD_Nonmem},{FcDensity_hit_ID_Nonmem},{FcDensity_hit_II_Nonmem}),...
    horzcat({FcDensity_hit_DD_InconInact},{FcDensity_hit_ID_InconInact},{FcDensity_hit_II_InconInact}),...
    horzcat({FcDensity_hit_DD_InconAct},{FcDensity_hit_ID_InconAct},{FcDensity_hit_II_InconAct}),...
    horzcat({FcDensity_hit_DD_ConInact},{FcDensity_hit_ID_ConInact},{FcDensity_hit_II_ConInact}),...
    horzcat({FcDensity_hit_DD_ConAct},{FcDensity_hit_ID_ConAct},{FcDensity_hit_II_ConAct}));

%% Figure
close all
C = [{[0 0 1]} {[0 0 0]} {[1 0 0]}];
figure('OuterPosition',[219 303 420 534]);
% non-memory pairs
if ~isempty(FcDensity_hit_DD_Nonmem)
    bar(1.1,mean(FcDensity_hit_DD_Nonmem),0.6,'facecolor',C{1},'edgecolor','none'); hold on
    if numel(FcDensity_hit_DD_Nonmem) > 1
        errorbar(1.1,mean(FcDensity_hit_DD_Nonmem),std(FcDensity_hit_DD_Nonmem)/sqrt(numel(FcDensity_hit_DD_Nonmem)),'color',C{1},'linestyle','none','marker','none'); hold on
    end
end
if ~isempty(FcDensity_hit_ID_Nonmem)
    bar(1.7,mean(FcDensity_hit_ID_Nonmem),0.6,'facecolor',C{2},'edgecolor','none'); hold on
    if numel(FcDensity_hit_ID_Nonmem) > 1
        errorbar(1.7,mean(FcDensity_hit_ID_Nonmem),std(FcDensity_hit_ID_Nonmem)/sqrt(numel(FcDensity_hit_ID_Nonmem)),'color',C{2},'linestyle','none','marker','none'); hold on
    end
end
if ~isempty(FcDensity_hit_II_Nonmem)
    bar(2.3,mean(FcDensity_hit_II_Nonmem),0.6,'facecolor',C{3},'edgecolor','none'); hold on
    if numel(FcDensity_hit_II_Nonmem) > 1
        errorbar(2.3,mean(FcDensity_hit_II_Nonmem),std(FcDensity_hit_II_Nonmem)/sqrt(numel(FcDensity_hit_II_Nonmem)),'color',C{3},'linestyle','none','marker','none'); hold on
    end
end
% incongruent inactive pairs
if ~isempty(FcDensity_hit_DD_InconInact)
    bar(3.1,mean(FcDensity_hit_DD_InconInact),0.6,'facecolor',C{1},'edgecolor','none'); hold on
    if numel(FcDensity_hit_DD_InconInact) > 1
        errorbar(3.1,mean(FcDensity_hit_DD_InconInact),std(FcDensity_hit_DD_InconInact)/sqrt(numel(FcDensity_hit_DD_InconInact)),'color',C{1},'linestyle','none','marker','none'); hold on
    end
end
if ~isempty(FcDensity_hit_ID_InconInact)
    bar(3.7,mean(FcDensity_hit_ID_InconInact),0.6,'facecolor',C{2},'edgecolor','none'); hold on
    if numel(FcDensity_hit_ID_InconInact) > 1
        errorbar(3.7,mean(FcDensity_hit_ID_InconInact),std(FcDensity_hit_ID_InconInact)/sqrt(numel(FcDensity_hit_ID_InconInact)),'color',C{2},'linestyle','none','marker','none'); hold on
    end
end
if ~isempty(FcDensity_hit_II_InconInact)
    bar(4.3,mean(FcDensity_hit_II_InconInact),0.6,'facecolor',C{3},'edgecolor','none'); hold on
    if numel(FcDensity_hit_II_InconInact) > 1
        errorbar(4.3,mean(FcDensity_hit_II_InconInact),std(FcDensity_hit_II_InconInact)/sqrt(numel(FcDensity_hit_II_InconInact)),'color',C{3},'linestyle','none','marker','none'); hold on
    end
end
% incongruent active pairs
if ~isempty(FcDensity_hit_DD_InconAct)
    bar(5.1,mean(FcDensity_hit_DD_InconAct),0.6,'facecolor',C{1},'edgecolor','none'); hold on
    if numel(FcDensity_hit_DD_InconAct) > 1
        errorbar(5.1,mean(FcDensity_hit_DD_InconAct),std(FcDensity_hit_DD_InconAct)/sqrt(numel(FcDensity_hit_DD_InconAct)),'color',C{1},'linestyle','none','marker','none'); hold on
    end
end
if ~isempty(FcDensity_hit_ID_InconAct)
    bar(5.7,mean(FcDensity_hit_ID_InconAct),0.6,'facecolor',C{2},'edgecolor','none'); hold on
    if numel(FcDensity_hit_ID_InconAct) > 1
        errorbar(5.7,mean(FcDensity_hit_ID_InconAct),std(FcDensity_hit_ID_InconAct)/sqrt(numel(FcDensity_hit_ID_InconAct)),'color',C{2},'linestyle','none','marker','none'); hold on
    end
end
if ~isempty(FcDensity_hit_II_InconAct)
    bar(6.3,mean(FcDensity_hit_II_InconAct),0.6,'facecolor',C{3},'edgecolor','none'); hold on
    if numel(FcDensity_hit_II_InconAct) > 1
        errorbar(6.3,mean(FcDensity_hit_II_InconAct),std(FcDensity_hit_II_InconAct)/sqrt(numel(FcDensity_hit_II_InconAct)),'color',C{3},'linestyle','none','marker','none'); hold on
    end
end
% congruent inactive pairs
if ~isempty(FcDensity_hit_DD_ConInact)
    bar(7.1,mean(FcDensity_hit_DD_ConInact),0.6,'facecolor',C{1},'edgecolor','none'); hold on
    if numel(FcDensity_hit_DD_ConInact) > 1
        errorbar(7.1,mean(FcDensity_hit_DD_ConInact),std(FcDensity_hit_DD_ConInact)/sqrt(numel(FcDensity_hit_DD_ConInact)),'color',C{1},'linestyle','none','marker','none'); hold on
    end
end
if ~isempty(FcDensity_hit_ID_ConInact)
    bar(7.7,mean(FcDensity_hit_ID_ConInact),0.6,'facecolor',C{2},'edgecolor','none'); hold on
    if numel(FcDensity_hit_ID_ConInact) > 1
        errorbar(7.7,mean(FcDensity_hit_ID_ConInact),std(FcDensity_hit_ID_ConInact)/sqrt(numel(FcDensity_hit_ID_ConInact)),'color',C{2},'linestyle','none','marker','none'); hold on
    end
end
if ~isempty(FcDensity_hit_II_ConInact)
    bar(8.3,mean(FcDensity_hit_II_ConInact),0.6,'facecolor',C{3},'edgecolor','none'); hold on
    if numel(FcDensity_hit_II_ConInact) > 1
        errorbar(8.3,mean(FcDensity_hit_II_ConInact),std(FcDensity_hit_II_ConInact)/sqrt(numel(FcDensity_hit_II_ConInact)),'color',C{3},'linestyle','none','marker','none'); hold on
    end
end
% congruent active pairs
if ~isempty(FcDensity_hit_DD_ConAct)
    bar(9.1,mean(FcDensity_hit_DD_ConAct),0.6,'facecolor',C{1},'edgecolor','none'); hold on
    if numel(FcDensity_hit_DD_ConAct) > 1
        errorbar(9.1,mean(FcDensity_hit_DD_ConAct),std(FcDensity_hit_DD_ConAct)/sqrt(numel(FcDensity_hit_DD_ConAct)),'color',C{1},'linestyle','none','marker','none'); hold on
    end
end
if ~isempty(FcDensity_hit_ID_ConAct)
    bar(9.7,mean(FcDensity_hit_ID_ConAct),0.6,'facecolor',C{2},'edgecolor','none'); hold on
    if numel(FcDensity_hit_ID_ConAct) > 1
        errorbar(9.7,mean(FcDensity_hit_ID_ConAct),std(FcDensity_hit_ID_ConAct)/sqrt(numel(FcDensity_hit_ID_ConAct)),'color',C{2},'linestyle','none','marker','none'); hold on
    end
end
if ~isempty(FcDensity_hit_II_ConAct)
    bar(10.3,mean(FcDensity_hit_II_ConAct),0.6,'facecolor',C{3},'edgecolor','none'); hold on
    if numel(FcDensity_hit_II_ConAct) > 1
        errorbar(10.3,mean(FcDensity_hit_II_ConAct),std(FcDensity_hit_II_ConAct)/sqrt(numel(FcDensity_hit_II_ConAct)),'color',C{3},'linestyle','none','marker','none'); hold on
    end
end
set(gca,'XTick',zeros(1,0),'xlim',[0.6 10.8]);
set(gca,'YTick',0:0.1:0.3,'YTickLabel',{'0','10','20','30'},'FontName','Arial','FontSize',16,'ylim',[0 0.3]);
ylabel('FC density (%)')
box off;
set(gcf,'Renderer','Painter'); saveas(gcf,'FcDensity_IncreaseDecreaseFrModulationOfAllTypesPairs_Local and Cross-region_hit','fig'); close;

