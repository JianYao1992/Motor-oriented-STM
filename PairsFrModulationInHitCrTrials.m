%% FR modulation of two neurons in neuronal pairs

clear; clc; close all;

%% Assignment
Group = 'ChR2Group';

%% Load result of FR modulation
load(['FRSiginCorrectTrialsRelativetoBaseline' Group]);

%% FR modulation of pairs
for bin = 1:4
    %% Pairs of non-memory neurons
    disp(['Analyze data in ' num2str(bin) 'th bin for pairs of non-memory neurons']);
    load(sprintf('FCofNonmemoryNeurons_PropertyPopulation_1msbin_%d_%d_%s.mat',bin,bin+1,Group));
    FRPairChains_go = zeros(size(Pair_chain,1),2);
    FRPairChains_nogo = zeros(size(Pair_chain,1),2);
    FR_conn_Chains_go = zeros(size(conn_chain_go,1),2);
    FR_conn_Chains_nogo = zeros(size(conn_chain_nogo,1),2);
    FRModulPairChains_go = zeros(size(Pair_chain,1),2);
    FRModulPairChains_nogo = zeros(size(Pair_chain,1),2);
    FRModul_conn_chain_go = zeros(size(conn_chain_go,1),2);
    FRModul_conn_chain_nogo = zeros(size(conn_chain_nogo,1),2);
    for iPair = 1:size(Pair_chain,1) % pairs
        Unit1 = Pair_chain(iPair,1);
        Unit2 = Pair_chain(iPair,2); 
        if ~isempty(FRSiginCorrectTrials{Unit1,1}) % Go trials
            FRPairChains_go(iPair,1) =  FRinCorrectTrials{Unit1,1}(bin);
            FRPairChains_go(iPair,2) = FRinCorrectTrials{Unit2,1}(bin); 
            FRModulPairChains_go(iPair,1) = FRSiginCorrectTrials{Unit1,1}(bin);
            FRModulPairChains_go(iPair,2) = FRSiginCorrectTrials{Unit2,1}(bin);
        else
            FRPairChains_go(iPair,1) = 0;
            FRPairChains_go(iPair,2) = 0;
            FRModulPairChains_go(iPair,1) = 0;
            FRModulPairChains_go(iPair,2) = 0;
        end
        if ~isempty(FRSiginCorrectTrials{Unit1,2}) % NoGo trials
            FRPairChains_nogo(iPair,1) =  FRinCorrectTrials{Unit1,2}(bin);
            FRPairChains_nogo(iPair,2) = FRinCorrectTrials{Unit2,2}(bin); 
            FRModulPairChains_nogo(iPair,1) = FRSiginCorrectTrials{Unit1,2}(bin); 
            FRModulPairChains_nogo(iPair,2) = FRSiginCorrectTrials{Unit2,2}(bin);
        else
            FRPairChains_nogo(iPair,1) = 0;
            FRPairChains_nogo(iPair,2) = 0;
            FRModulPairChains_nogo(iPair,1) = 0;
            FRModulPairChains_nogo(iPair,2) = 0;
        end
    end
    for iPair = 1:size(conn_chain_go,1) % FC pairs in Go trials
        Unit1 = conn_chain_go(iPair,1); 
        Unit2 = conn_chain_go(iPair,2);
        FR_conn_Chains_go(iPair,1) = FRinCorrectTrials{Unit1,1}(bin);
        FR_conn_Chains_go(iPair,2) = FRinCorrectTrials{Unit2,1}(bin); 
        FRModul_conn_chain_go(iPair,1) = FRSiginCorrectTrials{Unit1,1}(bin);
        FRModul_conn_chain_go(iPair,2) = FRSiginCorrectTrials{Unit2,1}(bin);
    end
    for iPair = 1:size(conn_chain_nogo,1) % FC pairs in NoGo trials
        Unit1 = conn_chain_nogo(iPair,1);
        Unit2 = conn_chain_nogo(iPair,2); 
        FR_conn_Chains_nogo(iPair,1) = FRinCorrectTrials{Unit1,2}(bin);
        FR_conn_Chains_nogo(iPair,2) = FRinCorrectTrials{Unit2,2}(bin); 
        FRModul_conn_chain_nogo(iPair,1) = FRSiginCorrectTrials{Unit1,2}(bin);
        FRModul_conn_chain_nogo(iPair,2) = FRSiginCorrectTrials{Unit2,2}(bin);
    end
    save(sprintf('PairsFRModulation_1msbin_nonsel_conn_chain_delay4s_%d_%d_%s.mat',bin,bin+1,Group),'FRPairChains_go','FRPairChains_nogo','FR_conn_Chains_go','FR_conn_Chains_nogo','FRModulPairChains_go','FRModulPairChains_nogo','FRModul_conn_chain_go','FRModul_conn_chain_nogo','-v7.3');
    clear Pair_chain conn_chain_go conn_chain_nogo
    
    %% Pairs of memory neurons
    disp(['Analyze data in ' num2str(bin) 'th bin for pairs of memory neurons']);
    load(sprintf('FCofMemoryNeurons_PropertyPopulation_1msbin_%d_%d_%s.mat',bin,bin+1,Group));
    FRPairChains_go = zeros(size(Pair_chain,1),2);
    FRPairChains_nogo = zeros(size(Pair_chain,1),2);
    FR_conn_Chains_go = zeros(size(conn_chain_go,1),2); 
    FR_conn_Chains_nogo = zeros(size(conn_chain_nogo,1),2);
    FRModulPairChains_go = zeros(size(Pair_chain,1),2);
    FRModulPairChains_nogo = zeros(size(Pair_chain,1),2);
    FRModul_conn_chain_go = zeros(size(conn_chain_go,1),2); 
    FRModul_conn_chain_nogo = zeros(size(conn_chain_nogo,1),2);
    for iPair = 1:size(Pair_chain,1) % pairs
        Unit1 = Pair_chain(iPair,1);
        Unit2 = Pair_chain(iPair,2); 
        if ~isempty(FRSiginCorrectTrials{Unit1,1}) % Go trials
            FRPairChains_go(iPair,1) = FRinCorrectTrials{Unit1,1}(bin);
            FRPairChains_go(iPair,2) = FRinCorrectTrials{Unit2,1}(bin); 
            FRModulPairChains_go(iPair,1) = FRSiginCorrectTrials{Unit1,1}(bin);
            FRModulPairChains_go(iPair,2) = FRSiginCorrectTrials{Unit2,1}(bin);
        else
            FRPairChains_go(iPair,1) = 0; 
            FRPairChains_go(iPair,2) = 0;
            FRModulPairChains_go(iPair,1) = 0;
            FRModulPairChains_go(iPair,2) = 0;
        end
        if ~isempty(FRSiginCorrectTrials{Unit1,2}) % NoGo trials
            FRPairChains_nogo(iPair,1) = FRinCorrectTrials{Unit1,2}(bin);
            FRPairChains_nogo(iPair,2) = FRinCorrectTrials{Unit2,2}(bin); 
            FRModulPairChains_nogo(iPair,1) = FRSiginCorrectTrials{Unit1,2}(bin);
            FRModulPairChains_nogo(iPair,2) = FRSiginCorrectTrials{Unit2,2}(bin);
        else
            FRPairChains_nogo(iPair,1) = 0;
            FRPairChains_nogo(iPair,2) = 0;
            FRModulPairChains_nogo(iPair,1) = 0;
            FRModulPairChains_nogo(iPair,2) = 0;
        end
    end
    for iPair = 1:size(conn_chain_go,1) % FC pairs in Go trials
        Unit1 = conn_chain_go(iPair,1);
        Unit2 = conn_chain_go(iPair,2);
        FR_conn_Chains_go(iPair,1) = FRinCorrectTrials{Unit1,1}(bin);
        FR_conn_Chains_go(iPair,2) = FRinCorrectTrials{Unit2,1}(bin); 
        FRModul_conn_chain_go(iPair,1) = FRSiginCorrectTrials{Unit1,1}(bin); 
        FRModul_conn_chain_go(iPair,2) = FRSiginCorrectTrials{Unit2,1}(bin);
    end
    for iPair = 1:size(conn_chain_nogo,1) % FC pairs in NoGo trials
        Unit1 = conn_chain_nogo(iPair,1);
        Unit2 = conn_chain_nogo(iPair,2); 
        FR_conn_Chains_nogo(iPair,1) = FRinCorrectTrials{Unit1,2}(bin);
        FR_conn_Chains_nogo(iPair,2) = FRinCorrectTrials{Unit2,2}(bin); 
        FRModul_conn_chain_nogo(iPair,1) = FRSiginCorrectTrials{Unit1,2}(bin); 
        FRModul_conn_chain_nogo(iPair,2) = FRSiginCorrectTrials{Unit2,2}(bin);
    end
    save(sprintf('PairsFRModulation_1msbin_selec_conn_chain_delay4s_%d_%d_%s.mat',bin,bin+1,Group),'FRPairChains_go','FRPairChains_nogo','FR_conn_Chains_go','FR_conn_Chains_nogo','FRModulPairChains_go','FRModulPairChains_nogo','FRModul_conn_chain_go','FRModul_conn_chain_nogo','-v7.3');
end
    
