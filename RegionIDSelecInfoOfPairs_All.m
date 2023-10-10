%% Information of region, ID, and selectivity direction of neurons of pairs
% all neuronal pairs, no matter memory or non-memory neuron

clear; clc; close all;

%% Assignment
Group = 'ChR2Group';
BinSize = 1;

%% Information of all pairs for each bin
for bin = -1:1:5
    disp(['Analyze ' num2str(bin) 'th bin data' ]);
    Pair_reg = []; Pair_chain = []; Pref_pair = []; Pair_mouse = [];
    reg_chain_go = []; conn_chain_go = []; pref_chain_go = []; mouse_chain_go = [];
    reg_chain_nogo = []; conn_chain_nogo = []; pref_chain_nogo = []; mouse_chain_nogo = [];
    stats = [];
    load(['TestFC_XCORR_stats_' num2str(bin) '_' num2str(bin+BinSize) '_1msbin_GoNoGo_' Group '.mat']);
    for iPair = 1:length(stats)
        disp(['Analyzing ' num2str(iPair) 'th Pair' ]);
        if strcmp(stats{iPair}.reg_su1,'mPFC') % region of unit 1
            tempreg1 = 1;
        else
            tempreg1 = 2;
        end
        if strcmp(stats{iPair}.reg_su2,'mPFC') % region of unit 2
            tempreg2 = 1;
        else
            tempreg2 = 2;
        end
        Pair_reg = [Pair_reg; tempreg1,tempreg2];
        Pair_chain = [Pair_chain; stats{iPair}.su1_clusterid, stats{iPair}.su2_clusterid];
        Pref_pair = [Pref_pair; stats{iPair}.prefered_sample_su1, stats{iPair}.prefered_sample_su2];
        Pair_mouse = [Pair_mouse; stats{iPair}.fileidx];
        if isfield(stats{iPair},'AIs1') && stats{iPair}.AIs1 ~= 0
            reg_chain_go = [reg_chain_go; tempreg1,tempreg2,stats{iPair}.AIs1]; % region of FC pair in context 1
            conn_chain_go = [conn_chain_go; stats{iPair}.su1_clusterid,stats{iPair}.su2_clusterid,stats{iPair}.AIs1]; % ID of neurons of FC pair in context 1
            pref_chain_go = [pref_chain_go; stats{iPair}.prefered_sample_su1,stats{iPair}.prefered_sample_su2]; % selectivity of neurons in FC pair in context 1
            mouse_chain_go = [mouse_chain_go; stats{iPair}.fileidx]; % ID of file of FC pair in context 1
        end
        if isfield(stats{iPair},'AIs2') && stats{iPair}.AIs2 ~= 0 
            reg_chain_nogo = [reg_chain_nogo; tempreg1,tempreg2,stats{iPair}.AIs2]; % region of FC pair in context 2
            conn_chain_nogo = [conn_chain_nogo; stats{iPair}.su1_clusterid,stats{iPair}.su2_clusterid,stats{iPair}.AIs2]; % ID of neurons of FC pair in context 2
            pref_chain_nogo = [pref_chain_nogo; stats{iPair}.prefered_sample_su1,stats{iPair}.prefered_sample_su2]; % selectivity of neurons in FC pair in context 2
            mouse_chain_nogo = [mouse_chain_nogo; stats{iPair}.fileidx]; % ID of file of FC pair in context 2
        end
    end
    save(['FCofNeurons_PropertyPopulation_1msbin_' num2str(bin) '_' num2str(bin+BinSize) '_GoNoGo_' Group '.mat'],'Pair_reg','Pair_chain','Pref_pair','Pair_mouse',...
        'reg_chain_go','conn_chain_go','pref_chain_go','mouse_chain_go','reg_chain_nogo','conn_chain_nogo','pref_chain_nogo','mouse_chain_nogo','-v7.3');
end


