%% Information of region, ID, and selectivity direction of neurons of pairs 
% switched neurons are removed
% neurons are divided into memory and non-memory neurons

clear; clc; close all;

%% Assignemnt
Group = 'ChR2Group'; 

%% Information of all pairs for each bin, excluding pairs with switched neurons
for bin = -1:1:5  % from baseline to response period
    disp(['Analyze ' num2str(bin) 'th bin data' ]);
    Pair_reg = []; Pair_chain = []; Pref_pair = []; Pair_mouse = []; Pair_trialnum = [];
    reg_chain_go = []; conn_chain_go = []; pref_chain_go = []; mouse_chain_go = []; trialnum_chain_go = [];
    reg_chain_nogo = []; conn_chain_nogo = []; pref_chain_nogo = []; mouse_chain_nogo = []; trialnum_chain_nogo = [];
    stats = [];
    load(sprintf('TestFC_XCORR_stats_%d_%d_1msbin_%s.mat',bin,bin+1,Group));
    for iPair = 1:length(stats)
        disp(['Analyze ' num2str(iPair) 'th Pair' ]);
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
        
        %% all pairs without switched neurons
        a = unique(stats{iPair}.prefered_sample_su1(3:6)); a(a==0) = []; % preference of unit 1
        b = unique(stats{iPair}.prefered_sample_su2(3:6)); b(b==0) = []; % preference of unit 2
        if isequal(max(a),min(a)) && isequal(max(b),min(b))
            Pair_reg = [Pair_reg; tempreg1,tempreg2];
            Pair_chain = [Pair_chain; stats{iPair}.su1_clusterid, stats{iPair}.su2_clusterid];
            Pref_pair = [Pref_pair; stats{iPair}.prefered_sample_su1, stats{iPair}.prefered_sample_su2];
            Pair_mouse = [Pair_mouse; stats{iPair}.fileidx];
            Pair_trialnum = [Pair_trialnum; horzcat(length(stats{iPair}.s1_trials),length(stats{iPair}.s2_trials))];
            if isfield(stats{iPair},'AIs1') && stats{iPair}.AIs1 ~= 0
                reg_chain_go = [reg_chain_go; tempreg1,tempreg2,stats{iPair}.AIs1]; % region of FC pair in context 1
                conn_chain_go = [conn_chain_go; stats{iPair}.su1_clusterid,stats{iPair}.su2_clusterid,stats{iPair}.AIs1]; % ID of neurons of FC pair in context 1
                pref_chain_go = [pref_chain_go; stats{iPair}.prefered_sample_su1,stats{iPair}.prefered_sample_su2]; % selectivity of neurons in FC pair in context 1
                mouse_chain_go = [mouse_chain_go; stats{iPair}.fileidx]; % ID of file of FC pair in context 1
                trialnum_chain_go = [trialnum_chain_go; length(stats{iPair}.s1_trials)];
            end
            if isfield(stats{iPair},'AIs2') && stats{iPair}.AIs2 ~= 0
                reg_chain_nogo = [reg_chain_nogo; tempreg1,tempreg2,stats{iPair}.AIs2]; % region of FC pair in context 2
                conn_chain_nogo = [conn_chain_nogo; stats{iPair}.su1_clusterid,stats{iPair}.su2_clusterid,stats{iPair}.AIs2]; % ID of neurons of FC pair in context 2
                pref_chain_nogo = [pref_chain_nogo; stats{iPair}.prefered_sample_su1,stats{iPair}.prefered_sample_su2]; % selectivity of neurons in FC pair in context 2
                mouse_chain_nogo = [mouse_chain_nogo; stats{iPair}.fileidx]; % ID of file of FC pair in context 2
                trialnum_chain_nogo = [trialnum_chain_nogo; length(stats{iPair}.s2_trials)];
            end
        end
    end
    Pair_reg1 = Pair_reg; Pair_chain1 = Pair_chain; Pref_pair1 = Pref_pair; Pair_mouse1 = Pair_mouse; Pair_trialnum1 = Pair_trialnum;
    reg_chain_go1 = reg_chain_go; conn_chain_go1 = conn_chain_go; pref_chain_go1 = pref_chain_go; mouse_chain_go1 = mouse_chain_go; trialnum_chain_go1 = trialnum_chain_go;
    reg_chain_nogo1 = reg_chain_nogo; conn_chain_nogo1 = conn_chain_nogo; pref_chain_nogo1 = pref_chain_nogo; mouse_chain_nogo1 = mouse_chain_nogo; trialnum_chain_nogo1 = trialnum_chain_nogo;
    
    %% Divide nonmemory and memory neuron
    % both neurons are memory neurons
    % //////possible pairs///////
    a = any(Pref_pair(:,3:6),2) & any(Pref_pair(:,10:13),2);
    Pair_reg = Pair_reg(a,:);
    Pair_chain = Pair_chain(a,:);
    Pref_pair = Pref_pair(a,:);
    Pair_mouse = Pair_mouse(a,:);
    Pair_trialnum = Pair_trialnum(a,:);
    % //////FC pairs in Go trials//////
    a = any(pref_chain_go(:,3:6),2) & any(pref_chain_go(:,10:13),2);
    reg_chain_go = reg_chain_go(a,:);
    conn_chain_go = conn_chain_go(a,:);
    pref_chain_go = pref_chain_go(a,:);
    mouse_chain_go = mouse_chain_go(a,:);
    trialnum_chain_go = trialnum_chain_go(a,:);
    % //////FC pairs in NoGo trials//////
    a = any(pref_chain_nogo(:,3:6),2) & any(pref_chain_nogo(:,10:13),2);
    reg_chain_nogo = reg_chain_nogo(a,:);
    conn_chain_nogo = conn_chain_nogo(a,:);
    pref_chain_nogo = pref_chain_nogo(a,:);
    mouse_chain_nogo = mouse_chain_nogo(a,:);
    trialnum_chain_nogo = trialnum_chain_nogo(a,:);
    save(sprintf('FCofMemoryNeurons_PropertyPopulation_1msbin_%d_%d_%s.mat',bin,bin+1,Group),'Pair_reg','Pair_chain','Pref_pair','Pair_mouse','Pair_trialnum',...
        'reg_chain_go','conn_chain_go','pref_chain_go','mouse_chain_go','trialnum_chain_go','reg_chain_nogo','conn_chain_nogo','pref_chain_nogo',...
        'mouse_chain_nogo','trialnum_chain_nogo','-v7.3');
    % both neurons are non-memory neurons
    a = all(Pref_pair1(:,[3 4 5 6 10 11 12 13])<1,2); 
    Pair_reg = Pair_reg1(a,:);
    Pair_chain = Pair_chain1(a,:);
    Pref_pair = Pref_pair1(a,:);
    Pair_mouse = Pair_mouse1(a,:);
    Pair_trialnum = Pair_trialnum1(a,:);
    a = all(pref_chain_go1(:,[3 4 5 6 10 11 12 13])<1,2); % two non-memory neurons in FC pairs in context 1
    reg_chain_go = reg_chain_go1(a,:);
    conn_chain_go = conn_chain_go1(a,:);
    pref_chain_go = pref_chain_go1(a,:);
    mouse_chain_go = mouse_chain_go1(a,:);
    trialnum_chain_go = trialnum_chain_go1(a,:);
    a = all(pref_chain_nogo1(:,[3 4 5 6 10 11 12 13])<1,2); % two non-memory neurons in FC pairs in context 2
    reg_chain_nogo = reg_chain_nogo1(a,:);
    conn_chain_nogo = conn_chain_nogo1(a,:);
    pref_chain_nogo = pref_chain_nogo1(a,:);
    mouse_chain_nogo = mouse_chain_nogo1(a,:);
    trialnum_chain_nogo = trialnum_chain_nogo1(a,:);
    save(sprintf('FCofNonmemoryNeurons_PropertyPopulation_1msbin_%d_%d_%s.mat',bin,bin+1,Group),'Pair_reg','Pair_chain','Pref_pair','Pair_mouse','Pair_trialnum',...
        'reg_chain_go','conn_chain_go','pref_chain_go','mouse_chain_go','trialnum_chain_go','reg_chain_nogo','conn_chain_nogo','pref_chain_nogo','mouse_chain_nogo','trialnum_chain_nogo','-v7.3');
    clearvars Pair_reg1 Pair_chain1 Pref_pair1 reg_chain_go1 conn_chain_go1 pref_chain_go1 reg_chain_nogo1 conn_chain_nogo1 pref_chain_nogo1
end


