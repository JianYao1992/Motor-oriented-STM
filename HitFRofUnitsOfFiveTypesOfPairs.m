%% FR of neurons of pairs

clear; clc; close all;

%% Assignment
Group = 'CtrlGroup';
region = 'Local and Cross-region';
TarReg1 = [1 2];
TarReg2 = [1 2];
DelayBinNum = 4;
WhichTrials = 'hit trials';
if strcmp(WhichTrials,'hit trials')
    TrialNumRange = [-1 1000];
elseif strcmp(WhichTrials,'CR trials')
    TrialNumRange = [60 80];
end

%% FR of neurons in non-memory pairs in the range of number of correct trials
FR_NonmemPair = cell(1,DelayBinNum);
for bin = 1:DelayBinNum % four 1-second bins of delay
    load(sprintf('FCofNonmemoryNeurons_PropertyPopulation_%d_%d_%s.mat',bin,bin+1,Group));
    load(sprintf('PairsFRModulation_nonsel_conn_chain_delay4s_%d_%d_%s.mat',bin,bin+1,Group));
    % pairs with enough number of trials
    if strcmp(WhichTrials,'hit trials')
        temp = ismember(Pair_reg(:,1),TarReg1) & ismember(Pair_reg(:,2),TarReg2) & Pair_trialnum(:,1)>=TrialNumRange(1) & Pair_trialnum(:,1)<=TrialNumRange(2);
        tempFR = FRPairChains_go(temp,:);
    elseif strcmp(WhichTrials,'CR trials')
        temp = Pair_chain(ismember(Pair_reg(:,1),TarReg1) & ismember(Pair_reg(:,2),TarReg2) & Pair_trialnum(:,2)>=TrialNumRange(1) & Pair_trialnum(:,2)<=TrialNumRange(2),:);
        tempFR = FRPairChains_nogo(temp,:);
    end
    Max_tempFR = max(tempFR,[],2);
    Min_tempFR = min(tempFR,[],2);
    FR_NonmemPair{bin} = horzcat(Min_tempFR,Max_tempFR);
    clear Pair_reg reg_chain_go reg_chain_nogo FRPairChains_go FRPairChains_nogo
end
%% FR of memory neurons in the range of number of correct trials
FR_CongActPair = cell(1,DelayBinNum);
FR_CongInactPair = cell(1,DelayBinNum);
FR_IncongActPair = cell(1,DelayBinNum);
FR_IncongInactPair = cell(1,DelayBinNum);
for bin = 1:DelayBinNum % four 1-second bins of delay
    load(sprintf('FCofMemoryNeurons_PropertyPopulation_%d_%d_%s.mat',bin,bin+1,Group));
    load(sprintf('PairsFRModulation_selec_conn_chain_delay4s_%d_%d_%s.mat',bin,bin+1,Group));
    % congruent active pair
    if strcmp(WhichTrials,'hit trials')
        temp = ismember(Pair_reg(:,1),TarReg1) & ismember(Pair_reg(:,2),TarReg2) & Pair_trialnum(:,1)>=TrialNumRange(1) & Pair_trialnum(:,1)<=TrialNumRange(2) & Pref_pair(:,2+bin)==Pref_pair(:,9+bin) & Pref_pair(:,2+bin)>0;
        tempFR = FRPairChains_go(temp,:);
    elseif strcmp(WhichTrials,'CR trials')
        temp = ismember(Pair_reg(:,1),TarReg1) & ismember(Pair_reg(:,2),TarReg2) & Pair_trialnum(:,2)>=TrialNumRange(1) & Pair_trialnum(:,2)<=TrialNumRange(2) & Pref_pair(:,2+bin)==Pref_pair(:,9+bin) & Pref_pair(:,2+bin) > 0;
        tempFR = FRPairChains_nogo(temp,:);
    end
    Max_tempFR = max(tempFR,[],2);
    Min_tempFR = min(tempFR,[],2);
    FR_CongActPair{bin} = horzcat(Min_tempFR,Max_tempFR);
    % congruent inactive pair
    if strcmp(WhichTrials,'hit trials')
        temp = ismember(Pair_reg(:,1),TarReg1) & ismember(Pair_reg(:,2),TarReg2) & Pair_trialnum(:,1)>=TrialNumRange(1) & Pair_trialnum(:,1)<=TrialNumRange(2) & max(Pref_pair(:,3:6),[],2)==max(Pref_pair(:,10:13),[],2) & any(Pref_pair(:,[bin+2 bin+9])<1,2);
        tempFR = FRPairChains_go(temp,:);
    elseif strcmp(WhichTrials,'CR trials')
        temp = ismember(Pair_reg(:,1),TarReg1) & ismember(Pair_reg(:,2),TarReg2) & Pair_trialnum(:,2)>=TrialNumRange(1) & Pair_trialnum(:,2)<=TrialNumRange(2) & max(Pref_pair(:,3:6),[],2)==max(Pref_pair(:,10:13),[],2) & any(Pref_pair(:,[bin+2 bin+9])<1,2);
        tempFR = FRPairChains_nogo(temp,:);
    end
    Max_tempFR = max(tempFR,[],2);
    Min_tempFR = min(tempFR,[],2);
    FR_CongInactPair{bin} = horzcat(Min_tempFR,Max_tempFR);
    % incongruent active pair
    if strcmp(WhichTrials,'hit trials')
        temp = ismember(Pair_reg(:,1),TarReg1) & ismember(Pair_reg(:,2),TarReg2) & Pair_trialnum(:,1)>=TrialNumRange(1) & Pair_trialnum(:,1)<=TrialNumRange(2) & Pref_pair(:,2+bin)>0 & Pref_pair(:,9+bin)>0 & Pref_pair(:,2+bin)~=Pref_pair(:,9+bin);
        tempFR = FRPairChains_go(temp,:);
    elseif strcmp(WhichTrials,'CR trials')
        temp = ismember(Pair_reg(:,1),TarReg1) & ismember(Pair_reg(:,2),TarReg2) & Pair_trialnum(:,2)>=TrialNumRange(1) & Pair_trialnum(:,2)<=TrialNumRange(2) & Pref_pair(:,2+bin)>0 & Pref_pair(:,9+bin)>0 & Pref_pair(:,2+bin)~=Pref_pair(:,9+bin);
        tempFR = FRPairChains_nogo(temp,:);
    end
    Max_tempFR = max(tempFR,[],2);
    Min_tempFR = min(tempFR,[],2);
    FR_IncongActPair{bin} = horzcat(Min_tempFR,Max_tempFR);
    % incongruent inactive pair
    if strcmp(WhichTrials,'hit trials')
        temp = ismember(Pair_reg(:,1),TarReg1) & ismember(Pair_reg(:,2),TarReg2) & Pair_trialnum(:,1)>=TrialNumRange(1) & Pair_trialnum(:,1)<=TrialNumRange(2) & max(Pref_pair(:,3:6),[],2)~=max(Pref_pair(:,10:13),[],2) & max(Pref_pair(:,3:6),[],2)>0 & max(Pref_pair(:,10:13),[],2)>0 & any(Pref_pair(:,[bin+2 bin+9])<1,2);
        tempFR = FRPairChains_go(temp,:);
    elseif strcmp(WhichTrials,'CR trials')
        temp = ismember(Pair_reg(:,1),TarReg1) & ismember(Pair_reg(:,2),TarReg2) & Pair_trialnum(:,2)>=TrialNumRange(1) & Pair_trialnum(:,2)<=TrialNumRange(2) & max(Pref_pair(:,3:6),[],2)~=max(Pref_pair(:,10:13),[],2) & max(Pref_pair(:,3:6),[],2)>0 & max(Pref_pair(:,10:13),[],2)>0 & any(Pref_pair(:,[bin+2 bin+9])<1,2);
        tempFR = FRPairChains_nogo(temp,:);
    end
    Max_tempFR = max(tempFR,[],2);
    Min_tempFR = min(tempFR,[],2);
    FR_IncongInactPair{bin} = horzcat(Min_tempFR,Max_tempFR);
end

%% Merge result
FR = vertcat(FR_NonmemPair,FR_IncongInactPair,FR_IncongActPair,FR_CongInactPair,FR_CongActPair);

%% Higher and lower FR of neurons of five types of pairs
data = []; HL = []; NP = []; ID = [];
for iRow = 1:size(FR,1) % pair type
    for iColumn = 1:size(FR,2) % higher vs. lower FR of neurons in a pair
        temp = mean(FR{iRow,iColumn});
        data = [data; temp(:)];
        NP = [NP; iRow*ones(size(FR{iRow,iColumn},2),1)];
        HL = [HL; [1;2]];
        ID = [ID; (iColumn+(iRow-1)*size(FR,2))*ones(size(FR{iRow,iColumn},2),1)];
    end
end

%% Tw-ANOVA-md
Frame = [data NP HL ID];
[SSQs, DFs, MSQs, Fs, Ps] = mixed_between_within_anova(Frame,0,sprintf('FR between neuronal pairs-%s-%s-%s',region,Group,WhichTrials));

%% Error bar plot
C = [{[0 0 1]} {[1 0 0]}];
AverFR = cellfun(@mean,FR,'UniformOutput',0);
NewAverFR = cell(1,size(AverFR,1));
for i = 1:size(AverFR,1) 
    tempNewData = AverFR(i,:);
    NewAverFR{i} = vertcat(tempNewData{:});
end
figure('Color','w','Position',[100,100,800,500]);
for i = 1:2
    tempNewAveFR = cellfun(@(x) x(:,i),NewAverFR,'UniformOutput',0);
    meanvalue = cellfun(@mean,tempNewAveFR);
    semvalue = cellfun(@(x) std(x)/sqrt(size(x,1)),tempNewAveFR);
    errorbar(1:numel(meanvalue),meanvalue,semvalue,'color',C{i},'marker','o','markerfacecolor',C{i},'markeredgecolor','none','markersize',5,'capsize',6);
    hold on
end
box off
set(gca,'XLim',[0.5 size(FR,1)+0.5],'YLim',[0 12]);
title(sprintf('Pnp = %d, Phl = %d',Ps{1},Ps{3}));
set(gcf,'Render','Painter'); saveas(gcf,sprintf('FR between neuronal pairs-%s-%s-%s',region,Group,WhichTrials),'fig'); close all;



