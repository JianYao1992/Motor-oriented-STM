%% FC density of neuronal pairs controlled for FR

clear; clc; close all;

%% Assignment
Group = 'CtrlGroup';
region = 'Local and Cross-region';
WhichTrials = 'hit trials';
if strcmp(region,'Within mPFC')
    TarReg1 = 1; TarReg2 = 1;
elseif strcmp(region,'Within aAIC')
    TarReg1 = 2; TarReg2 = 2;
elseif strcmp(region,'mPFC-aAIC')
    TarReg1 = 2; TarReg2 = 1;
elseif strcmp(region,'Local and Cross-region')
    TarReg1 = [1 2]; TarReg2 = [1 2];
end
PossPairNumThes = 5;
DelayBinNum = 4;
TargetFR = [4.5 9.5 14.5];
TriNumRange_hit = [-1 1000];
TriNumRange_CR = [60 80];

%% FC density /// controlled for FR ///
if strcmp(WhichTrials,'hit trials')
    fprintf('Calculate FC density when %d <= trial num <= %d\n',TriNumRange_hit(1),TriNumRange_hit(2));
else
    fprintf('Calculate FC density when %d <= trial num <= %d\n',TriNumRange_CR(1),TriNumRange_CR(2));
end
FcDensity_Nonmem = cell(1,numel(TargetFR));
FcDensity_ConAct = cell(1,numel(TargetFR));
FcDensity_ConInact = cell(1,numel(TargetFR));
FcDensity_InconAct = cell(1,numel(TargetFR));
FcDensity_InconInact = cell(1,numel(TargetFR));
for i = 1:length(TargetFR)
    FRRange = [TargetFR(i)-(TargetFR(2)-TargetFR(1))/2 TargetFR(i)+(TargetFR(2)-TargetFR(1))/2];
    FcDensity_Nonmem{i} = zeros(2,DelayBinNum);
    %% FC density of pairs of non-memory neurons
    for bin = 1:DelayBinNum   % four 1-second bins of delay
        load(sprintf('FCofNonmemoryNeurons_PropertyPopulation_1msbin_%d_%d_%s.mat',bin,bin+1,Group));
        load(sprintf('PairsFRModulation_1msbin_nonsel_conn_chain_delay4s_%d_%d_%s.mat',bin,bin+1,Group));
        if strcmp(WhichTrials,'hit trials')
            FcDensity_Nonmem{i}(1,bin) = nnz(ismember(reg_chain_go(:,1),TarReg1) & ismember(reg_chain_go(:,2),TarReg2) & trialnum_chain_go>=TriNumRange_hit(1) & trialnum_chain_go<=TriNumRange_hit(2) & FR_conn_Chains_go(:,1)>=FRRange(1) & FR_conn_Chains_go(:,1)<FRRange(2) & FR_conn_Chains_go(:,2)>=FRRange(1) & FR_conn_Chains_go(:,2)<FRRange(2));
            FcDensity_Nonmem{i}(2,bin) = nnz(ismember(Pair_reg(:,1),TarReg1) & ismember(Pair_reg(:,2),TarReg2) & Pair_trialnum(:,1)>=TriNumRange_hit(1) & Pair_trialnum(:,1)<=TriNumRange_hit(2) & FRPairChains_go(:,1)>=FRRange(1) & FRPairChains_go(:,1)<FRRange(2) & FRPairChains_go(:,2)>=FRRange(1) & FRPairChains_go(:,2)<FRRange(2));
        elseif strcmp(WhichTrials,'CR trials')
            FcDensity_Nonmem{i}(1,bin) = nnz(ismember(reg_chain_nogo(:,1),TarReg1) & ismember(reg_chain_nogo(:,2),TarReg2) & trialnum_chain_nogo>=TriNumRange_CR(1) & trialnum_chain_nogo<=TriNumRange_CR(2) & FR_conn_Chains_nogo(:,1)>=FRRange(1) & FR_conn_Chains_nogo(:,1)<FRRange(2) & FR_conn_Chains_nogo(:,2)>=FRRange(1) & FR_conn_Chains_nogo(:,2)<FRRange(2));
            FcDensity_Nonmem{i}(2,bin) = nnz(ismember(Pair_reg(:,1),TarReg1) & ismember(Pair_reg(:,2),TarReg2) & Pair_trialnum(:,2)>=TriNumRange_CR(1) & Pair_trialnum(:,2)<=TriNumRange_CR(2) & FRPairChains_nogo(:,1)>=FRRange(1) & FRPairChains_nogo(:,1)<FRRange(2) & FRPairChains_nogo(:,2)>=FRRange(1) & FRPairChains_nogo(:,2)<FRRange(2));
        end
        clear Pair_reg reg_chain_go reg_chain_nogo FRPairChains_go FRPairChains_nogo FR_conn_Chains_go FR_conn_Chains_nogo
    end
    
    %% FC density of pairs of memory neurons
    FcDensity_ConAct{i} = zeros(2,DelayBinNum);
    FcDensity_ConInact{i} = zeros(2,DelayBinNum);
    FcDensity_InconAct{i} = zeros(2,DelayBinNum);
    FcDensity_InconInact{i} = zeros(2,DelayBinNum);
    for bin = 1:DelayBinNum % four 1-second bins of delay
        load(sprintf('FCofMemoryNeurons_PropertyPopulation_1msbin_%d_%d_%s.mat',bin,bin+1,Group));
        load(sprintf('PairsFRModulation_1msbin_selec_conn_chain_delay4s_%d_%d_%s.mat',bin,bin+1,Group));
        % congruent active pairs
        if strcmp(WhichTrials,'hit trials')
            FcDensity_ConAct{i}(1,bin) = nnz(ismember(reg_chain_go(:,1),TarReg1) & ismember(reg_chain_go(:,2),TarReg2) & trialnum_chain_go>=TriNumRange_hit(1) & trialnum_chain_go<=TriNumRange_hit(2) & pref_chain_go(:,2+bin)==pref_chain_go(:,9+bin) & pref_chain_go(:,2+bin)>0 & FR_conn_Chains_go(:,1)>=FRRange(1) & FR_conn_Chains_go(:,1)<FRRange(2) & FR_conn_Chains_go(:,2)>=FRRange(1) & FR_conn_Chains_go(:,2)<FRRange(2));
            FcDensity_ConAct{i}(2,bin) = nnz(ismember(Pair_reg(:,1),TarReg1) & ismember(Pair_reg(:,2),TarReg2) & Pair_trialnum(:,1)>=TriNumRange_hit(1) & Pair_trialnum(:,1)<=TriNumRange_hit(2) & Pref_pair(:,2+bin)==Pref_pair(:,9+bin) & Pref_pair(:,2+bin)>0 & FRPairChains_go(:,1)>=FRRange(1) & FRPairChains_go(:,1)<FRRange(2) & FRPairChains_go(:,2)>=FRRange(1) & FRPairChains_go(:,2)<FRRange(2));
        else
            FcDensity_ConAct{i}(1,bin) = nnz(ismember(reg_chain_nogo(:,1),TarReg1) & ismember(reg_chain_nogo(:,2),TarReg2) & trialnum_chain_nogo>=TriNumRange_CR(1) & trialnum_chain_nogo<=TriNumRange_CR(2) & pref_chain_nogo(:,2+bin)==pref_chain_nogo(:,9+bin) & pref_chain_nogo(:,2+bin)>0 & FR_conn_Chains_nogo(:,1)>=FRRange(1) & FR_conn_Chains_nogo(:,1)<FRRange(2) & FR_conn_Chains_nogo(:,2)>=FRRange(1) & FR_conn_Chains_nogo(:,2)<FRRange(2));
            FcDensity_ConAct{i}(2,bin) = nnz(ismember(Pair_reg(:,1),TarReg1) & ismember(Pair_reg(:,2),TarReg2) & Pair_trialnum(:,2)>=TriNumRange_CR(1) & Pair_trialnum(:,2)<=TriNumRange_CR(2) & Pref_pair(:,2+bin)==Pref_pair(:,9+bin) & Pref_pair(:,2+bin)>0 & FRPairChains_nogo(:,1)>=FRRange(1) & FRPairChains_nogo(:,1)<FRRange(2) & FRPairChains_nogo(:,2)>=FRRange(1) & FRPairChains_nogo(:,2)<FRRange(2));
        end
        % congruent inactive pairs
        if strcmp(WhichTrials,'hit trials')
            FcDensity_ConInact{i}(1,bin) = nnz(ismember(reg_chain_go(:,1),TarReg1) & ismember(reg_chain_go(:,2),TarReg2) & trialnum_chain_go >= TriNumRange_hit(1) & trialnum_chain_go <= TriNumRange_hit(2) & max(pref_chain_go(:,3:6),[],2) == max(pref_chain_go(:,10:13),[],2) & any(pref_chain_go(:,[bin+2 bin+9])<1,2) & FR_conn_Chains_go(:,1) >= FRRange(1) & FR_conn_Chains_go(:,1) < FRRange(2) & FR_conn_Chains_go(:,2) >= FRRange(1) & FR_conn_Chains_go(:,2) < FRRange(2));
            FcDensity_ConInact{i}(2,bin) = nnz(ismember(Pair_reg(:,1),TarReg1) & ismember(Pair_reg(:,2),TarReg2) & Pair_trialnum(:,1) >= TriNumRange_hit(1) & Pair_trialnum(:,1) <= TriNumRange_hit(2) & max(Pref_pair(:,3:6),[],2) == max(Pref_pair(:,10:13),[],2) & any(Pref_pair(:,[bin+2 bin+9])<1,2) & FRPairChains_go(:,1) >= FRRange(1) & FRPairChains_go(:,1) < FRRange(2) & FRPairChains_go(:,2) >= FRRange(1) & FRPairChains_go(:,2) < FRRange(2));
        else
            FcDensity_ConInact{i}(1,bin) = nnz(ismember(reg_chain_nogo(:,1),TarReg1) & ismember(reg_chain_nogo(:,2),TarReg2) & trialnum_chain_nogo >= TriNumRange_CR(1) & trialnum_chain_nogo <= TriNumRange_CR(2) & max(pref_chain_nogo(:,3:6),[],2) == max(pref_chain_nogo(:,10:13),[],2) & any(pref_chain_nogo(:,[bin+2 bin+9])<1,2) & FR_conn_Chains_nogo(:,1) >= FRRange(1) & FR_conn_Chains_nogo(:,1) < FRRange(2) & FR_conn_Chains_nogo(:,2) >= FRRange(1) & FR_conn_Chains_nogo(:,2) < FRRange(2));
            FcDensity_ConInact{i}(2,bin) = nnz(ismember(Pair_reg(:,1),TarReg1) & ismember(Pair_reg(:,2),TarReg2) & Pair_trialnum(:,2) >= TriNumRange_CR(1) & Pair_trialnum(:,2) <= TriNumRange_CR(2) & max(Pref_pair(:,3:6),[],2) == max(Pref_pair(:,10:13),[],2) & any(Pref_pair(:,[bin+2 bin+9])<1,2) & FRPairChains_nogo(:,1) >= FRRange(1) & FRPairChains_nogo(:,1) < FRRange(2) & FRPairChains_nogo(:,2) >= FRRange(1) & FRPairChains_nogo(:,2) < FRRange(2));
        end
        % incongruent active pairs
        if strcmp(WhichTrials,'hit trials')
            FcDensity_InconAct{i}(1,bin) = nnz(ismember(reg_chain_go(:,1),TarReg1) & ismember(reg_chain_go(:,2),TarReg2) & trialnum_chain_go>=TriNumRange_hit(1) & trialnum_chain_go<=TriNumRange_hit(2) & pref_chain_go(:,2+bin)>0 & pref_chain_go(:,9+bin)>0 & pref_chain_go(:,2+bin)~=pref_chain_go(:,9+bin) & FR_conn_Chains_go(:,1)>=FRRange(1) & FR_conn_Chains_go(:,1)<FRRange(2) & FR_conn_Chains_go(:,2)>=FRRange(1) & FR_conn_Chains_go(:,2)<FRRange(2));
            FcDensity_InconAct{i}(2,bin) = nnz(ismember(Pair_reg(:,1),TarReg1) & ismember(Pair_reg(:,2),TarReg2) & Pair_trialnum(:,1)>=TriNumRange_hit(1) & Pair_trialnum(:,1)<=TriNumRange_hit(2) & Pref_pair(:,2+bin)>0 & Pref_pair(:,9+bin)>0 & Pref_pair(:,2+bin)~=Pref_pair(:,9+bin) & FRPairChains_go(:,1)>=FRRange(1) & FRPairChains_go(:,1)<FRRange(2) & FRPairChains_go(:,2)>=FRRange(1) & FRPairChains_go(:,2)<FRRange(2));
        else
            FcDensity_InconAct{i}(1,bin) = nnz(ismember(reg_chain_nogo(:,1),TarReg1) & ismember(reg_chain_nogo(:,2),TarReg2) & trialnum_chain_nogo>=TriNumRange_CR(1) & trialnum_chain_nogo<=TriNumRange_CR(2) & pref_chain_nogo(:,2+bin)>0 & pref_chain_nogo(:,9+bin)>0 & pref_chain_nogo(:,2+bin)~=pref_chain_nogo(:,9+bin) & FR_conn_Chains_nogo(:,1)>=FRRange(1) & FR_conn_Chains_nogo(:,1)<FRRange(2) & FR_conn_Chains_nogo(:,2)>=FRRange(1) & FR_conn_Chains_nogo(:,2)<FRRange(2));
            FcDensity_InconAct{i}(2,bin) = nnz(ismember(Pair_reg(:,1),TarReg1) & ismember(Pair_reg(:,2),TarReg2) & Pair_trialnum(:,2)>=TriNumRange_CR(1) & Pair_trialnum(:,2)<=TriNumRange_CR(2) & Pref_pair(:,2+bin)>0 & Pref_pair(:,9+bin)>0 & Pref_pair(:,2+bin)~=Pref_pair(:,9+bin) & FRPairChains_nogo(:,1)>=FRRange(1) & FRPairChains_nogo(:,1)<FRRange(2) & FRPairChains_nogo(:,2)>=FRRange(1) & FRPairChains_nogo(:,2)<FRRange(2));
        end
        % incongruent inactive pairs
        if strcmp(WhichTrials,'hit trials')
            FcDensity_InconInact{i}(1,bin) = nnz(ismember(reg_chain_go(:,1),TarReg1) & ismember(reg_chain_go(:,2),TarReg2) & trialnum_chain_go>=TriNumRange_hit(1) & trialnum_chain_go<=TriNumRange_hit(2) & max(pref_chain_go(:,3:6),[],2)~=max(pref_chain_go(:,10:13),[],2) & max(pref_chain_go(:,3:6),[],2)>0 & max(pref_chain_go(:,10:13),[],2)>0 & any(pref_chain_go(:,[bin+2 bin+9])<1,2) & FR_conn_Chains_go(:,1)>=FRRange(1) & FR_conn_Chains_go(:,1)<FRRange(2) & FR_conn_Chains_go(:,2)>=FRRange(1) & FR_conn_Chains_go(:,2)<FRRange(2));
            FcDensity_InconInact{i}(2,bin) = nnz(ismember(Pair_reg(:,1),TarReg1) & ismember(Pair_reg(:,2),TarReg2) & Pair_trialnum(:,1)>=TriNumRange_hit(1) & Pair_trialnum(:,1)<=TriNumRange_hit(2) & max(Pref_pair(:,3:6),[],2)~=max(Pref_pair(:,10:13),[],2) & max(Pref_pair(:,3:6),[],2)>0 & max(Pref_pair(:,10:13),[],2)>0 & any(Pref_pair(:,[bin+2 bin+9])<1,2) & FRPairChains_go(:,1) >= FRRange(1) & FRPairChains_go(:,1) < FRRange(2) & FRPairChains_go(:,2) >= FRRange(1) & FRPairChains_go(:,2) < FRRange(2));
        else
            FcDensity_InconInact{i}(1,bin) = nnz(ismember(reg_chain_nogo(:,1),TarReg1) & ismember(reg_chain_nogo(:,2),TarReg2) & trialnum_chain_nogo>=TriNumRange_CR(1) & trialnum_chain_nogo<=TriNumRange_CR(2) & max(pref_chain_nogo(:,3:6),[],2)~=max(pref_chain_nogo(:,10:13),[],2) & max(pref_chain_nogo(:,3:6),[],2)>0 & max(pref_chain_nogo(:,10:13),[],2)>0 & any(pref_chain_nogo(:,[bin+2 bin+9])<1,2) & FR_conn_Chains_nogo(:,1)>=FRRange(1) & FR_conn_Chains_nogo(:,1)<FRRange(2) & FR_conn_Chains_nogo(:,2)>=FRRange(1) & FR_conn_Chains_nogo(:,2)<FRRange(2));
            FcDensity_InconInact{i}(2,bin) = nnz(ismember(Pair_reg(:,1),TarReg1) & ismember(Pair_reg(:,2),TarReg2) & Pair_trialnum(:,2)>=TriNumRange_CR(1) & Pair_trialnum(:,2)<=TriNumRange_CR(2) & max(Pref_pair(:,3:6),[],2)~=max(Pref_pair(:,10:13),[],2) & max(Pref_pair(:,3:6),[],2)>0 & max(Pref_pair(:,10:13),[],2)>0 & any(Pref_pair(:,[bin+2 bin+9])<1,2) & FRPairChains_nogo(:,1)>=FRRange(1) & FRPairChains_nogo(:,1)<FRRange(2) & FRPairChains_nogo(:,2)>=FRRange(1) & FRPairChains_nogo(:,2)<FRRange(2));
        end
    end
end

%% Delete data without enough pairs number and assign activity level and type of neuronal pairs
FcDensity = vertcat(FcDensity_Nonmem,FcDensity_InconInact,FcDensity_InconAct,FcDensity_ConInact,FcDensity_ConAct);
data = []; 
PairsGroup = []; 
FR_ID = [];
for iRow = 1:size(FcDensity,1) % type of neuronal pairs
    for iColumn = 1:size(FcDensity,2) % activity level
        FcDensity{iRow,iColumn}(:,FcDensity{iRow,iColumn}(2,:)<PossPairNumThes) = [];
        PairsGroup = [PairsGroup; iRow*ones(size(FcDensity{iRow,iColumn},2),1)];
        FR_ID = [FR_ID; iColumn*ones(size(FcDensity{iRow,iColumn},2),1)];
        data = [data; FcDensity{iRow,iColumn}(1,:)'./FcDensity{iRow,iColumn}(2,:)'];
    end
end

%% Two-way ANOVA
[p,tbl] = anovan(data,{PairsGroup FR_ID},'model',2,'varnames',{'PairsGroup','FR_ID'},'display','on');

%% Error bar plot
C = [{[0 0 0]} {[0 0 1]} {[1 0 0]}];
figure('Color','w','Position',[100,100,800,500]);
for i = 1:numel(TargetFR) % FR
    % non-memory pairs
    tempdata = FcDensity{1,i};
    tempdata = tempdata(:,tempdata(2,:)>=PossPairNumThes);
    if ~isempty(tempdata)
        tempdata = tempdata(1,:)./tempdata(2,:);
        AverFcDensity_Nonmem = mean(tempdata);
        if numel(tempdata)>1
            errorbar(1,AverFcDensity_Nonmem,std(tempdata)/sqrt(numel(tempdata)),'color',C{i},'marker','o','markerfacecolor',C{i},'markeredgecolor','none','markersize',5,'capsize',6);
            hold on
        end
    else
        AverFcDensity_Nonmem = 0;
    end
    % incongruent inactive pairs
    tempdata = FcDensity{2,i};
    tempdata = tempdata(:,tempdata(2,:)>=PossPairNumThes);
    if ~isempty(tempdata)
        tempdata = tempdata(1,:)./tempdata(2,:);
        AverFcDensity_InconInact = mean(tempdata);
        if numel(tempdata)>1
            errorbar(2,AverFcDensity_InconInact,std(tempdata)/sqrt(numel(tempdata)),'color',C{i},'marker','o','markerfacecolor',C{i},'markeredgecolor','none','markersize',5,'capsize',6);
            hold on
        end
    else
        AverFcDensity_InconInact = 0;
    end
    % incongruent active pairs
    tempdata = FcDensity{3,i};
    tempdata = tempdata(:,tempdata(2,:)>=PossPairNumThes);
    if ~isempty(tempdata)
        tempdata = tempdata(1,:)./tempdata(2,:);
        AverFcDensity_InconAct = mean(tempdata);
        if numel(tempdata)>1
            errorbar(3,AverFcDensity_InconAct,std(tempdata)/sqrt(numel(tempdata)),'color',C{i},'marker','o','markerfacecolor',C{i},'markeredgecolor','none','markersize',5,'capsize',6);
            hold on
        end
    else
        AverFcDensity_InconAct = 0;
    end
    % congruent inactive pairs
    tempdata = FcDensity{4,i};
    tempdata = tempdata(:,tempdata(2,:)>=PossPairNumThes);
    if ~isempty(tempdata)
        tempdata = tempdata(1,:)./tempdata(2,:);
        AverFcDensity_ConInact = mean(tempdata);
        if numel(tempdata)>1
            errorbar(4,AverFcDensity_ConInact,std(tempdata)/sqrt(numel(tempdata)),'color',C{i},'marker','o','markerfacecolor',C{i},'markeredgecolor','none','markersize',5,'capsize',6);
            hold on
        end
    else
        AverFcDensity_ConInact = 0;
    end
    % congruent active pairs
    tempdata = FcDensity{5,i};
    tempdata = tempdata(:,tempdata(2,:)>=PossPairNumThes);
    if ~isempty(tempdata)
        tempdata = tempdata(1,:)./tempdata(2,:);
        AverFcDensity_ConAct = mean(tempdata);
        if numel(tempdata)>1
            errorbar(5,AverFcDensity_ConAct,std(tempdata)/sqrt(numel(tempdata)),'color',C{i},'marker','o','markerfacecolor',C{i},'markeredgecolor','none','markersize',5,'capsize',6);
            hold on
        end
    else
        AverFcDensity_ConAct = 0;
    end
    plot(1:5,horzcat(AverFcDensity_Nonmem,AverFcDensity_InconInact,AverFcDensity_InconAct,AverFcDensity_ConInact,AverFcDensity_ConAct),'color',C{i},'marker','none');
    hold on
end
if strcmp(WhichTrials,'hit trials')
    TrialNumRange = TriNumRange_hit;
elseif strcmp(WhichTrials,'CR trials')
    TrialNumRange = TriNumRange_CR;
end
set(gcf,'Render','Painter');
saveas(gcf,sprintf('FC density-1msbin-FR %d to %d Hz-TrialNumRange %d to %d-%s-%s-%s',TargetFR(1)-(TargetFR(2)-TargetFR(1))/2,TargetFR(end)+(TargetFR(2)-TargetFR(1))/2,...
    TrialNumRange(1),TrialNumRange(2),WhichTrials,region,Group),'fig');
close all;
save(sprintf('FC density-1msbin-FR %d to %d Hz-TrialNumRange %d to %d-%s-%s-%s',TargetFR(1)-(TargetFR(2)-TargetFR(1))/2,TargetFR(end)+(TargetFR(2)-TargetFR(1))/2,...
    TrialNumRange(1),TrialNumRange(2),WhichTrials,region,Group),'FcDensity','data','FR_ID','PairsGroup','tbl','-v7.3');