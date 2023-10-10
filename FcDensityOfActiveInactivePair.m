%% FC density of active and inactive pairs

clear; clc; close all;

%% Assignment
DelayBinsNum = 4;
PairNumThres = 60;

%% Session-based FC density in each bin of delay
FCDensity_Act = cell(1 ,DelayBinsNum);
FCDensity_Inact = cell(1,DelayBinsNum);
for bin = 1:DelayBinsNum
    load(sprintf('FCofNeurons_PropertyPopulation_1msbin_%d_%d_CtrlGroup',bin,bin+1));
    % session ID of neuronal pairs
    tempSess_ActPair = Pair_mouse(Pref_pair(:,2+bin)>0 & Pref_pair(:,9+bin)>0);
    tempSess_InactPair = Pair_mouse(max(Pref_pair(:,3:6),[],2)>0 & max(Pref_pair(:,10:13),[],2)>0 & Pref_pair(:,2+bin).*Pref_pair(:,9+bin)==0);
    tempSess_ActFC = mouse_chain_go(pref_chain_go(:,2+bin)>0 & pref_chain_go(:,9+bin)>0);
    tempSess_InactFC = mouse_chain_go(max(pref_chain_go(:,3:6),[],2)>0 & max(pref_chain_go(:,10:13),[],2)>0 & pref_chain_go(:,2+bin).*pref_chain_go(:,9+bin)==0);
    % target session ID above threshold
    [C,~,ic] = unique(tempSess_ActPair); % active pairs
    tempPairNum_Act = accumarray(ic,1);
    C = C(tempPairNum_Act>=PairNumThres);
    tempPairNum_Act = tempPairNum_Act(tempPairNum_Act>=PairNumThres);
    for iSess = 1:numel(C)
        FCDensity_Act{bin} = [FCDensity_Act{bin}; nnz(tempSess_ActFC==C(iSess))/tempPairNum_Act(iSess)];
    end
    [C,~,ic] = unique(tempSess_InactPair); % inactive pairs
    tempPairNum_Inact = accumarray(ic,1);
    C = C(tempPairNum_Inact>=PairNumThres);
    tempPairNum_Inact = tempPairNum_Inact(tempPairNum_Inact>=PairNumThres);
    for iSess = 1:numel(C)
        FCDensity_Inact{bin} = [FCDensity_Inact{bin}; nnz(tempSess_InactFC==C(iSess))/tempPairNum_Inact(iSess)];
    end
    clear Pref_pair pref_chain_go Pair_mouse mouse_chain_go
end

%% Two-way ANOVA
FCDensity = vertcat(FCDensity_Act,FCDensity_Inact);
data = [];
PairType = [];
BinID = [];
for iType = 1:size(FCDensity,1)
    for iBin = 1:size(FCDensity,2)
        data = [data; FCDensity{iType,iBin}];
        PairType = [PairType; iType*ones(size(FCDensity{iType,iBin}))];
        BinID = [BinID; iBin*ones(size(FCDensity{iType,iBin}))];
    end
end
[p,tbl] = anovan(data,{PairType BinID},'model',2,'varnames',{'Pair','Time'});

%% Error bar plot
figure('position',[200 200 300 500]);
AverDensity_Act = cellfun(@mean,FCDensity_Act);
Sem_Act = cellfun(@(x) std(x)/sqrt(numel(x)),FCDensity_Act);
errorbar(1:numel(AverDensity_Act),AverDensity_Act,Sem_Act,'r','marker','o','markerfacecolor','r','markeredgecolor','none'); hold on
AverDensity_Inact = cellfun(@mean,FCDensity_Inact);
Sem_Inact = cellfun(@(x) std(x)/sqrt(numel(x)),FCDensity_Inact);
errorbar(1:numel(AverDensity_Inact),AverDensity_Inact,Sem_Inact,'b','marker','o','markerfacecolor','b','markeredgecolor','none');
set(gca,'XTick',1:1:DelayBinsNum,'xlim',[0.8 DelayBinsNum+0.2],'FontName','Arial','FontSize',16,'ylim',[0 1.25*max(max(AverDensity_Act),max(AverDensity_Inact))]);
ylabel('FC density (%)','FontSize',18,'FontName','Arial'); box off;
set(gcf,'Renderer','Painter'); saveas(gcf,'Compare FC density between active and inactiveFC_1msbin','fig'); close all;
save('Compare FC density between active and inactiveFC1msbin','tbl','-v7.3');


