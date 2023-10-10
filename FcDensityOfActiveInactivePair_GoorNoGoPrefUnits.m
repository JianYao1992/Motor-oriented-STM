%% Compare FC density within Go-preferred and within NoGo-preferred neurons of active or inactive pairs

clear; clc; close all;

%% Assignment
DelayBinsNum = 4;
PairNumThres = 10;

%% Session-based FC density in each bin of delay
FCDensity_GoPrefUnits_ActPair = cell(1,DelayBinsNum);
FCDensity_GoPrefUnits_InactPair = cell(1,DelayBinsNum);
FCDensity_NoGoPrefUnits_ActPair = cell(1,DelayBinsNum);
FCDensity_NoGoPrefUnits_InactPair = cell(1,DelayBinsNum);
for bin = 1:DelayBinsNum
    load(sprintf('FCofNeurons_PropertyPopulation_%d_%d_CtrlGroup',bin,bin+1));
    % session ID of neuronal pairs of Go-preferred neurons
    tempSess_GoPrefUnits_ActPair = Pair_mouse(Pref_pair(:,2+bin)>0 & Pref_pair(:,9+bin)>0 & max(Pref_pair(:,3:6),[],2)==1 & max(Pref_pair(:,10:13),[],2)==1);
    tempSess_GoPrefUnits_ActFC = mouse_chain_go(pref_chain_go(:,2+bin)>0 & pref_chain_go(:,9+bin)>0 & max(pref_chain_go(:,3:6),[],2)==1 & max(pref_chain_go(:,10:13),[],2)==1);
    tempSess_GoPrefUnits_InactPair = Pair_mouse(max(Pref_pair(:,3:6),[],2)==1 & max(Pref_pair(:,10:13),[],2)==1 & Pref_pair(:,2+bin).*Pref_pair(:,9+bin)==0);
    tempSess_GoPrefUnits_InactFC = mouse_chain_go(max(pref_chain_go(:,3:6),[],2)==1 & max(pref_chain_go(:,10:13),[],2)==1 & pref_chain_go(:,2+bin).*pref_chain_go(:,9+bin)==0);
    % session ID of neuronal pairs of NoGo-preferred neurons
    tempSess_NoGoPrefUnits_ActPair = Pair_mouse(Pref_pair(:,2+bin)==2 & Pref_pair(:,9+bin)==2 & ~(any(Pref_pair(:,3:6)==1,2)|any(Pref_pair(:,10:13)==1,2)));
    tempSess_NoGoPrefUnits_ActFC = mouse_chain_go(pref_chain_go(:,2+bin)==2 & pref_chain_go(:,9+bin)==2 & ~(any(pref_chain_go(:,3:6)==1,2)|any(pref_chain_go(:,10:13)==1,2)));
    tempSess_NoGoPrefUnits_InactPair = Pair_mouse(max(Pref_pair(:,3:6),[],2)==2 & max(Pref_pair(:,10:13),[],2)==2 & Pref_pair(:,2+bin).*Pref_pair(:,9+bin)==0 & ~(any(Pref_pair(:,3:6)==1,2)|any(Pref_pair(:,10:13)==1,2)));
    tempSess_NoGoPrefUnits_InactFC = mouse_chain_go(max(pref_chain_go(:,3:6),[],2)==2 & max(pref_chain_go(:,10:13),[],2)==2 & pref_chain_go(:,2+bin).*pref_chain_go(:,9+bin)==0 & ~(any(pref_chain_go(:,3:6)==1,2)|any(pref_chain_go(:,10:13)==1,2)));
    % target session ID above threshold //////active pairs of Go-preferred neurons//////
    [C,~,ic] = unique(tempSess_GoPrefUnits_ActPair); 
    tempPairNum = accumarray(ic,1);
    C = C(tempPairNum>=PairNumThres);
    tempPairNum = tempPairNum(tempPairNum>=PairNumThres);
    for iSess = 1:numel(C)
        FCDensity_GoPrefUnits_ActPair{bin} = [FCDensity_GoPrefUnits_ActPair{bin}; nnz(tempSess_GoPrefUnits_ActFC==C(iSess))/tempPairNum(iSess)];
    end
    % target session ID above threshold //////inactive pairs of Go-preferred neurons//////
    [C,~,ic] = unique(tempSess_GoPrefUnits_InactPair); 
    tempPairNum = accumarray(ic,1);
    C = C(tempPairNum>=PairNumThres);
    tempPairNum = tempPairNum(tempPairNum>=PairNumThres);
    for iSess = 1:numel(C)
        FCDensity_GoPrefUnits_InactPair{bin} = [FCDensity_GoPrefUnits_InactPair{bin}; nnz(tempSess_GoPrefUnits_InactFC==C(iSess))/tempPairNum(iSess)];
    end
    % target session ID above threshold //////active pairs of NoGo-preferred neurons//////
    [C,~,ic] = unique(tempSess_NoGoPrefUnits_ActPair); 
    tempPairNum = accumarray(ic,1);
    C = C(tempPairNum>=PairNumThres);
    tempPairNum = tempPairNum(tempPairNum>=PairNumThres);
    for iSess = 1:numel(C)
        FCDensity_NoGoPrefUnits_ActPair{bin} = [FCDensity_NoGoPrefUnits_ActPair{bin}; nnz(tempSess_NoGoPrefUnits_ActFC==C(iSess))/tempPairNum(iSess)];
    end
    % target session ID above threshold //////inactive pairs of NoGo-preferred neurons//////
    [C,~,ic] = unique(tempSess_NoGoPrefUnits_InactPair); 
    tempPairNum = accumarray(ic,1);
    C = C(tempPairNum>=PairNumThres);
    tempPairNum = tempPairNum(tempPairNum>=PairNumThres);
    for iSess = 1:numel(C)
        FCDensity_NoGoPrefUnits_InactPair{bin} = [FCDensity_NoGoPrefUnits_InactPair{bin}; nnz(tempSess_NoGoPrefUnits_InactFC==C(iSess))/tempPairNum(iSess)];
    end
    clear Pref_pair pref_chain_go Pair_mouse mouse_chain_go
end

%% Two-way ANOVA ///between active pairs of Go- and NoGo- preferred neurons///
FCDensity_ActPair = vertcat(FCDensity_GoPrefUnits_ActPair,FCDensity_NoGoPrefUnits_ActPair);
data = [];
PairType = [];
BinID = [];
for iType = 1:size(FCDensity_ActPair,1)
    for iBin = 1:size(FCDensity_ActPair,2)
        data = [data; FCDensity_ActPair{iType,iBin}];
        PairType = [PairType; iType*ones(size(FCDensity_ActPair{iType,iBin}))];
        BinID = [BinID; iBin*ones(size(FCDensity_ActPair{iType,iBin}))];
    end
end
[p_ActPair,tbl_ActPair] = anovan(data,{PairType BinID},'model',2,'varnames',{'Pair','Time'});

%% Figure ///error bar plot, between Go- and NoGo-preferred neurons in Active pairs///
figure('position',[200 200 300 500]);
AverDensity_GoPref_Act = cellfun(@mean,FCDensity_GoPrefUnits_ActPair);
Sem_GoPref_Act = cellfun(@(x) std(x)/sqrt(numel(x)),FCDensity_GoPrefUnits_ActPair);
errorbar(1:numel(AverDensity_GoPref_Act),AverDensity_GoPref_Act,Sem_GoPref_Act,'r','marker','o','markerfacecolor','r','markeredgecolor','none'); hold on
AverDensity_NoGoPref_Act = cellfun(@mean,FCDensity_NoGoPrefUnits_ActPair);
Sem_NoGoPref_Act = cellfun(@(x) std(x)/sqrt(numel(x)),FCDensity_NoGoPrefUnits_ActPair);
errorbar(1:numel(AverDensity_NoGoPref_Act),AverDensity_NoGoPref_Act,Sem_NoGoPref_Act,'r','marker','o','linestyle','--','markerfacecolor','r','markeredgecolor','none');
set(gca,'XTick',1:1:DelayBinsNum,'xlim',[0.8 DelayBinsNum+0.2],'FontName','Arial','FontSize',16,'ylim',[0 1.25*max(max(AverDensity_GoPref_Act),max(AverDensity_NoGoPref_Act))]);
ylabel('FC density (%)','FontSize',18,'FontName','Arial'); box off;
set(gcf,'Renderer','Painter'); saveas(gcf,'Compare FC density between Go- and NoGo-preferred neurons_Active','fig'); close all;
save('Compare FC density between Go- and NoGo-preferred neurons_Active','tbl_ActPair','-v7.3');

%% Two-way ANOVA ///between inactive pairs of Go- and NoGo- preferred neurons///
FCDensity_InactPair = vertcat(FCDensity_GoPrefUnits_InactPair,FCDensity_NoGoPrefUnits_InactPair);
data = [];
PairType = [];
BinID = [];
for iType = 1:size(FCDensity_InactPair,1)
    for iBin = 1:size(FCDensity_InactPair,2)
        data = [data; FCDensity_InactPair{iType,iBin}];
        PairType = [PairType; iType*ones(size(FCDensity_InactPair{iType,iBin}))];
        BinID = [BinID; iBin*ones(size(FCDensity_InactPair{iType,iBin}))];
    end
end
[p_InactPair,tbl_InactPair] = anovan(data,{PairType BinID},'model',2,'varnames',{'Pair','Time'});

%% Figure ///error bar plot, between Go- and NoGo-preferred neurons in Inactive pairs///
figure('position',[200 200 300 500]);
AverDensity_GoPref_Inact = cellfun(@mean,FCDensity_GoPrefUnits_InactPair);
Sem_GoPref_Inact = cellfun(@(x) std(x)/sqrt(numel(x)),FCDensity_GoPrefUnits_InactPair);
errorbar(1:numel(AverDensity_GoPref_Inact),AverDensity_GoPref_Inact,Sem_GoPref_Inact,'b','marker','o','markerfacecolor','b','markeredgecolor','none'); hold on
AverDensity_NoGoPref_Inact = cellfun(@mean,FCDensity_NoGoPrefUnits_InactPair);
Sem_NoGoPref_Inact = cellfun(@(x) std(x)/sqrt(numel(x)),FCDensity_NoGoPrefUnits_InactPair);
errorbar(1:numel(AverDensity_NoGoPref_Inact),AverDensity_NoGoPref_Inact,Sem_NoGoPref_Inact,'b','marker','o','linestyle','--','markerfacecolor','b','markeredgecolor','none');
set(gca,'XTick',1:1:DelayBinsNum,'xlim',[0.8 DelayBinsNum+0.2],'FontName','Arial','FontSize',16,'ylim',[0 1.25*max(max(AverDensity_GoPref_Inact),max(AverDensity_NoGoPref_Inact))]);
ylabel('FC density (%)','FontSize',18,'FontName','Arial'); box off;
set(gcf,'Renderer','Painter'); saveas(gcf,'Compare FC density between Go- and NoGo-preferred neurons_Inactive','fig'); close all;
save('Compare FC density between Go- and NoGo-preferred neurons_Inactive','tbl_InactPair','-v7.3');
