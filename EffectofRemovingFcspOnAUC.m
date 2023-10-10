%% Effect of removing FCSP events on STM-encoding ability of following neurons
%% //////ecreased STM-encoding ability following removing FCSP events//////

clear; clc; close all;

%% Assignment
Group = 'CtrlGroup';
Reg = 'Within aAIC';
IsSemOrCiPlot = 2; % 1: SEM; 2: 95% CI
DelayBinID = 1:1:4;
BaseLen = 1;
OdorLen = 1;
MaxNum_LU = 6;
NumThres_FCpair = 5;

%% Target FC neuronal pairs
FcPairs_Go_MemtoMem = []; % each cell stands for FC neuronal pairs with different numbers of leading neurons
FcPairs_NoGo_MemtoMem = [];
FcPairs_Go_NonmemtoMem = [];
FcPairs_NoGo_NonmemtoMem = [];
for i = 1:numel(DelayBinID)
    %% FC neuronal pairs in target region
    load(sprintf('FCofNeurons_PropertyPopulation_1msbin_%d_%d_%s.mat',DelayBinID(i),DelayBinID(i)+1,Group));
    if strcmp(Reg,'Within mPFC')
        temp_Go = ismember(reg_chain_go(:,1),1) & ismember(reg_chain_go(:,2),1);
        temp_NoGo = ismember(reg_chain_nogo(:,1),1) & ismember(reg_chain_nogo(:,2),1);
    elseif strcmp(Reg,'Within aAIC')
        temp_Go = ismember(reg_chain_go(:,1),2) & ismember(reg_chain_go(:,2),2);
        temp_NoGo = ismember(reg_chain_nogo(:,1),2) & ismember(reg_chain_nogo(:,2),2);
    elseif strcmp(Reg,'mPFC-aAIC')
        temp_Go = ismember(reg_chain_go(:,1),2) & ismember(reg_chain_go(:,2),1) & reg_chain_go(:,3)<0;
        temp_NoGo = ismember(reg_chain_nogo(:,1),2) & ismember(reg_chain_nogo(:,2),1) & reg_chain_nogo(:,3)<0;
    elseif strcmp(Reg,'aAIC-mPFC')
        temp_Go = ismember(reg_chain_go(:,1),2) & ismember(reg_chain_go(:,2),1) & reg_chain_go(:,3)>0;
        temp_NoGo = ismember(reg_chain_nogo(:,1),2) & ismember(reg_chain_nogo(:,2),1) & reg_chain_nogo(:,3)>0;
    end
    conn_chain_go = conn_chain_go(temp_Go,:);
    pref_chain_go = pref_chain_go(temp_Go,:);
    reg_chain_go = reg_chain_go(temp_Go,:);
    conn_chain_nogo = conn_chain_nogo(temp_NoGo,:);
    pref_chain_nogo = pref_chain_nogo(temp_NoGo,:);
    reg_chain_nogo = reg_chain_nogo(temp_NoGo,:);
    
    %% Leading to following direction
    for iPair = 1:size(conn_chain_go,1) % hit trials
        if conn_chain_go(iPair,3) < 0
            temp = conn_chain_go(iPair,1); conn_chain_go(iPair,1) = conn_chain_go(iPair,2); conn_chain_go(iPair,2) = temp;
            conn_chain_go(iPair,3) = -1*conn_chain_go(iPair,3);
            temp = pref_chain_go(iPair,1:7); pref_chain_go(iPair,1:7) = pref_chain_go(iPair,8:14); pref_chain_go(iPair,8:14) = temp;
            temp = reg_chain_go(iPair,1); reg_chain_go(iPair,1) = reg_chain_go(iPair,2); reg_chain_go(iPair,2) = temp;
            reg_chain_go(iPair,3) = -1*reg_chain_go(iPair,3);
        end
    end
    for iPair = 1:size(conn_chain_nogo,1) % CR trials
        if conn_chain_nogo(iPair,3) < 0
            temp = conn_chain_nogo(iPair,1); conn_chain_nogo(iPair,1) = conn_chain_nogo(iPair,2); conn_chain_nogo(iPair,2) = temp;
            conn_chain_nogo(iPair,3) = -1*conn_chain_nogo(iPair,3);
            temp = pref_chain_nogo(iPair,1:7); pref_chain_nogo(iPair,1:7) = pref_chain_nogo(iPair,8:14); pref_chain_nogo(iPair,8:14) = temp;
            temp = reg_chain_nogo(iPair,1); reg_chain_nogo(iPair,1) = reg_chain_nogo(iPair,2); reg_chain_nogo(iPair,2) = temp;
            reg_chain_nogo(iPair,3) = -1*reg_chain_nogo(iPair,3);
        end
    end
    
    %% Divide FC neuronal pairs into memory to memory and non-memory to memory pairs
    conn_Go_MemtoMem = conn_chain_go(pref_chain_go(:,BaseLen+OdorLen+DelayBinID(i))>0 & pref_chain_go(:,BaseLen+OdorLen+DelayBinID(i)+size(pref_chain_go,2)/2)>0,:);
    conn_NoGo_MemtoMem = conn_chain_nogo(pref_chain_nogo(:,BaseLen+OdorLen+DelayBinID(i))>0 & pref_chain_nogo(:,BaseLen+OdorLen+DelayBinID(i)+size(pref_chain_nogo,2)/2)>0,:);
    conn_Go_NonmemtoMem = conn_chain_go(pref_chain_go(:,BaseLen+OdorLen+DelayBinID(i))==0 & pref_chain_go(:,BaseLen+OdorLen+DelayBinID(i)+size(pref_chain_go,2)/2)>0,:);
    conn_NoGo_NonmemtoMem = conn_chain_nogo(pref_chain_nogo(:,BaseLen+OdorLen+DelayBinID(i))==0 & pref_chain_nogo(:,BaseLen+OdorLen+DelayBinID(i)+size(pref_chain_nogo,2)/2)>0,:);
    
    %% Divide FC neuronal pairs into groups with different numbers of leading neurons
    tempFCpair_Go_MemtoMem = GetPairAlignedWithUniqueFu(conn_Go_MemtoMem,DelayBinID(i),MaxNum_LU);
    tempFCpair_NoGo_MemtoMem = GetPairAlignedWithUniqueFu(conn_NoGo_MemtoMem,DelayBinID(i),MaxNum_LU);
    tempFCpair_Go_NonmemtoMem = GetPairAlignedWithUniqueFu(conn_Go_NonmemtoMem,DelayBinID(i),MaxNum_LU);
    tempFCpair_NoGo_NonmemtoMem = GetPairAlignedWithUniqueFu(conn_NoGo_NonmemtoMem,DelayBinID(i),MaxNum_LU);
    FcPairs_Go_MemtoMem = [FcPairs_Go_MemtoMem; tempFCpair_Go_MemtoMem];
    FcPairs_NoGo_MemtoMem = [FcPairs_NoGo_MemtoMem; tempFCpair_NoGo_MemtoMem];
    FcPairs_Go_NonmemtoMem = [FcPairs_Go_NonmemtoMem; tempFCpair_Go_NonmemtoMem];
    FcPairs_NoGo_NonmemtoMem = [FcPairs_NoGo_NonmemtoMem; tempFCpair_NoGo_NonmemtoMem];
end
FCpair_MemtoMem(1,:) = SumFCofDelayBins(FcPairs_Go_MemtoMem);
FCpair_MemtoMem(2,:) = SumFCofDelayBins(FcPairs_NoGo_MemtoMem);
FCpair_NonmemtoMem(1,:) = SumFCofDelayBins(FcPairs_Go_NonmemtoMem);
FCpair_NonmemtoMem(2,:) = SumFCofDelayBins(FcPairs_NoGo_NonmemtoMem);

%% ID of mPFC and aAIC neurons
load(sprintf('AllNeuronsAUCValue_%s.mat',Group));

%% Load trial-based spike rasters for mPFC and aAIC neurons
load(sprintf('NeuronOriginandSpikeRasterInformationfor%s.mat',Group));

%% FR of leading neurons
FR_MemtoMem = ComputeAllLuFR(FCpair_MemtoMem,Reg,ID_mPFC,mPFCSpikeTimeinIndividualTrial,mPFCUnitTrialMark,ID_aAIC,aAICSpikeTimeinIndividualTrial,aAICUnitTrialMark,BaseLen,OdorLen);
FR_NonmemtoMem = ComputeAllLuFR(FCpair_NonmemtoMem,Reg,ID_mPFC,mPFCSpikeTimeinIndividualTrial,mPFCUnitTrialMark,ID_aAIC,aAICSpikeTimeinIndividualTrial,aAICUnitTrialMark,BaseLen,OdorLen);
for i = 1:size(FR_MemtoMem,1)
    for j = 1:size(FR_MemtoMem,2)
        tempID = find(FR_MemtoMem{i,j}<60);
        FR_MemtoMem{i,j} = FR_MemtoMem{i,j}(tempID);
        FCpair_MemtoMem{i,j} = FCpair_MemtoMem{i,j}(tempID);
    end
end
for i = 1:size(FR_NonmemtoMem,1)
    for j = 1:size(FR_NonmemtoMem,2)
        tempID = find(FR_NonmemtoMem{i,j}<60);
        FR_NonmemtoMem{i,j} = FR_NonmemtoMem{i,j}(tempID);
        FCpair_NonmemtoMem{i,j} = FCpair_NonmemtoMem{i,j}(tempID);
    end
end

%% Unique FC pair
UniqFC_MemtoMem = GetUniqueFC(FCpair_MemtoMem);
UniqFC_NonmemtoMem = GetUniqueFC(FCpair_NonmemtoMem);

%% Effect of removing FCSP events on STM-encoding ability of following neurons
% FC of memory to memoruy pairs
AllFuAUC_MemtoMem = GetAUCofallFUsAfterRemovingFC(UniqFC_MemtoMem,ID_mPFC,mPFCSpikeTimeinIndividualTrial,mPFCUnitTrialMark,ID_aAIC,aAICSpikeTimeinIndividualTrial,aAICUnitTrialMark);
% FC of non-memory to memory pairs
AllFuAUC_NonmemtoMem = GetAUCofallFUsAfterRemovingFC(UniqFC_NonmemtoMem,ID_mPFC,mPFCSpikeTimeinIndividualTrial,mPFCUnitTrialMark,ID_aAIC,aAICSpikeTimeinIndividualTrial,aAICUnitTrialMark);

%% FC pair groups with enough number
% FCnum_MemtoMem = cellfun(@(x) size(x,1),AllFuAUC_MemtoMem);
% temp_MemtoMem = find(FCnum_MemtoMem<NumThres_FCpair);
% FCnum_NonmemtoMem = cellfun(@(x) size(x,1),AllFuAUC_NonmemtoMem);
% temp_NonmemtoMem = find(FCnum_NonmemtoMem<NumThres_FCpair);
% MaxID = min([min(temp_MemtoMem) min(temp_NonmemtoMem)]) - 1;
MaxID = 2;
AllFuAUC_MemtoMem = AllFuAUC_MemtoMem(:,1:MaxID);
AllFuAUC_NonmemtoMem = AllFuAUC_NonmemtoMem(:,1:MaxID);
FCpair_MemtoMem = FCpair_MemtoMem(:,1:MaxID);
FCpair_NonmemtoMem = FCpair_NonmemtoMem(:,1:MaxID);
ChangedValue_MemtoMem = cellfun(@(x) (x(:,2)-x(:,1))./x(:,1),AllFuAUC_MemtoMem,'UniformOutput',0);
AverChangedValue_MemtoMem = cellfun(@mean,ChangedValue_MemtoMem);
ChangedValue_NonmemtoMem = cellfun(@(x) (x(:,2)-x(:,1))./x(:,1),AllFuAUC_NonmemtoMem,'UniformOutput',0);
AverChangedValue_NonmemtoMem = cellfun(@mean,ChangedValue_NonmemtoMem);

%% Plot change in the selectivity, following removing FCSP events
if IsSemOrCiPlot == 1
    % two-way ANOVA
    data_MemtoMem = vertcat(ChangedValue_MemtoMem{:});
    data_NonmemtoMem = vertcat(ChangedValue_NonmemtoMem{:});
    data = vertcat(data_MemtoMem(:),data_NonmemtoMem(:));
    FCtype = [];
    RemovedLuNum = [];
    FCtype = vertcat(ones(numel(data_MemtoMem),1),2*ones(numel(data_NonmemtoMem),1));
    temp = vertcat(ChangedValue_MemtoMem,ChangedValue_NonmemtoMem);
    for iRow = 1:size(temp,1)
        for iColumn = 1:size(temp,2)
            RemovedLuNum = [RemovedLuNum; iColumn*ones(numel(temp{iRow,iColumn}),1)];
        end
    end
    [p,tb1] = anovan(data,{FCtype RemovedLuNum},'model',2,'varnames',{'FCtype','RemovedLUnum'},'display','on');
    % plot figure
    SemChangedValue_MemtoMem = cellfun(@(x) std(x)/sqrt(numel(x)),ChangedValue_MemtoMem);
    SemChangedValue_NonmemtoMem = cellfun(@(x) std(x)/sqrt(numel(x)),ChangedValue_NonmemtoMem);
    for i=1:numel(AverChangedValue_MemtoMem)
        bar(1.1+1.8*(i-1),AverChangedValue_MemtoMem(i),0.6,'facecolor',[0 0 0],'edgecolor',[0 0 0]); hold on
        bar(1.9+1.8*(i-1),AverChangedValue_NonmemtoMem(i),0.6,'facecolor','none','edgecolor',[0 0 0]); hold on
        errorbar(1.1+1.8*(i-1),AverChangedValue_MemtoMem(i),SemChangedValue_MemtoMem(i),'k','marker','none'); hold on
        errorbar(1.9+1.8*(i-1),AverChangedValue_NonmemtoMem(i),SemChangedValue_NonmemtoMem(i),'k','marker','none'); hold on
    end
    set(gca,'XTick',1.5:1.8:2.5+1.8*(numel(AverChangedValue_MemtoMem)-1),'XTickLabel',num2cell(1:1:6),'XLim',[0.5 2.5+1.8*(numel(AverChangedValue_MemtoMem)-1)],'YTick',-1:0.05:0,'YTickLabel',num2cell(-100:5:0),'YLim',[min(data) 0]);
    set(gcf,'Render','Painter'); saveas(gcf,sprintf('Decreased FU AUC with removing FCSP-related spikes-%s-%s',Reg,Group),'fig'); close all;
    
    %% Save result
    save(sprintf('Decreased FU AUC with removing FCSP-related spikes-%s-%s',Reg,Group),'ChangedValue_MemtoMem','ChangedValue_NonmemtoMem','FCpair_MemtoMem','FCpair_NonmemtoMem','FR_MemtoMem','FR_NonmemtoMem','tb1','-v7.3');
else
    figure('position',[300 300 300 600]);
    for i = 1:numel(AverChangedValue_MemtoMem)
        bar(1.1+1.8*(i-1),AverChangedValue_MemtoMem(i),0.6,'facecolor',[0 0 0],'edgecolor',[0 0 0]); hold on
        bar(1.9+1.8*(i-1),AverChangedValue_NonmemtoMem(i),0.6,'facecolor','none','edgecolor',[0 0 0]); hold on
        % bootstrap test
        BsTimes = 1000;
        BsSize_1 = numel(ChangedValue_MemtoMem{i});
        if BsSize_1 >=50
            BsSize_1 = 9*floor(BsSize_1/10);
        else
            BsSize_1 = BsSize_1-2*ceil(BsSize_1/25);
        end
        BsSize_2 = numel(ChangedValue_NonmemtoMem{i});
        if BsSize_2 >=50
            BsSize_2 = 9*floor(BsSize_2/10);
        else
            BsSize_2 = BsSize_2-2*ceil(BsSize_2/25);
        end
        BsValue = cell(1,2);
        for k = 1:BsTimes
            RandID = randperm(numel(ChangedValue_MemtoMem{i}));
            RandID = RandID(1:BsSize_1);
            BsValue{1,1} = [BsValue{1,1}; mean(ChangedValue_MemtoMem{i}(RandID))];
            RandID = randperm(numel(ChangedValue_NonmemtoMem{i}));
            RandID = RandID(1:BsSize_2);
            BsValue{1,2} = [BsValue{1,2}; mean(ChangedValue_NonmemtoMem{i}(RandID))];
        end
        CI_MemtoMem = prctile(BsValue{1,1},[2.5 97.5],1);
        errorbar(1.1+1.8*(i-1),AverChangedValue_MemtoMem(i),AverChangedValue_MemtoMem(i)-CI_MemtoMem(1),CI_MemtoMem(2)-AverChangedValue_MemtoMem(i),'k','marker','none'); hold on
        CI_NonmemtoMem = prctile(BsValue{1,2},[2.5 97.5],1);
        errorbar(1.9+1.8*(i-1),AverChangedValue_NonmemtoMem(i),AverChangedValue_NonmemtoMem(i)-CI_NonmemtoMem(1),CI_NonmemtoMem(2)-AverChangedValue_NonmemtoMem(i),'k','marker','none'); hold on
        diffvalue = BsValue{1,1} - BsValue{1,2};
        if mean(diffvalue) <= 0
            tempP = nnz(diffvalue>=0)/BsTimes;
        else
            tempP = nnz(diffvalue<=0)/BsTimes;
        end
        text(1.1+1.8*(i-1),AverChangedValue_MemtoMem(i)/2,['n = ' num2str(numel(ChangedValue_MemtoMem{i}))],'color','r'); hold on
        text(1.5+1.8*(i-1),AverChangedValue_MemtoMem(i),['p = ' num2str(tempP)],'color','r'); hold on
        text(1.9+1.8*(i-1),AverChangedValue_NonmemtoMem(i)/2,['n = ' num2str(numel(ChangedValue_NonmemtoMem{i}))],'color','r'); hold on
    end
    box off
    title('Bootstrap test without Bonferroni correction');
    set(gca,'XLim',[0.6 2.4+1.8*(numel(AverChangedValue_NonmemtoMem)-1)],'YTick',-1:0.05:1,'YTickLabel',num2cell(-100:5:100),'YLim',[-0.25 0]);
    set(gcf,'Render','Painter'); saveas(gcf,sprintf('Decreased FU AUC with removing FCSP-related spikes-%s-%s-1msbin',Reg,Group),'fig'); close all;
    
    %% Save result
    save(sprintf('Decreased FU AUC with removing FCSP-related spikes-%s-%s-1msbin',Reg,Group),'ChangedValue_MemtoMem','ChangedValue_NonmemtoMem','FCpair_MemtoMem','FCpair_NonmemtoMem','FR_MemtoMem','FR_NonmemtoMem','-v7.3');
end

function UniqFC = GetUniqueFC(FCpair)

UniqFC = cell(1,size(FCpair,2));
for i = 1:size(FCpair,2)
    temp = FCpair(:,i);
    temp = vertcat(temp{:});
    for j = 1:size(temp,1)
        CurrFC = temp{j};
        RestFC = temp(setdiff(1:size(temp,1),j),:);
        SameFcId = [];
        for k = 1:size(RestFC,1)
            if isequal(CurrFC,RestFC{k})
                SameFcId = [SameFcId k];
            end
        end
        if isempty(SameFcId)
            UniqFC{i} = [UniqFC{i}; {CurrFC}];
        end
    end
end
end

function SumFCpair = SumFCofDelayBins(FCpair)

SumFCpair = cell(1,size(FCpair,2));
for i = 1:size(FCpair,2)
    temp = FCpair(:,i);
    SumFCpair{i} = vertcat(temp{:});
end
end

function TarPair = GetPairAlignedWithUniqueFu(ConnChain,BinID,MaxLuNum)

FCpair = cell(1,MaxLuNum);
UniqFU = unique(ConnChain(:,2),'stable');
for iFU = 1:numel(UniqFU)
    tempLUnum = nnz(ConnChain(:,2)==UniqFU(iFU));
    for iLUnum = 1:tempLUnum
        if iLUnum <= MaxLuNum
            FCpair{iLUnum}{end+1,1} = horzcat(BinID*ones(tempLUnum,1),ConnChain(ConnChain(:,2)==UniqFU(iFU),1:2));
        end
    end
end
TarPair = cell(size(FCpair));
for i = 1:numel(FCpair) % different numbers of leading neurons
    for j = 1:size(FCpair{i},1) % based on same following neurons
        temp = nchoosek(1:size(FCpair{i}{j},1),i);
        for k = 1:size(temp,1)
            TarPair{i}{end+1,1} = FCpair{i}{j}(temp(k,:),:);
        end
    end
end
end

function FR = ComputeAllLuFR(FCpair,Reg,mPFCUnitsID,mPFCUnitsSpikeTime,mPFCUnitsTrialMark,aAICUnitsID,aAICUnitsSpikeTime,aAICUnitsTrialMark,BaseLen,SampLen)

FR = cell(size(FCpair));
for i = 1:size(FCpair,2) % different leading neurons number
    for j = 1:size(FCpair,1) % Go or NoGo trials
        for k = 1:size(FCpair{j,i},1)
            BinID = FCpair{j,i}{k}(1,1);
            ID_LU = FCpair{j,i}{k}(:,2);
            IDinReg = [];
            if strcmp(Reg,'Within mPFC') || strcmp(Reg,'mPFC-aAIC')
                for iUnit = 1:numel(ID_LU)
                    IDinReg = [IDinReg; find(mPFCUnitsID==ID_LU(iUnit))];
                end
                SpikeTime = mPFCUnitsSpikeTime;
                TrialMark = mPFCUnitsTrialMark;
            elseif strcmp(Reg,'Within aAIC')
                for iUnit = 1:numel(ID_LU)
                    IDinReg = [IDinReg; find(aAICUnitsID==ID_LU(iUnit))];
                end
                SpikeTime = aAICUnitsSpikeTime;
                TrialMark = aAICUnitsTrialMark;
            end
            tempMaxFR = UnitsMaxFR(j,BaseLen,SampLen,BinID,IDinReg,SpikeTime,TrialMark);
            FR{j,i} = [FR{j,i}; tempMaxFR];
        end
    end
end
end

function MaxFR = UnitsMaxFR(TrialType,BaseLen,SampLen,DelayBinID,UnitsID,SpikeRasterTime,TrialMark)

TrialMark = TrialMark{UnitsID(1)};
MaxFR = [];
for iUnit = 1:numel(UnitsID)
    tempUnitID = UnitsID(iUnit);
    if TrialType == 1
        tempTrialID = find(TrialMark(:,2)==1);
    else
        tempTrialID = find(TrialMark(:,2)==2);
    end
    tempSpike = SpikeRasterTime{tempUnitID}(:,tempTrialID);
    SuFR = [];
    for iTrial = 1:numel(tempSpike)
        TempSingleTrialRaster = tempSpike{iTrial};
        TempSingleTrialRaster(TempSingleTrialRaster<BaseLen+SampLen+DelayBinID-1 | TempSingleTrialRaster>=BaseLen+SampLen+DelayBinID) = [];
        SuFR = [SuFR; numel(TempSingleTrialRaster)];
    end
    MaxFR = [MaxFR; mean(SuFR)];
end
MaxFR = max(MaxFR);
end



function AllFuAuc = GetAUCofallFUsAfterRemovingFC(FCpair,mPFCUnitsID,mPFCUnitsSpikeTime,mPFCUnitsTrialMark,aAICUnitsID,aAICUnitsSpikeTime,aAICUnitsTrialMark)
AllFuAuc = cell(size(FCpair));
for i = 1:numel(FCpair) % different numbers of leading neurons
    for j = 1:size(FCpair{i},1) % based on same following neurons
        fprintf('Process %dth one of total %d FC pairs with %d leading neurons\n',j,size(FCpair{i},1),i);
        ID_Bin = FCpair{i}{j}(1,1);
        ID_LU = FCpair{i}{j}(:,2);
        ID_FU = FCpair{i}{j}(1,3);
        % judge region ID of leading and following neurons
        IDinReg_LU = [];
        if ismember(ID_LU(1,1),mPFCUnitsID) % leading neuron
            SpikeInfo_LU = mPFCUnitsSpikeTime;
            for iUnit = 1:numel(ID_LU)
                IDinReg_LU = [IDinReg_LU; find(mPFCUnitsID==ID_LU(iUnit))];
            end
        else
            SpikeInfo_LU = aAICUnitsSpikeTime;
            for iUnit = 1:numel(ID_LU)
                IDinReg_LU = [IDinReg_LU; find(aAICUnitsID==ID_LU(iUnit))];
            end
        end
        if ismember(ID_FU,mPFCUnitsID) % following neuron
            SpikeInfo_FU = mPFCUnitsSpikeTime;
            IDinReg_FU = find(mPFCUnitsID==ID_FU);
            TrialMark = mPFCUnitsTrialMark{IDinReg_FU};
        else
            SpikeInfo_FU = aAICUnitsSpikeTime;
            IDinReg_FU = find(aAICUnitsID==ID_FU);
            TrialMark = aAICUnitsTrialMark{IDinReg_FU};
        end
        [tempAUC_wiFC,tempAUC_woFC] = DecreaFUauROCofRemoveFCSP(ID_Bin,IDinReg_LU,IDinReg_FU,SpikeInfo_LU,SpikeInfo_FU,TrialMark);
        AllFuAuc{i} = [AllFuAuc{i}; horzcat(tempAUC_wiFC,tempAUC_woFC)];
    end
end
end

function [AUC_wiFC,AUC_woFC] = DecreaFUauROCofRemoveFCSP(BinID,IDinReg_PreUnit,IDinReg_PostUnit,RG_PreUnit,RG_PostUnit,TrialMark)

BaseLen = 4;
OdorLen = 1;
FR_wiFC_Go = [];
FR_woFC_Go = [];
FR_wiFC_NoGo = [];
FR_woFC_NoGo = [];

% raster timestamps of target units
temp = RG_PreUnit(IDinReg_PreUnit,:);
RG_PreUnit = [];
for i = 1:numel(temp)
    RG_PreUnit = [RG_PreUnit; temp{i}];
end
temp = RG_PostUnit(IDinReg_PostUnit,:);
RG_PostUnit = [];
for i = 1:numel(temp)
    RG_PostUnit = [RG_PostUnit; temp{i}];
end

% trial ID of Go and NoGo trials
TrialID_Go = find(TrialMark(:,2)==1);
TrialID_NoGo = find(TrialMark(:,2)==2);
% raster timestamps of target units in Go and NoGo trials
RG_PreUnit_Go = RG_PreUnit(:,TrialID_Go);
RG_PreUnit_NoGo = RG_PreUnit(:,TrialID_NoGo);
RG_PostUnit_Go = RG_PostUnit(:,TrialID_Go);
RG_PostUnit_NoGo = RG_PostUnit(:,TrialID_NoGo);
for iTrialType = 1:2
    % removing FCSP-related spikes from activities of following neurons in Go (iTrialType==1) and NoGo (iTrialType==2) trials
    if iTrialType == 1
        RG_PostUnit = RG_PostUnit_Go;
        RG_PreUnit = RG_PreUnit_Go;
    elseif iTrialType == 2
        RG_PostUnit = RG_PostUnit_NoGo;
        RG_PreUnit = RG_PreUnit_NoGo;
    end
    for iTrial = 1:numel(RG_PostUnit)
        % raster timestamp of postunit
        tempRG_PostUnit = RG_PostUnit(:,iTrial);
        tempRG_PostUnit = vertcat(tempRG_PostUnit{:});
        tempRG_PostUnit = tempRG_PostUnit(tempRG_PostUnit>BaseLen+OdorLen+BinID-1 & tempRG_PostUnit<=BaseLen+OdorLen+BinID);
        % raster timestamp of preunit
        tempRG_PreUnit = RG_PreUnit(:,iTrial);
        tempRG_PreUnit = vertcat(tempRG_PreUnit{:});
        tempRG_PreUnit = tempRG_PreUnit(tempRG_PreUnit>=BaseLen+OdorLen+BinID-1 & tempRG_PreUnit<BaseLen+OdorLen+BinID);
        % removing FCSP-related spikes of following neurons
        KeptSpkId = [];
        for iSpike = 1:numel(tempRG_PostUnit)
            temp = find(tempRG_PostUnit(iSpike)-tempRG_PreUnit<=0.005 & tempRG_PostUnit(iSpike)-tempRG_PreUnit>0.002);
            if isempty(temp)
                KeptSpkId = [KeptSpkId; iSpike];
            end
        end
        if iTrialType == 1
            FR_wiFC_Go = [FR_wiFC_Go; numel(tempRG_PostUnit)];
            FR_woFC_Go = [FR_woFC_Go; numel(tempRG_PostUnit(KeptSpkId))];
        elseif iTrialType == 2
            FR_wiFC_NoGo = [FR_wiFC_NoGo; numel(tempRG_PostUnit)];
            FR_woFC_NoGo = [FR_woFC_NoGo; numel(tempRG_PostUnit(KeptSpkId))];
        end
    end
end
AUC_wiFC = CalculateSingleUnitAUCinOneBin(FR_wiFC_Go,FR_wiFC_NoGo);
AUC_woFC = CalculateSingleUnitAUCinOneBin(FR_woFC_Go,FR_woFC_NoGo);
AUC_wiFC = 2*(max([0.5 AUC_wiFC])-0.5);
AUC_woFC = 2*(max([0.5 AUC_woFC])-0.5);
end

function SingleUnitAUCvalue = CalculateSingleUnitAUCinOneBin(FR_S1,FR_S2)

% calculate decision variable (DV)
[DV_S1,DV_S2] = DecisionVariableCalculation(FR_S1,FR_S2);
% FPR and TPR
[TPR,FPR] = TprFprCalculation(DV_S1,DV_S2);
% calculate auROC
SingleUnitAUCvalue = AucAnalysis(FPR,TPR);
end
