
%% Find example FC neuronal pair to show that FCSP events transfer information from leading to following neurons.

clear; clc; close all;

%% Assignment
Group = 'CtrlGroup';
DelayBinID = 1:1:4;
FCContext = 'hit'; % 'Hit': FC pair in hit trials; 'CR': in CR trials; 'Combine': combine pair in Hit and CR trials.
BaseLen = 1; BaseLen_RG = 4; OdorLen = 1; DelayLen = 4;
Thres_Selec = 0.20;
Thres_DecrSelec = -0.01;
Thres_PropDecrSelec = -0.01;
Thres_PropFcsp = 0.05;
ShownTrialNum = 2;

%% Load FC result, get target FC
FCpair_MemtoMem = [];
FCpair_NonmemtoMem = [];
for i = 1:numel(DelayBinID)
    load(sprintf('FCofNeurons_PropertyPopulation_%d_%d_%s.mat',DelayBinID(i),DelayBinID(i)+1,Group));
    temp_Go = ismember(reg_chain_go(:,1),[1 2]) & ismember(reg_chain_go(:,2),[1 2]);
    temp_NoGo = ismember(reg_chain_nogo(:,1),[1 2]) & ismember(reg_chain_nogo(:,2),[1 2]);
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
    
    %% Target FC
    if strcmp(FCContext,'hit')
        conn_chain = conn_chain_go;
        pref_chain = pref_chain_go;
    elseif strcmp(FCContext,'CR')
        conn_chain = conn_chain_nogo;
        pref_chain = pref_chain_nogo;
    elseif strcmp(FCContext,'Combine')
        conn_chain = vertcat(conn_chain_go,conn_chain_nogo);
        pref_chain = vertcat(pref_chain_go,pref_chain_nogo);
        [~,ia,~] = unique(conn_chain(:,1:2),'rows','stable');
        conn_chain = conn_chain(ia,:);
        pref_chain = pref_chain(ia,:);
    end
    
    %% Divide FC into memory to memory and non-memory to memory neuronal pairs
    conn_chain_MemtoMem = conn_chain(pref_chain(:,BaseLen+OdorLen+DelayBinID(i))>0 & pref_chain(:,BaseLen+OdorLen+DelayBinID(i)+size(pref_chain,2)/2)>0,1:2);
    conn_chain_NonmemtoMem = conn_chain(pref_chain(:,BaseLen+OdorLen+DelayBinID(i))==0 & pref_chain(:,BaseLen+OdorLen+DelayBinID(i)+size(pref_chain,2)/2)>0,1:2);
    FCpair_MemtoMem = [FCpair_MemtoMem; horzcat(DelayBinID(i)*ones(size(conn_chain_MemtoMem,1),1),conn_chain_MemtoMem)];
    FCpair_NonmemtoMem = [FCpair_NonmemtoMem; horzcat(DelayBinID(i)*ones(size(conn_chain_NonmemtoMem,1),1),conn_chain_NonmemtoMem)];
end

%% Load ID of mPFC and aAIC neurons
load(sprintf('AllNeuronsAUCValue_%s.mat',Group));

%% Load trial-based timestamps of spike rasters for mPFC and aAIC neurons
load(sprintf('NeuronOriginandSpikeRasterInformationfor%s.mat',Group));

%% Selectivity of following neurons with and without FCSP, and proportion of FCSP-related spikes in spikes of following neurons
% FC of memory to memory pairs
[AllFuSelec_MemtoMem,AllFuPropFcsp_MemtoMem] = GetSelecOfAllFuAfterRemovingFcsp(FCpair_MemtoMem,ID_mPFC,mPFCSpikeTimeinIndividualTrial,mPFCUnitTrialMark,ID_aAIC,aAICSpikeTimeinIndividualTrial,...
    aAICUnitTrialMark);
SelecDiff_MemtoMem = abs(AllFuSelec_MemtoMem(:,2)) - abs(AllFuSelec_MemtoMem(:,1));
PropSelecDiff_MemtoMem = SelecDiff_MemtoMem./abs(AllFuSelec_MemtoMem(:,1));
% FC of non-memory to memory pairs
[AllFuSelec_NonmemtoMem,AllFuPropFcsp_NonmemtoMem] = GetSelecOfAllFuAfterRemovingFcsp(FCpair_NonmemtoMem,ID_mPFC,mPFCSpikeTimeinIndividualTrial,mPFCUnitTrialMark,ID_aAIC,aAICSpikeTimeinIndividualTrial,...
    aAICUnitTrialMark);
SelecDiff_NonmemtoMem = abs(AllFuSelec_NonmemtoMem(:,2)) - abs(AllFuSelec_NonmemtoMem(:,1));
PropSelecDiff_NonmemtoMem = SelecDiff_NonmemtoMem./abs(AllFuSelec_NonmemtoMem(:,1));

%% Example FC pairs meeting demands
% memory to memory neuronal pairs
tempID_MemtoMem = find(AllFuSelec_MemtoMem(:,1)>=Thres_Selec & SelecDiff_MemtoMem<=Thres_DecrSelec & PropSelecDiff_MemtoMem<=Thres_PropDecrSelec & AllFuPropFcsp_MemtoMem>=Thres_PropFcsp);
FCpair_MemtoMem = FCpair_MemtoMem(tempID_MemtoMem,:);
AllFuSelec_MemtoMem = AllFuSelec_MemtoMem(tempID_MemtoMem,:);
AllFuPropFcsp_MemtoMem = AllFuPropFcsp_MemtoMem(tempID_MemtoMem,:);
SelecDiff_MemtoMem = SelecDiff_MemtoMem(tempID_MemtoMem,:);
PropSelecDiff_MemtoMem = PropSelecDiff_MemtoMem(tempID_MemtoMem,:);
% non-memory to memory neuronal pairs
tempID_NonmemtoMem = find(AllFuSelec_NonmemtoMem(:,1)>=Thres_Selec & SelecDiff_NonmemtoMem<=Thres_DecrSelec & PropSelecDiff_NonmemtoMem<=Thres_PropDecrSelec & AllFuPropFcsp_NonmemtoMem>=Thres_PropFcsp);
FCpair_NonmemtoMem = FCpair_NonmemtoMem(tempID_NonmemtoMem,:);
AllFuSelec_NonmemtoMem = AllFuSelec_NonmemtoMem(tempID_NonmemtoMem,:);
AllFuPropFcsp_NonmemtoMem = AllFuPropFcsp_NonmemtoMem(tempID_NonmemtoMem,:);
SelecDiff_NonmemtoMem = SelecDiff_NonmemtoMem(tempID_NonmemtoMem,:);
PropSelecDiff_NonmemtoMem = PropSelecDiff_NonmemtoMem(tempID_NonmemtoMem,:);

%% Raster plot of two neurons in example FC pairs
% coding to coding
PlotPopulationPairsSpikeRaster('MemtoMem',FCpair_MemtoMem,AllFuSelec_MemtoMem,AllFuPropFcsp_MemtoMem,PropSelecDiff_MemtoMem,ID_mPFC,mPFCSpikeTimeinIndividualTrial,...
    mPFCUnitTrialMark,ID_aAIC,aAICSpikeTimeinIndividualTrial,aAICUnitTrialMark,BaseLen_RG,OdorLen,ShownTrialNum);
% non-coding to coding
PlotPopulationPairsSpikeRaster('NonmemtoMem',FCpair_NonmemtoMem,AllFuSelec_NonmemtoMem,AllFuPropFcsp_NonmemtoMem,PropSelecDiff_NonmemtoMem,ID_mPFC,mPFCSpikeTimeinIndividualTrial,...
    mPFCUnitTrialMark,ID_aAIC,aAICSpikeTimeinIndividualTrial,aAICUnitTrialMark,BaseLen_RG,OdorLen,ShownTrialNum);

function PlotPopulationPairsSpikeRaster(PairType,Pair,AllFuSelec,AllFuFcspProp,PropofDecrSelec,ID_mPFC,mPFCSpikeTime,mPFCTrialMark,ID_aAIC,aAICSpikeTime,aAICTrialMark,BaseLen,OdorLen,ShownTrialNum)
for iPair = 1:size(Pair,1)
    tempBinID = Pair(iPair,1);
    tempLU = Pair(iPair,2);
    tempFU = Pair(iPair,3);
    tempFuSelec_wiFCSP = AllFuSelec(iPair,1);
    tempFuSelec_woFCSP = AllFuSelec(iPair,2);
    tempFuFcspProp = AllFuFcspProp(iPair,1);
    tempDecSelecProp = PropofDecrSelec(iPair,1);
    if ismember(tempLU,ID_mPFC)
        tempReg_LU = 'mPFC';
        tempIDinReg_LU = find(ID_mPFC==tempLU);
        tempRG_LU = mPFCSpikeTime{tempIDinReg_LU};
    elseif ismember(tempLU,ID_aAIC)
        tempReg_LU = 'aAIC';
        tempIDinReg_LU = find(ID_aAIC==tempLU);
        tempRG_LU = aAICSpikeTime{tempIDinReg_LU};
    end
    if ismember(tempFU,ID_mPFC)
        tempReg_FU = 'mPFC';
        tempIDinReg_FU = find(ID_mPFC==tempFU);
        tempRG_FU = mPFCSpikeTime{tempIDinReg_FU};
        tempTrialMark = mPFCTrialMark{tempIDinReg_FU};
    elseif ismember(tempFU,ID_aAIC)
        tempReg_FU = 'aAIC';
        tempIDinReg_FU = find(ID_aAIC==tempFU);
        tempRG_FU = aAICSpikeTime{tempIDinReg_FU};
        tempTrialMark = aAICTrialMark{tempIDinReg_FU};
    end
    tempID_GoTrial = find(tempTrialMark(:,2)==1);
    tempRG_Go_LU = tempRG_LU(:,tempID_GoTrial);
    tempRG_Go_FU = tempRG_FU(:,tempID_GoTrial);
    tempID_NoGoTrial = find(tempTrialMark(:,2)==2);
    tempRG_NoGo_LU = tempRG_LU(:,tempID_NoGoTrial);
    tempRG_NoGo_FU = tempRG_FU(:,tempID_NoGoTrial);
    PlotIndividualPairSpikeRaster(tempReg_LU,tempIDinReg_LU,tempRG_Go_LU,tempRG_NoGo_LU,tempReg_FU,tempIDinReg_FU,tempRG_Go_FU,tempRG_NoGo_FU,PairType,...
        tempFuSelec_wiFCSP,tempFuSelec_woFCSP,tempFuFcspProp,tempDecSelecProp,BaseLen,OdorLen,tempBinID,ShownTrialNum);
end
end

function PlotIndividualPairSpikeRaster(Reg_PreUnit,IDinReg_PreUnit,SpikeTime_PreUnit_S1,SpikeTime_PreUnit_S2,Reg_PostUnit,IDinReg_PostUnit,SpikeTime_PostUnit_S1,...
    SpikeTime_PostUnit_S2,PairType,PostUnitSelec_wiFCSP,PostUnitSelec_woFCSP,PostUnitFcspProp,DecSelecProp,BaseLen,SampLen,BinID,ShownTrialNum)

C1 = [0.5 0.5 0.5];
C2 = [220 39 114]/255;
C3 = [39 148 211]/255;
FcPeriod = BaseLen+SampLen+BinID;

%% FR in S1 trials
FcspNum_S1 = [];
FR_wiFCSP_S1 = [];
FR_woFCSP_S1 = [];
for iTrial = 1:length(SpikeTime_PreUnit_S1)
    FCspiketime = [];
    for iSpike = 1:length(SpikeTime_PreUnit_S1{iTrial})
        if ~isempty(find(SpikeTime_PostUnit_S1{iTrial}(:,1)-SpikeTime_PreUnit_S1{iTrial}(iSpike)<=0.01 & SpikeTime_PostUnit_S1{iTrial}(:,1)-SpikeTime_PreUnit_S1{iTrial}(iSpike)>0.002))
            FCspiketime = [FCspiketime SpikeTime_PreUnit_S1{iTrial}(iSpike)];
        end
    end
    FCspiketime(FCspiketime<FcPeriod-1 | FCspiketime>=FcPeriod) = [];
    FcspNum_S1 = [FcspNum_S1; numel(FCspiketime)];
    KeptSpkId = [];
    for iSpike = 1:length(SpikeTime_PostUnit_S1{iTrial})
        temp = find(SpikeTime_PostUnit_S1{iTrial}(iSpike)-SpikeTime_PreUnit_S1{iTrial}<=0.01 & SpikeTime_PostUnit_S1{iTrial}(iSpike)-SpikeTime_PreUnit_S1{iTrial}>0.002);
        if isempty(temp)
            KeptSpkId = [KeptSpkId; iSpike];
        end
    end
    tempRG_wiFCSP = SpikeTime_PostUnit_S1{iTrial};
    tempRG_woFCSP = tempRG_wiFCSP(KeptSpkId);
    tempRG_wiFCSP(tempRG_wiFCSP <= FcPeriod-1 | tempRG_wiFCSP > FcPeriod) = [];
    tempRG_woFCSP(tempRG_woFCSP <= FcPeriod-1 | tempRG_woFCSP > FcPeriod) = [];
    FR_wiFCSP_S1 = [FR_wiFCSP_S1; numel(tempRG_wiFCSP)];
    FR_woFCSP_S1 = [FR_woFCSP_S1; numel(tempRG_woFCSP)];
end

%% FR in S2 trials
FR_wiFCSP_S2 = [];
FR_woFCSP_S2 = [];
for iTrial = 1:length(SpikeTime_PreUnit_S2)
    KeptSpkId = [];
    for iSpike = 1:length(SpikeTime_PostUnit_S2{iTrial})
        temp = find(SpikeTime_PostUnit_S2{iTrial}(iSpike)-SpikeTime_PreUnit_S2{iTrial}<=0.01 & SpikeTime_PostUnit_S2{iTrial}(iSpike)-SpikeTime_PreUnit_S2{iTrial}>0.002);
        if isempty(temp)
            KeptSpkId = [KeptSpkId; iSpike];
        end
    end
    tempRG_wiFCSP = SpikeTime_PostUnit_S2{iTrial};
    tempRG_woFCSP = tempRG_wiFCSP(KeptSpkId);
    tempRG_wiFCSP(tempRG_wiFCSP<=FcPeriod-1 | tempRG_wiFCSP>FcPeriod) = [];
    tempRG_woFCSP(tempRG_woFCSP<=FcPeriod-1 | tempRG_woFCSP>FcPeriod) = [];
    FR_wiFCSP_S2 = [FR_wiFCSP_S2; numel(tempRG_wiFCSP)];
    FR_woFCSP_S2 = [FR_woFCSP_S2; numel(tempRG_woFCSP)];
end

%% Rastergram
if nnz(FcspNum_S1>=2) >= ShownTrialNum
    % ID of S1 trials
    [~,sortid] = sortrows(FcspNum_S1,'descend');
    TrialID_S1 = sortid(1:ShownTrialNum);
    % ID of S2 trials
    TrialID_S2 = randperm(numel(SpikeTime_PreUnit_S2));
    TrialID_S2 = TrialID_S2(1:ShownTrialNum);
    % spike raster
    figure('Color','w','Position',[500,200,300,700])
    SpikeData_PreUnit = horzcat(SpikeTime_PreUnit_S1(:,TrialID_S1),SpikeTime_PreUnit_S2(:,TrialID_S2));
    SpikeData_PostUnit = horzcat(SpikeTime_PostUnit_S1(:,TrialID_S1),SpikeTime_PostUnit_S2(:,TrialID_S2));
    SpikeData = vertcat(SpikeData_PreUnit,SpikeData_PostUnit);
    for iTrial = 1:size(SpikeData,2)
        for iSpike = 1:numel(SpikeData{1,iTrial}) % preunit
            if ~isempty(find(SpikeData{2,iTrial}-SpikeData{1,iTrial}(iSpike)<=0.01 & SpikeData{2,iTrial}-SpikeData{1,iTrial}(iSpike)>0.002))
                plot([SpikeData{1,iTrial}(iSpike) SpikeData{1,iTrial}(iSpike)], [1.5+3*(size(SpikeData,2)-iTrial) 2.5+3*(size(SpikeData,2)-iTrial)],'color',C2);
            else
                plot([SpikeData{1,iTrial}(iSpike) SpikeData{1,iTrial}(iSpike)], [1.5+3*(size(SpikeData,2)-iTrial) 2.5+3*(size(SpikeData,2)-iTrial)],'color',C1);
            end
            hold on
        end
        for iSpike = 1:numel(SpikeData{2,iTrial}) % postunit
            if ~isempty(find(SpikeData{2,iTrial}(iSpike)-SpikeData{1,iTrial}<=0.01 & SpikeData{2,iTrial}(iSpike)-SpikeData{1,iTrial}>0.002))
                plot([SpikeData{2,iTrial}(iSpike) SpikeData{2,iTrial}(iSpike)], [0.5+3*(size(SpikeData,2)-iTrial) 1.5+3*(size(SpikeData,2)-iTrial)],'color',C3);
            else
                plot([SpikeData{2,iTrial}(iSpike) SpikeData{2,iTrial}(iSpike)], [0.5+3*(size(SpikeData,2)-iTrial) 1.5+3*(size(SpikeData,2)-iTrial)],'color',C1);
            end
            hold on
        end
        text(FcPeriod-0.2,1.5+3*(size(SpikeData,2)-iTrial),num2str(nnz(SpikeData{2,iTrial}>FcPeriod-1 & SpikeData{2,iTrial}<=FcPeriod)),'color','red'); hold on
    end
    box off;
    title(sprintf('Selectivity with FCSP is %d, FR_S1 = %d, FR_S2 = %d\n selectivity without FCSP is %d, FR_S1_woFCSP = %d, FR_S2_woFCSP = %d\n selectivity change is %d, FCSP proportion is %d',PostUnitSelec_wiFCSP,mean(FR_wiFCSP_S1),mean(FR_wiFCSP_S2),...
        PostUnitSelec_woFCSP,mean(FR_woFCSP_S1),mean(FR_woFCSP_S2),DecSelecProp,PostUnitFcspProp));
    set(gca,'XTick',FcPeriod-1:1:FcPeriod,'XTickLabel',{'','','','',''},'YTick',1.5:3:3*size(SpikeData,2),'YTickLabel',num2cell(1:1:size(SpikeData,2)),'FontName','Arial','FontSize',16,'XLim',[FcPeriod-1 FcPeriod],'YLim',[0 3*size(SpikeData,2)]);
    set(gcf,'Render','Painter'); saveas(gcf,sprintf('Remove FCSP case-%s-%s#%d-%s#%d_at delay %dth s',PairType,Reg_PreUnit,IDinReg_PreUnit,Reg_PostUnit,IDinReg_PostUnit,BinID),'fig'); close all;
end
end

function [FuSelec,FuPropFcsp] = GetSelecOfAllFuAfterRemovingFcsp(FCpair,mPFCUnitsID,mPFCUnitsSpikeTime,mPFCUnitsTrialMark,aAICUnitsID,aAICUnitsSpikeTime,aAICUnitsTrialMark)
FuSelec = [];
FuPropFcsp = [];
for iPair = 1:size(FCpair,1) 
    fprintf('Analyze %dth one of total %d FC pairs\n',iPair,size(FCpair,1));
    ID_Bin = FCpair(iPair,1);
    ID_LU = FCpair(iPair,2);
    ID_FU = FCpair(iPair,3);
    % regions of leading and following neurons
    if ismember(ID_LU,mPFCUnitsID) % leading neuron
        SpikeInfo_LU = mPFCUnitsSpikeTime;
        IDinReg_LU = find(mPFCUnitsID==ID_LU);
    else
        SpikeInfo_LU = aAICUnitsSpikeTime;
        IDinReg_LU = find(aAICUnitsID==ID_LU);
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
    [tempSelec_wiFC,tempSelec_woFC,tempFcspProp] = CalculateFuDecrSelec(ID_Bin,IDinReg_LU,IDinReg_FU,SpikeInfo_LU,SpikeInfo_FU,TrialMark);
    FuSelec = [FuSelec; horzcat(tempSelec_wiFC,tempSelec_woFC)];
    FuPropFcsp = [FuPropFcsp; tempFcspProp];
end
end

function [Selec_wiFC,Selec_woFC,FcspProp] = CalculateFuDecrSelec(BinID,IDinReg_PreUnit,IDinReg_PostUnit,RG_PreUnit,RG_PostUnit,TrialMark)

BaseLen = 4;
OdorLen = 1;
FR_wiFC_Go = [];
FR_woFC_Go = [];
FR_wiFC_NoGo = [];
FR_woFC_NoGo = [];
FcspProp = [];

%% Raster timestamps of target neurons
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

%% ID of Go and NoGo trials
TrialID_Go = find(TrialMark(:,2)==1);
TrialID_NoGo = find(TrialMark(:,2)==2);
% raster timestamps of target neurons in Go and NoGo trials
RG_PreUnit_Go = RG_PreUnit(:,TrialID_Go);
RG_PreUnit_NoGo = RG_PreUnit(:,TrialID_NoGo);
RG_PostUnit_Go = RG_PostUnit(:,TrialID_Go);
RG_PostUnit_NoGo = RG_PostUnit(:,TrialID_NoGo);
for iTrialType = 1:2
    % remove FCSP-related spikes from spikes of following neurons in Go (iTrialType=1) and NoGo (iTrialType=2) trials
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
        % remove FCSP-related spikes of postunit
        KeptSpkId = [];
        for iSpike = 1:numel(tempRG_PostUnit)
            temp = find(tempRG_PostUnit(iSpike)-tempRG_PreUnit<=0.01 & tempRG_PostUnit(iSpike)-tempRG_PreUnit>0.002);
            if isempty(temp)
                KeptSpkId = [KeptSpkId; iSpike];
            end
        end
        if iTrialType == 1
            FR_wiFC_Go = [FR_wiFC_Go; numel(tempRG_PostUnit)];
            FR_woFC_Go = [FR_woFC_Go; numel(tempRG_PostUnit(KeptSpkId))];
            if numel(tempRG_PostUnit) > 0
                FcspProp = [FcspProp; (numel(tempRG_PostUnit)-numel(tempRG_PostUnit(KeptSpkId)))/numel(tempRG_PostUnit)];
            elseif numel(tempRG_PostUnit) == 0
                FcspProp = [FcspProp; 0];
            end
        elseif iTrialType == 2
            FR_wiFC_NoGo = [FR_wiFC_NoGo; numel(tempRG_PostUnit)];
            FR_woFC_NoGo = [FR_woFC_NoGo; numel(tempRG_PostUnit(KeptSpkId))];
        end
    end
end
Selec_wiFC = (mean(FR_wiFC_Go)-mean(FR_wiFC_NoGo))/(mean(FR_wiFC_Go)+mean(FR_wiFC_NoGo));
Selec_woFC = (mean(FR_woFC_Go)-mean(FR_woFC_NoGo))/(mean(FR_woFC_Go)+mean(FR_woFC_NoGo));
FcspProp = mean(FcspProp);
end
