%% Effect of removing FCSP events on STM-encoding ability of following neurons
%% /// decreased STM-encoding ability following removing FCSP events ///
%% /// controlling for FR of leading neurons ///

clear; clc; close all;

%% Assignment
Group = 'CtrlGroup';
Reg = 'Within mPFC';
FRrange = [8 10];
NumThres_FCpair = 5;
IsPlot = 0;

%% Load result
load(sprintf('Decreased FU selectivity with removing FCSP-related spikes-%s-%s',Reg,Group));
load(sprintf('AllNeuronsAUCValue_%s.mat',Group));
load(sprintf('NeuronOriginandSpikeRasterInformationfor%s.mat',Group));

%% FC pairs with FR of LU in the range
MaxNum_LU = min(horzcat(size(FR_MemtoMem,2),size(FR_NonmemtoMem,2),size(FCpair_MemtoMem,2),size(FCpair_NonmemtoMem,2)));
FR_MemtoMem = FR_MemtoMem(:,1:MaxNum_LU);
FR_NonmemtoMem = FR_NonmemtoMem(:,1:MaxNum_LU);
FCpair_MemtoMem = FCpair_MemtoMem(:,1:MaxNum_LU);
FCpair_NonmemtoMem = FCpair_NonmemtoMem(:,1:MaxNum_LU);
for i = 1:size(FR_MemtoMem,1)
    for j = 1:size(FR_MemtoMem,2)
        temp = find(FR_MemtoMem{i,j}>=FRrange(1) & FR_MemtoMem{i,j}<=FRrange(2));
        FR_MemtoMem{i,j} = FR_MemtoMem{i,j}(temp);
        FCpair_MemtoMem{i,j} = FCpair_MemtoMem{i,j}(temp);
        temp = find(FR_NonmemtoMem{i,j}>=FRrange(1) & FR_NonmemtoMem{i,j}<=FRrange(2));
        FR_NonmemtoMem{i,j} = FR_NonmemtoMem{i,j}(temp);
        FCpair_NonmemtoMem{i,j} = FCpair_NonmemtoMem{i,j}(temp);
    end
end

%% Unique FC neuronal pair
[UniqFC_MemtoMem,UniqFR_MemtoMem] = GetUniqueFC(FCpair_MemtoMem,FR_MemtoMem);
[UniqFC_NonmemtoMem,UniqFR_NonmemtoMem] = GetUniqueFC(FCpair_NonmemtoMem,FR_NonmemtoMem);

%% Effect of removing FCSP events on STM-encoding ability of following neurons
% FC neuronal pairs: memory to memory
AllFuSelec_MemtoMem = GetSelecOfFuWithRemovingFC(UniqFC_MemtoMem,ID_mPFC,mPFCSpikeTimeinIndividualTrial,mPFCUnitTrialMark,ID_aAIC,aAICSpikeTimeinIndividualTrial,...
    aAICUnitTrialMark);
% FC neuronal pairs: non-memory to memory
AllFuSelec_NonmemtoMem = GetSelecOfFuWithRemovingFC(UniqFC_NonmemtoMem,ID_mPFC,mPFCSpikeTimeinIndividualTrial,mPFCUnitTrialMark,ID_aAIC,aAICSpikeTimeinIndividualTrial,...
    aAICUnitTrialMark);

%% Plot change in the STM-encoding ability, following removing FCSP events
% FC pair groups with enough number
FCnum_MemtoMem = cellfun(@(x) size(x,1),AllFuSelec_MemtoMem);
temp_MemtoMem = find(FCnum_MemtoMem<NumThres_FCpair);
FCnum_NonmemtoMem = cellfun(@(x) size(x,1),AllFuSelec_NonmemtoMem);
temp_NonmemtoMem = find(FCnum_NonmemtoMem<NumThres_FCpair);
MaxID = min([min(temp_MemtoMem) min(temp_NonmemtoMem)]);
if MaxID >= 2
    MaxID = MaxID - 1;
else
    MaxID = 1;
end
AllFuSelec_MemtoMem = AllFuSelec_MemtoMem(:,1:MaxID);
AllFuSelec_NonmemtoMem = AllFuSelec_NonmemtoMem(:,1:MaxID);
UniqFC_MemtoMem = UniqFC_MemtoMem(:,1:MaxID);
UniqFC_NonmemtoMem = UniqFC_NonmemtoMem(:,1:MaxID);
UniqFR_MemtoMem = UniqFR_MemtoMem(:,1:MaxID);
UniqFR_NonmemtoMem = UniqFR_NonmemtoMem(:,1:MaxID);
ChangedValue_MemtoMem = cellfun(@(x) (abs(x(:,2))-abs(x(:,1)))./abs(x(:,1)),AllFuSelec_MemtoMem,'UniformOutput',0);
ChangedValue_NonmemtoMem = cellfun(@(x) (abs(x(:,2))-abs(x(:,1)))./abs(x(:,1)),AllFuSelec_NonmemtoMem,'UniformOutput',0);
UniqFC_MemtoMem = cellfun(@(x,y) x(y<=1,:),UniqFC_MemtoMem,ChangedValue_MemtoMem,'UniformOutput',false);
UniqFR_MemtoMem = cellfun(@(x,y) x(y<=1,:),UniqFR_MemtoMem,ChangedValue_MemtoMem,'UniformOutput',false);
ChangedValue_MemtoMem = cellfun(@(x) x(x<=1,:),ChangedValue_MemtoMem,'UniformOutput',false);
UniqFC_NonmemtoMem = cellfun(@(x,y) x(y<=1,:),UniqFC_NonmemtoMem,ChangedValue_NonmemtoMem,'UniformOutput',false);
UniqFR_NonmemtoMem = cellfun(@(x,y) x(y<=1,:),UniqFR_NonmemtoMem,ChangedValue_NonmemtoMem,'UniformOutput',false);
ChangedValue_NonmemtoMem = cellfun(@(x) x(x<=1,:),ChangedValue_NonmemtoMem,'UniformOutput',false);
if IsPlot == 1
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
    [p,tbl] = anovan(data,{FCtype RemovedLuNum},'model',2,'varnames',{'FCtype','RemovedLUnum'},'display','on');
    % mean¡ÀSEM
    MeanChangedValue_MemtoMem = cellfun(@mean,ChangedValue_MemtoMem);
    SemChangedValue_MemtoMem = cellfun(@(x) std(x)/sqrt(numel(x)),ChangedValue_MemtoMem);
    MeanChangedValue_NonmemtoMem = cellfun(@mean,ChangedValue_NonmemtoMem);
    SemChangedValue_NonmemtoMem = cellfun(@(x) std(x)/sqrt(numel(x)),ChangedValue_NonmemtoMem);
    for i=1:numel(MeanChangedValue_MemtoMem)
        bar(1.1+1.8*(i-1),MeanChangedValue_MemtoMem(i),0.6,'facecolor',[0 0 0],'edgecolor',[0 0 0]); hold on
        bar(1.9+1.8*(i-1),MeanChangedValue_NonmemtoMem(i),0.6,'facecolor','none','edgecolor',[0 0 0]); hold on
        errorbar(1.1+1.8*(i-1),MeanChangedValue_MemtoMem(i),SemChangedValue_MemtoMem(i),'k','marker','none'); hold on
        errorbar(1.9+1.8*(i-1),MeanChangedValue_NonmemtoMem(i),SemChangedValue_NonmemtoMem(i),'k','marker','none'); hold on
    end
    set(gca,'XTick',1.5:1.8:2.5+1.8*(numel(MeanChangedValue_MemtoMem)-1),'XLim',[0.5 2.5+1.8*(numel(MeanChangedValue_MemtoMem)-1)],'YTick',-1:0.05:0,'YTickLabel',num2cell(-100:5:0),'YLim',[min(data) 0]);
    set(gcf,'Render','Painter'); saveas(gcf,sprintf('Decreased FU selectivity with removing FCSP-related spikes-FRrange-%d-%d-%s-CtrlGroup',FRrange(1),FRrange(2),Reg),'fig'); close all;
end
save(sprintf('Decreased FU selectivity with removing FCSP-related spikes-FRrange-%d-%dHz-%s-%s',FRrange(1),FRrange(2),Reg,Group),'ChangedValue_MemtoMem','ChangedValue_NonmemtoMem','UniqFC_MemtoMem','UniqFC_NonmemtoMem','UniqFR_MemtoMem','UniqFR_NonmemtoMem','tbl','-v7.3');

function [UniqFC,UniqFR] = GetUniqueFC(FCpair,FR)

UniqFC = cell(1,size(FCpair,2));
UniqFR = cell(1,size(FCpair,2));
for i = 1:size(FCpair,2)
    tempFC = FCpair(:,i);
    tempFC = vertcat(tempFC{:});
    tempFR = FR(:,i);
    tempFR = vertcat(tempFR{:});
    for j = 1:size(tempFC,1)
        CurrFC = tempFC{j};
        CurrFR = tempFR(j);
        RestFC = tempFC(setdiff(1:size(tempFC,1),j),:);
        SameFcId = [];
        for k = 1:size(RestFC,1)
            if isequal(CurrFC,RestFC{k})
                SameFcId = [SameFcId k];
            end
        end
        if isempty(SameFcId)
            UniqFC{i} = [UniqFC{i}; {CurrFC}];
            UniqFR{i} = [UniqFR{i}; CurrFR];
        end
    end
end
end

function AllFuSelec = GetSelecOfFuWithRemovingFC(UniqFCpair,mPFCUnitsID,mPFCUnitsSpikeTime,mPFCUnitsTrialMark,aAICUnitsID,aAICUnitsSpikeTime,aAICUnitsTrialMark)

AllFuSelec = cell(size(UniqFCpair));
for i = 1:numel(UniqFCpair) % different numbers of leading neurons
    for j = 1:size(UniqFCpair{i},1)
        fprintf('Process %dth one of total %d FC pairs with %d leading neurons\n',j,size(UniqFCpair{i},1),i);
        ID_Bin = UniqFCpair{i}{j}(1,1);
        ID_LU = UniqFCpair{i}{j}(:,2);
        ID_FU = UniqFCpair{i}{j}(1,3);
        
        %% Region of leading and following neurons
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
        [tempSelec_wiFC,tempSelec_woFC] = DecrFuSelecWithRemovingFCSP(ID_Bin,IDinReg_LU,IDinReg_FU,SpikeInfo_LU,SpikeInfo_FU,TrialMark);
        AllFuSelec{i} = [AllFuSelec{i}; horzcat(tempSelec_wiFC,tempSelec_woFC)];
    end
end
end

function [Selec_wiFC,Selec_woFC] = DecrFuSelecWithRemovingFCSP(BinID,IDinReg_PreUnit,IDinReg_PostUnit,RG_PreUnit,RG_PostUnit,TrialMark)

BaseLen = 4;
OdorLen = 1;
FR_wiFC_Go = [];
FR_woFC_Go = [];
FR_wiFC_NoGo = [];
FR_woFC_NoGo = [];

%% Spike rasters of target neurons
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

%% Spike rasters in Go and NoGo trials
RG_PreUnit_Go = RG_PreUnit(:,TrialID_Go);
RG_PreUnit_NoGo = RG_PreUnit(:,TrialID_NoGo);
RG_PostUnit_Go = RG_PostUnit(:,TrialID_Go);
RG_PostUnit_NoGo = RG_PostUnit(:,TrialID_NoGo);
for iTrialType = 1:2
    % remove FCSP events from spikes of following neurons in Go (iTrialType==1) and NoGo (iTrialType==2) trials
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
        % removing FCSP events of following neurons
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
        elseif iTrialType == 2
            FR_wiFC_NoGo = [FR_wiFC_NoGo; numel(tempRG_PostUnit)];
            FR_woFC_NoGo = [FR_woFC_NoGo; numel(tempRG_PostUnit(KeptSpkId))];
        end
    end
end
Selec_wiFC = (mean(FR_wiFC_Go)-mean(FR_wiFC_NoGo))/(mean(FR_wiFC_Go)+mean(FR_wiFC_NoGo));
Selec_woFC = (mean(FR_woFC_Go)-mean(FR_woFC_NoGo))/(mean(FR_woFC_Go)+mean(FR_woFC_NoGo));
end

% function AllFuAuc = GetAUCofallFUsAfterRemovingFC(FCpair,mPFCUnitsID,mPFCUnitsSpikeTime,mPFCUnitsTrialMark,aAICUnitsID,aAICUnitsSpikeTime,aAICUnitsTrialMark)
% TarFCpair = cell(size(FCpair));
% AllFuAuc = cell(size(FCpair));
% for i = 1:numel(FCpair) % different numbers of leading neurons
%     for j = 1:size(FCpair{i},1) % based on same following neurons
%         fprintf('Process %dth one of total %d FC pairs with %d leading neurons\n',j,size(FCpair{i},1),i);
%         temp = nchoosek(1:size(FCpair{i}{j},1),i);
%         for k = 1:size(temp,1)
%             TarFCpair{i}{end+1,1} = FCpair{i}{j}(temp(k,:),:);
%             ID_Bin = TarFCpair{i}{end}(1,1);
%             ID_LU = TarFCpair{i}{end}(:,2);
%             ID_FU = TarFCpair{i}{end}(1,3);
%             % judge region ID of leading and following neurons
%             IDinReg_LU = [];
%             if ismember(ID_LU(1,1),mPFCUnitsID) % leading neuron
%                 SpikeInfo_LU = mPFCUnitsSpikeTime;
%                 for iUnit = 1:numel(ID_LU)
%                     IDinReg_LU = [IDinReg_LU; find(mPFCUnitsID==ID_LU(iUnit))];
%                 end
%             else
%                 SpikeInfo_LU = aAICUnitsSpikeTime;
%                 for iUnit = 1:numel(ID_LU)
%                     IDinReg_LU = [IDinReg_LU; find(aAICUnitsID==ID_LU(iUnit))];
%                 end
%             end
%             if ismember(ID_FU,mPFCUnitsID) % following neuron
%                 SpikeInfo_FU = mPFCUnitsSpikeTime;
%                 IDinReg_FU = find(mPFCUnitsID==ID_FU);
%                 TrialMark = mPFCUnitsTrialMark{IDinReg_FU};
%             else
%                 SpikeInfo_FU = aAICUnitsSpikeTime;
%                 IDinReg_FU = find(aAICUnitsID==ID_FU);
%                 TrialMark = aAICUnitsTrialMark{IDinReg_FU};
%             end
%             [tempAUC_wiFC,tempAUC_woFC] = DecreaFUauROCofRemoveFCSP(ID_Bin,IDinReg_LU,IDinReg_FU,SpikeInfo_LU,SpikeInfo_FU,TrialMark);
%             AllFuAuc{i} = [AllFuAuc{i}; horzcat(tempAUC_wiFC,tempAUC_woFC)];
%         end
%     end
% end
% end
% 
% function [AUC_wiFC,AUC_woFC] = DecreaFUauROCofRemoveFCSP(BinID,IDinReg_PreUnit,IDinReg_PostUnit,RG_PreUnit,RG_PostUnit,TrialMark)
% 
% BaseLen = 4;
% OdorLen = 1;
% FR_wiFC_Go = [];
% FR_woFC_Go = [];
% FR_wiFC_NoGo = [];
% FR_woFC_NoGo = [];
% 
% % raster timestamps of target units
% temp = RG_PreUnit(IDinReg_PreUnit,:);
% RG_PreUnit = [];
% for i = 1:numel(temp)
%     RG_PreUnit = [RG_PreUnit; temp{i}];
% end
% temp = RG_PostUnit(IDinReg_PostUnit,:);
% RG_PostUnit = [];
% for i = 1:numel(temp)
%     RG_PostUnit = [RG_PostUnit; temp{i}];
% end
% 
% % trial ID of Go and NoGo trials
% TrialID_Go = find(TrialMark(:,2)==1);
% TrialID_NoGo = find(TrialMark(:,2)==2);
% % raster timestamps of target units in Go and NoGo trials
% RG_PreUnit_Go = RG_PreUnit(:,TrialID_Go);
% RG_PreUnit_NoGo = RG_PreUnit(:,TrialID_NoGo);
% RG_PostUnit_Go = RG_PostUnit(:,TrialID_Go);
% RG_PostUnit_NoGo = RG_PostUnit(:,TrialID_NoGo);
% for iTrialType = 1:2
%     % removing FCSP-related spikes from activities of following neurons in Go (iTrialType==1) and NoGo (iTrialType==2) trials
%     if iTrialType == 1
%         RG_PostUnit = RG_PostUnit_Go;
%         RG_PreUnit = RG_PreUnit_Go;
%     elseif iTrialType == 2
%         RG_PostUnit = RG_PostUnit_NoGo;
%         RG_PreUnit = RG_PreUnit_NoGo;
%     end
%     for iTrial = 1:numel(RG_PostUnit)
%         % raster timestamp of postunit
%         tempRG_PostUnit = RG_PostUnit(:,iTrial);
%         tempRG_PostUnit = vertcat(tempRG_PostUnit{:});
%         tempRG_PostUnit = tempRG_PostUnit(tempRG_PostUnit>BaseLen+OdorLen+BinID-1 & tempRG_PostUnit<=BaseLen+OdorLen+BinID);
%         % raster timestamp of preunit
%         tempRG_PreUnit = RG_PreUnit(:,iTrial);
%         tempRG_PreUnit = vertcat(tempRG_PreUnit{:});
%         tempRG_PreUnit = tempRG_PreUnit(tempRG_PreUnit>=BaseLen+OdorLen+BinID-1 & tempRG_PreUnit<BaseLen+OdorLen+BinID);
%         % removing FCSP-related spikes of following neurons
%         KeptSpkId = [];
%         for iSpike = 1:numel(tempRG_PostUnit)
%             temp = find(tempRG_PostUnit(iSpike)-tempRG_PreUnit<=0.01 & tempRG_PostUnit(iSpike)-tempRG_PreUnit>0.002);
%             if isempty(temp)
%                 KeptSpkId = [KeptSpkId; iSpike];
%             end
%         end
%         if iTrialType == 1
%             FR_wiFC_Go = [FR_wiFC_Go; numel(tempRG_PostUnit)];
%             FR_woFC_Go = [FR_woFC_Go; numel(tempRG_PostUnit(KeptSpkId))];
%         elseif iTrialType == 2
%             FR_wiFC_NoGo = [FR_wiFC_NoGo; numel(tempRG_PostUnit)];
%             FR_woFC_NoGo = [FR_woFC_NoGo; numel(tempRG_PostUnit(KeptSpkId))];
%         end
%     end
% end
% AUC_wiFC = CalculateSingleUnitAUCinOneBin(FR_wiFC_Go,FR_wiFC_NoGo);
% AUC_woFC = CalculateSingleUnitAUCinOneBin(FR_woFC_Go,FR_woFC_NoGo);
% end
% 
% function SingleUnitAUCvalue = CalculateSingleUnitAUCinOneBin(FR_S1,FR_S2)
% 
% % calculate decision variable (DV)
% [DV_S1,DV_S2] = CalculateDV(FR_S1,FR_S2);
% % FPR and TPR
% [TPR,FPR] = CalculateTPRFPRValue(DV_S1,DV_S2);
% % calculate auROC
% SingleUnitAUCvalue = CalculateAUCValue(FPR,TPR);
% end
