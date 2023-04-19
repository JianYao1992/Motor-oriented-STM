%% Rasters of neurons in a FC (transient) pair over the time course of delay

clear; clc; close all;

%% Assignment
BaseLen = 4;
OdorLen = 1;
DelayLen = 4;
DelayPeriod = BaseLen+OdorLen:1:BaseLen+OdorLen+DelayLen;
PairType = 'Con Act';
ShownTrialNum = 3;
C1 = [0.5 0.5 0.5];
C2 = [220 39 114]/255;
C3 = [39 148 211]/255;

%% Enter path
CurrPath = uigetdir;
cd(CurrPath);

%% Load result of FC in hit trials
load('FcTemporalPattern_Local+Cross-region_Within memory neurons_CtrlGroup_Hit trials.mat');

%% ID of mPFC and aAIC neurons
load('UnitsSummary_CtrlGroup.mat');
ID_mPFC = [];
ID_aAIC = [];
% mPFC
for iUnit = 1:size(UnitsSummary.mPFC,1)
    ID_mPFC = [ID_mPFC UnitsSummary.mPFC{iUnit,1}(1)];
end
% aAIC
for iUnit = 1:size(UnitsSummary.aAIC,1)
    ID_aAIC = [ID_aAIC UnitsSummary.aAIC{iUnit,1}(1)];
end

%% Load result of spike rasters
load('NeuronOriginandSpikeRasterInformationforCtrlGroup.mat');

%% Integrate result of mPFC and aAIC
% information of neurons
ID = [ID_mPFC'; ID_aAIC'];
ID = num2cell(ID);
AllUnitsInfo = [UnitsSummary.mPFC; UnitsSummary.aAIC];
AllUnitsInfo = [ID AllUnitsInfo];
SortedAllUnitsInfo = sortrows(AllUnitsInfo,1);
% trial marker
UnitTrialMark = horzcat(ID,[mPFCUnitTrialMark;aAICUnitTrialMark]);
SortedUnitTrialMark = sortrows(UnitTrialMark,1);
% spike raster
SpikeTimeIndividualTrial = horzcat(ID,[mPFCSpikeTimeinIndividualTrial;aAICSpikeTimeinIndividualTrial]);
SortedSpikeTimeIndividualTrial = sortrows(SpikeTimeIndividualTrial,1);
clear UnitsSummary AllUnitsInfo UnitTrialMark SpikeTimeIndividualTrial

%% Find example case
PairCase = cell(0);
for iPair = 1:size(FCinEachBin,1)
    disp([num2str(iPair) ' th pair of total' num2str(size(CouplinginEachBin,1)) 'pairs']);
    Unit1 = FCinEachBin(iPair,1);
    Unit2 = FCinEachBin(iPair,2);
    if ismember(Unit1,ID_mPFC)
        Reg_Unit1 = 'mPFC';
        IDinReg_Unit1 = find(ID_mPFC==Unit1);
    else
        Reg_Unit1 = 'aAIC';
        IDinReg_Unit1 = find(ID_aAIC==Unit1);
    end
    if ismember(Unit2,ID_mPFC)
        Reg_Unit2 = 'mPFC';
        IDinReg_Unit2 = find(ID_mPFC==Unit2);
    else
        Reg_Unit2 = 'aAIC';
        IDinReg_Unit2 = find(ID_aAIC==Unit2);
    end
    CouplingPeriod = 2 + find(FCinEachBin(iPair,3:6)>0);
    for iTime = CouplingPeriod
        switch PairType
            case 'Con Act'
                IsFCinTarPair = (SortedAllUnitsInfo{Unit1,4}(iTime)*SortedAllUnitsInfo{Unit2,4}(iTime)==1 | SortedAllUnitsInfo{Unit1,4}(iTime)*SortedAllUnitsInfo{Unit2,4}(iTime)==4)...
                    & ((any(SortedAllUnitInfo{Unit1,4}(:,3:6)==1,2)&any(SortedAllUnitInfo{Unit1,4}(:,3:6)==2,2)) | (any(SortedAllUnitInfo{Unit2,4}(:,3:6)==1,2)&any(SortedAllUnitInfo{Unit2,4}(:,3:6)==2,2)));
            case 'Con Inact'
                IsFCinTarPair = (SortedAllUnitsInfo{Unit1,4}(iTime)*SortedAllUnitsInfo{Unit2,4}(iTime)==0 & max(SortedAllUnitsInfo{Unit1,4}(3:6)) == max(SortedAllUnitsInfo{Unit2,4}(3:6)))...
                    & ((any(SortedAllUnitInfo{Unit1,4}(:,3:6)==1,2)&any(SortedAllUnitInfo{Unit1,4}(:,3:6)==2,2)) | (any(SortedAllUnitInfo{Unit2,4}(:,3:6)==1,2)&any(SortedAllUnitInfo{Unit2,4}(:,3:6)==2,2)));
            case 'Incon Act'
                IsFCinTarPair = (SortedAllUnitsInfo{Unit1,4}(iTime)*SortedAllUnitsInfo{Unit2,4}(iTime)==2) & ((any(SortedAllUnitInfo{Unit1,4}(:,3:6)==1,2)&any(SortedAllUnitInfo{Unit1,4}(:,3:6)==2,2)) | (any(SortedAllUnitInfo{Unit2,4}(:,3:6)==1,2)&any(SortedAllUnitInfo{Unit2,4}(:,3:6)==2,2)));
            case 'Incon Inact'
                IsFCinTarPair = (SortedAllUnitsInfo{Unit1,4}(iTime)*SortedAllUnitsInfo{Unit2,4}(iTime)==0 & max(SortedAllUnitsInfo{Unit1,4}(3:6)) ~= max(SortedAllUnitsInfo{Unit2,4}(3:6)))...
                    & ((any(SortedAllUnitInfo{Unit1,4}(:,3:6)==1,2)&any(SortedAllUnitInfo{Unit1,4}(:,3:6)==2,2)) | (any(SortedAllUnitInfo{Unit2,4}(:,3:6)==1,2)&any(SortedAllUnitInfo{Unit2,4}(:,3:6)==2,2)));
        end
        if IsFCinTarPair == 1 % belong to one type of neuronal pair
            %% Number of FCSP events in example trials
            CouplingPeriod = BaseLen + OdorLen + find(FCinEachBin(iPair,3:6)>0);
            TrialMark = SortedUnitTrialMark{Unit1,2};
            TrialID_hit = find(TrialMark(:,4) == 1);
            SpikeTime_Unit1 = SortedSpikeTimeIndividualTrial{Unit1,2};
            SpikeTime_Unit2 = SortedSpikeTimeIndividualTrial{Unit2,2};
            SpikeTime_Unit1 = SpikeTime_Unit1(:,TrialID_hit);
            SpikeTime_Unit2 = SpikeTime_Unit2(:,TrialID_hit);
            ExampleTrialID = [];
            FcspNum = [];
            for itr = 1:length(SpikeTime_Unit1) % hit trial
                FcSpikeTime = [];
                for iSpike = 1:length(SpikeTime_Unit1{itr})
                    if ~isempty(find(SpikeTime_Unit2{itr}(:,1) - SpikeTime_Unit1{itr}(iSpike) <= 0.01 & SpikeTime_Unit2{itr}(:,1) - SpikeTime_Unit1{itr}(iSpike) > 0.002))
                        FcSpikeTime = [FcSpikeTime SpikeTime_Unit1{itr}(iSpike)];
                    end
                end
                FcSpikeTime(FcSpikeTime < BaseLen+OdorLen | FcSpikeTime >= BaseLen+OdorLen+DelayLen) = [];
                if ~isempty(FcSpikeTime)
                    NonFcPeriod = setdiff(DelayPeriod,CouplingPeriod);
                    tempCouplingNum = [];
                    for iBin = NonFcPeriod
                        tempCouplingNum = [tempCouplingNum length(find(FcSpikeTime >= iBin -1 & FcSpikeTime < iBin))];
                    end
                    if all(tempCouplingNum == 0)
                        tempCouplingNum = [];
                        for iBin = CouplingPeriod
                            tempCouplingNum = [tempCouplingNum length(find(FcSpikeTime >= iBin -1 & FcSpikeTime < iBin))];
                        end
                        if all(tempCouplingNum >= 1)
                            ExampleTrialID = [ExampleTrialID itr];
                            FcspNum = [FcspNum mean(tempCouplingNum)];
                        end
                    end
                end
            end
            if numel(FcspNum) >= 3
                [~,sortid] = sortrows(FcspNum(:),'descend');
                ExampleTrialID = ExampleTrialID(sortid(1:ShownTrialNum));
                
                %% Plot spike rasters of example trials
                figure('Color','w','Position',[500,500,700,250])
                for itr = 1:ShownTrialNum
                    % unit 1
                    for iSpike = 1:length(SpikeTime_Unit1{ExampleTrialID(itr)})
                        if ~isempty(find(SpikeTime_Unit2{ExampleTrialID(itr)}(:,1) - SpikeTime_Unit1{ExampleTrialID(itr)}(iSpike) <= 0.01 & SpikeTime_Unit2{ExampleTrialID(itr)}(:,1) - SpikeTime_Unit1{ExampleTrialID(itr)}(iSpike) > 0.002))
                            plot([SpikeTime_Unit1{ExampleTrialID(itr)}(iSpike) SpikeTime_Unit1{ExampleTrialID(itr)}(iSpike)], [1.5+3*(ShownTrialNum-itr) 2.5+3*(ShownTrialNum-itr)],'color',C2);
                        else
                            plot([SpikeTime_Unit1{ExampleTrialID(itr)}(iSpike) SpikeTime_Unit1{ExampleTrialID(itr)}(iSpike)], [1.5+3*(ShownTrialNum-itr) 2.5+3*(ShownTrialNum-itr)],'color',C1);
                        end
                        hold on
                    end
                    % unit 2
                    for iSpike = 1:length(SpikeTime_Unit2{ExampleTrialID(itr)})
                        if ~isempty(find(SpikeTime_Unit2{ExampleTrialID(itr)}(iSpike) - SpikeTime_Unit1{ExampleTrialID(itr)}(:,1) <= 0.01 & SpikeTime_Unit2{ExampleTrialID(itr)}(iSpike) - SpikeTime_Unit1{ExampleTrialID(itr)}(:,1) > 0.002))
                            plot([SpikeTime_Unit2{ExampleTrialID(itr)}(iSpike) SpikeTime_Unit2{ExampleTrialID(itr)}(iSpike)], [0.5+3*(ShownTrialNum-itr) 1.5+3*(ShownTrialNum-itr)],'color',C3);
                        else
                            plot([SpikeTime_Unit2{ExampleTrialID(itr)}(iSpike) SpikeTime_Unit2{ExampleTrialID(itr)}(iSpike)], [0.5+3*(ShownTrialNum-itr) 1.5+3*(ShownTrialNum-itr)],'color',C1);
                        end
                        hold on
                    end
                end
                set(gca,'XTick',BaseLen:1:BaseLen+OdorLen+DelayLen,'XTickLabel',num2cell(BaseLen:1:BaseLen+OdorLen+DelayLen),'YTick',zeros(1,0),'FontName','Arial','FontSize',16,'XLim',[BaseLen+OdorLen BaseLen+OdorLen+DelayLen],'YLim',[0.25 2.75+3*(ShownTrialNum-1)]);
                set(gcf,'Render','Painter'); saveas(gcf,sprintf('FCSP in %s#%d-%s#%d-at%dth second of delay-%s',Reg_Unit1,IDinReg_Unit1,Reg_Unit2,IDinReg_Unit2,iTime-2,PairType),'fig'); close all;
            end
        end
    end
end




