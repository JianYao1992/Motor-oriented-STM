%% Example spike rasters showing following memory neurons in FC with STM-encoding leading neurons

clear; clc; close all;

%% Assignment
Type = 'Leading memory unit';
if strcmp(Type,'Leading memory unit')
    Pref = 1;
elseif strcmp(Type,'Leading non-memory unit')
    Pref = 0;
end
ShownTriNum_hit = 3;
IsConsideringCRTrial = 0;
if IsConsideringCRTrial == 1
    ShownTriNum_CR = 3;
elseif IsConsideringCRTrial == 0
    ShownTriNum_CR = 0;
end
% event duration
BaseLen = 1;
BaseLen_SpikeRaster = 4;
SampOdorLen = 1;
DelayLen = 4;
DelayPeriod = BaseLen_SpikeRaster+SampOdorLen+1:BaseLen_SpikeRaster+SampOdorLen+DelayLen;
C1 = [0.5 0.5 0.5];
C2 = [1 0 0];
C3 = [0 0 1];

%% Enter path
CurrPath = uigetdir;
cd(CurrPath);

%% ID of mPFC and aAIC neurons
load('UnitsSummary_CtrlGroup.mat');
% mPFC
ID_mPFC = [];
for iUnit = 1:size(UnitsSummary.mPFC,1) 
    ID_mPFC = [ID_mPFC UnitsSummary.mPFC{iUnit,1}(1)];
end
% aAIC
ID_aAIC = [];
for iUnit = 1:size(UnitsSummary.aAIC,1) 
    ID_aAIC = [ID_aAIC UnitsSummary.aAIC{iUnit,1}(1)];
end

%% Sort information of neurons
load('NeuronOriginandSpikeRasterInformationforCtrlGroup.mat');
ID = [ID_mPFC'; ID_aAIC'];
ID = num2cell(ID);
AllUnitsInfo = [UnitsSummary.mPFC; UnitsSummary.aAIC];
AllUnitsInfo = [ID AllUnitsInfo];
SortedAllUnitsInfo = sortrows(AllUnitsInfo,1);
UnitTrialMark = horzcat(ID,[mPFCUnitTrialMark; aAICUnitTrialMark]);
SortedUnitTrialMark  = sortrows(UnitTrialMark,1);
SpikeTimeIndividualTrial = horzcat(ID,[mPFCSpikeTimeinIndividualTrial; aAICSpikeTimeinIndividualTrial]);
SortedSpikeTimeIndividualTrial = sortrows(SpikeTimeIndividualTrial,1);
clear AllUnitsInfo UnitTrialMark SpikeTimeIndividualTrial

%% Load result of FC
load('FcTemporalPattern_Local+Cross-region_all_CtrlGroup_hit trials.mat');
FCinEachBin_hit = FCinEachBin;
load('FcTemporalPattern_Local+Cross-region_all_CtrlGroup_CR trials.mat');
FCinEachBin_CR = FCinEachBin;

%% Selectivity patterns of neuronal pairs showing FC in hit trials
SelecPattern = [];
for iFC = 1:size(FCinEachBin_hit,1)
    tempUnit1 = FCinEachBin_hit(iFC,1);
    tempUnit2 = FCinEachBin_hit(iFC,2);
    tempSelecPattern_Unit1 = SortedAllUnitsInfo{tempUnit1,4};
    tempSelecPattern_Unit2 = SortedAllUnitsInfo{tempUnit2,4};
    SelecPattern = [SelecPattern; horzcat(tempSelecPattern_Unit1,tempSelecPattern_Unit2)];
end

%% Find neuronal pairs reaching criteria
for iFC = 1:size(FCinEachBin_hit,1)
    fprintf('Analyze %dth FC of total %d FC\n',iFC,size(FCinEachBin_hit,1));
    tempUnit1 = FCinEachBin_hit(iFC,1);
    tempUnit2 = FCinEachBin_hit(iFC,2);
    if ismember(tempUnit1,ID_mPFC)
        IDinReg_unit1 = find(ID_mPFC == tempUnit1);
        Reg_unit1 = 'mPFC';
    elseif ismember(tempUnit1,ID_aAIC)
        IDinReg_unit1 = find(ID_aAIC == tempUnit1);
        Reg_unit1 = 'aAIC';
    end
    if ismember(tempUnit2,ID_mPFC)
        IDinReg_unit2 = find(ID_mPFC == tempUnit2);
        Reg_unit2 = 'mPFC';
    elseif ismember(tempUnit2,ID_aAIC)
        IDinReg_unit2 = find(ID_aAIC == tempUnit2);
        Reg_unit2 = 'aAIC';
    end
    tempFCpattern = FCinEachBin_hit(iFC,3:6);
    tempFCperiod = find(tempFCpattern>0);
    for iTime = tempFCperiod
        if SelecPattern(iFC,iTime+2) == Pref && SelecPattern(iFC,iTime+9) == Pref
            if IsConsideringCRTrial == 1
                IsFCinCRTrial = find(FCinEachBin_CR(:,1) == tempUnit1 & FCinEachBin_CR(:,2) == tempUnit2);
            elseif IsConsideringCRTrial == 0
                IsFCinCRTrial = [];
            end
            if isempty(IsFCinCRTrial) || (~isempty(IsFCinCRTrial) && FCinEachBin_CR(1,2+iTime) == 0)
                tempTrialMark = SortedUnitTrialMark{tempUnit1,2};
                tempTrialID_hit = find(tempTrialMark(:,4) == 1);
                tempTrialID_CR = find(tempTrialMark(:,4) == 4);
                tempSpikeTime_Unit1 = SortedSpikeTimeIndividualTrial{tempUnit1,2};
                tempSpikeTime_Unit2 = SortedSpikeTimeIndividualTrial{tempUnit2,2};
                % spike time of two units in hit trials
                tempSpikeTime_Unit1_hit = tempSpikeTime_Unit1(:,tempTrialID_hit);
                tempSpikeTime_Unit2_hit = tempSpikeTime_Unit2(:,tempTrialID_hit);
                % spike time of two units in CR trials
                tempSpikeTime_Unit1_CR = tempSpikeTime_Unit1(:,tempTrialID_CR);
                tempSpikeTime_Unit2_CR = tempSpikeTime_Unit2(:,tempTrialID_CR);
                % number of trial-based FCSP events in hit trials
                FcspNum_hit = [];
                for itr = 1:length(tempSpikeTime_Unit1_hit)
                    FCspiketime = [];
                    for iSpike = 1:length(tempSpikeTime_Unit1_hit{itr})
                        if ~isempty(find(tempSpikeTime_Unit2_hit{itr}(:,1) - tempSpikeTime_Unit1_hit{itr}(iSpike) <= 0.01 & tempSpikeTime_Unit2_hit{itr}(:,1) - tempSpikeTime_Unit1_hit{itr}(iSpike) > 0.002))
                            FCspiketime = [FCspiketime tempSpikeTime_Unit1_hit{itr}(iSpike)];
                        end
                    end
                    FCspiketime(FCspiketime < BaseLen_SpikeRaster+SampOdorLen+iTime-1 | FCspiketime >= BaseLen_SpikeRaster+SampOdorLen+iTime) = [];
                    FcspNum_hit = [FcspNum_hit; numel(FCspiketime)];
                end
                % number of trial-based FCSP events in CR trials
                FcspNum_CR = [];
                for itr = 1:length(tempSpikeTime_Unit1_CR)
                    FCspiketime = [];
                    for iSpike = 1:length(tempSpikeTime_Unit1_CR{itr})
                        if ~isempty(find(tempSpikeTime_Unit2_CR{itr}(:,1) - tempSpikeTime_Unit1_CR{itr}(iSpike) <= 0.01 & tempSpikeTime_Unit2_CR{itr}(:,1) - tempSpikeTime_Unit1_CR{itr}(iSpike) > 0.002))
                            FCspiketime = [FCspiketime tempSpikeTime_Unit1_CR{itr}(iSpike)];
                        end
                    end
                    FCspiketime(FCspiketime < BaseLen_SpikeRaster+SampOdorLen+iTime-1 | FCspiketime >= BaseLen_SpikeRaster+SampOdorLen+iTime) = [];
                    FcspNum_CR = [FcspNum_CR; numel(FCspiketime)];
                end
                if nnz(FcspNum_hit>=2) >= ShownTriNum_hit && nnz(FcspNum_CR==0) >= ShownTriNum_CR
                 %% Raster plot
                    % ID of hit trials
                    [~,sortid] = sortrows(FcspNum_hit,'descend');
                    SelectedTriID_hit = sortid(1:ShownTriNum_hit);
                    % ID of CR trials
                    id = find(FcspNum_CR == 0);
                    id = id(randperm(numel(id)));
                    SelectedTriID_CR = id(1:ShownTriNum_CR);
                    % spike raster
                    figure('Color','w','Position',[500,500,700,300])
                    subplot(1,2,1); % hit trial
                    for itr = 1:ShownTriNum_hit
                        for iSpike = 1:length(tempSpikeTime_Unit1_hit{SelectedTriID_hit(itr)}) % leading neuron
                            if ~isempty(find(tempSpikeTime_Unit2_hit{SelectedTriID_hit(itr)}(:,1) - tempSpikeTime_Unit1_hit{SelectedTriID_hit(itr)}(iSpike) <= 0.01 & tempSpikeTime_Unit2_hit{SelectedTriID_hit(itr)}(:,1) - tempSpikeTime_Unit1_hit{SelectedTriID_hit(itr)}(iSpike) > 0.002))
                                plot([tempSpikeTime_Unit1_hit{SelectedTriID_hit(itr)}(iSpike) tempSpikeTime_Unit1_hit{SelectedTriID_hit(itr)}(iSpike)], [1.5+3*(ShownTriNum_hit-itr) 2.5+3*(ShownTriNum_hit-itr)],'color',C2);
                            else
                                plot([tempSpikeTime_Unit1_hit{SelectedTriID_hit(itr)}(iSpike) tempSpikeTime_Unit1_hit{SelectedTriID_hit(itr)}(iSpike)], [1.5+3*(ShownTriNum_hit-itr) 2.5+3*(ShownTriNum_hit-itr)],'color',C1);
                            end
                            hold on
                        end
                        for iSpike = 1:length(tempSpikeTime_Unit2_hit{SelectedTriID_hit(itr)}) % following neuron
                            if ~isempty(find(tempSpikeTime_Unit2_hit{SelectedTriID_hit(itr)}(iSpike) - tempSpikeTime_Unit1_hit{SelectedTriID_hit(itr)}(:,1) <= 0.01 & tempSpikeTime_Unit2_hit{SelectedTriID_hit(itr)}(iSpike) - tempSpikeTime_Unit1_hit{SelectedTriID_hit(itr)}(:,1) > 0.002))
                                plot([tempSpikeTime_Unit2_hit{SelectedTriID_hit(itr)}(iSpike) tempSpikeTime_Unit2_hit{SelectedTriID_hit(itr)}(iSpike)], [0.5+3*(ShownTriNum_hit-itr) 1.5+3*(ShownTriNum_hit-itr)],'color',C3);
                            else
                                plot([tempSpikeTime_Unit2_hit{SelectedTriID_hit(itr)}(iSpike) tempSpikeTime_Unit2_hit{SelectedTriID_hit(itr)}(iSpike)], [0.5+3*(ShownTriNum_hit-itr) 1.5+3*(ShownTriNum_hit-itr)],'color',C1);
                            end
                            hold on
                        end
                    end
                    set(gca,'XTick',BaseLen_SpikeRaster+SampOdorLen:1:BaseLen_SpikeRaster+SampOdorLen+DelayLen,'XTickLabel',{'','','','',''},'YTick',zeros(1,0),'FontName','Arial','FontSize',16,'XLim',[BaseLen_SpikeRaster+SampOdorLen+iTime-1 BaseLen_SpikeRaster+SampOdorLen+iTime],'YLim',[0.25 2.75+3*(ShownTriNum_hit-1)]);
                    if IsConsideringCRTrial == 1
                        subplot(1,2,2); % CR trial
                        for itr = 1:ShownTriNum_CR
                            for iSpike = 1:length(tempSpikeTime_Unit1_CR{SelectedTriID_CR(itr)}) % leading unit
                                if ~isempty(find(tempSpikeTime_Unit2_CR{SelectedTriID_CR(itr)}(:,1) - tempSpikeTime_Unit1_CR{SelectedTriID_CR(itr)}(iSpike) <= 0.01 & tempSpikeTime_Unit2_CR{SelectedTriID_CR(itr)}(:,1) - tempSpikeTime_Unit1_CR{SelectedTriID_CR(itr)}(iSpike) > 0.002))
                                    plot([tempSpikeTime_Unit1_CR{SelectedTriID_CR(itr)}(iSpike) tempSpikeTime_Unit1_CR{SelectedTriID_CR(itr)}(iSpike)], [1.5+3*(ShownTriNum_CR-itr) 2.5+3*(ShownTriNum_CR-itr)],'color',C2);
                                else
                                    plot([tempSpikeTime_Unit1_CR{SelectedTriID_CR(itr)}(iSpike) tempSpikeTime_Unit1_CR{SelectedTriID_CR(itr)}(iSpike)], [1.5+3*(ShownTriNum_CR-itr) 2.5+3*(ShownTriNum_CR-itr)],'color',C1);
                                end
                                hold on
                            end
                            for iSpike = 1:length(tempSpikeTime_Unit2_CR{SelectedTriID_CR(itr)}) % following unit
                                if ~isempty(find(tempSpikeTime_Unit2_CR{SelectedTriID_CR(itr)}(iSpike) - tempSpikeTime_Unit1_CR{SelectedTriID_CR(itr)}(:,1) <= 0.01 & tempSpikeTime_Unit2_CR{SelectedTriID_CR(itr)}(iSpike) - tempSpikeTime_Unit1_CR{SelectedTriID_CR(itr)}(:,1) > 0.002))
                                    plot([tempSpikeTime_Unit2_CR{SelectedTriID_CR(itr)}(iSpike) tempSpikeTime_Unit2_CR{SelectedTriID_CR(itr)}(iSpike)], [0.5+3*(ShownTriNum_CR-itr) 1.5+3*(ShownTriNum_CR-itr)],'color',C3);
                                else
                                    plot([tempSpikeTime_Unit2_CR{SelectedTriID_CR(itr)}(iSpike) tempSpikeTime_Unit2_CR{SelectedTriID_CR(itr)}(iSpike)], [0.5+3*(ShownTriNum_CR-itr) 1.5+3*(ShownTriNum_CR-itr)],'color',C1);
                                end
                                hold on
                            end
                        end
                        set(gca,'XTick',BaseLen_SpikeRaster+SampOdorLen:1:BaseLen_SpikeRaster+SampOdorLen+DelayLen,'XTickLabel',{'','','','',''},'YTick',zeros(1,0),'FontName','Arial','FontSize',16,'XLim',[BaseLen_SpikeRaster+SampOdorLen+iTime-1 BaseLen_SpikeRaster+SampOdorLen+iTime],'YLim',[0.25 2.75+3*(ShownTriNum_CR-1)]);
                    end
                    set(gcf,'Render','Painter'); saveas(gcf,sprintf('%s-FCSP in %s#%d-%s#%d_at delay %dth s',Type,Reg_unit1,IDinReg_unit1,Reg_unit2,IDinReg_unit2,iTime),'fig'); close all;
                end
            end
        end
    end
end









