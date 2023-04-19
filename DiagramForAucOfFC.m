%% Example neuronal pairs with STM-encoding FCSP events

clear; clc; close all;

%% Assignment
Type = 'leading coding unit';
if strcmp(Type,'leading coding unit')
    SamplePrefer = 1;
elseif strcmp(Type,'leading non-coding unit')
    SamplePrefer = 0;
end
ShownTriNum_Go = 3;
IsConsiderNoGoTrial = 1;
if IsConsiderNoGoTrial == 1
    ShownTriNum_NoGo = 3;
elseif IsConsiderNoGoTrial == 0
    ShownTriNum_NoGo = 0;
end
% event duration
BaseLen = 1;
BaseLen_SpikeRaster = 4;
SampOdorLen = 1;
DelayLen = 4;
DelayPeriod = BaseLen_SpikeRaster+SampOdorLen+1:BaseLen_SpikeRaster+SampOdorLen+DelayLen;
C1 = [0.5 0.5 0.5];
C2 = [220 39 114]/255;
C3 = [39 148 211]/255;

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
UnitTrialMark = horzcat(ID,[mPFCUnitTrialMark;aAICUnitTrialMark]);
SortedUnitTrialMark  = sortrows(UnitTrialMark,1);
SpikeTimeIndivTrial = horzcat(ID,[mPFCSpikeTimeinIndividualTrial; aAICSpikeTimeinIndividualTrial]);
SortedSpikeTimeIndividualTrial = sortrows(SpikeTimeIndivTrial,1);
clear AllUnitsInfo UnitTrialMark SpikeTimeIndivTrial

%% Filter FC neuronal pairs, plot rastergram showing FCSP events
for iBin = 1:DelayLen
    load(['FCofNeurons_PropertyPopulation_' num2str(bin) '_' num2str(bin+1) '_GoNogo' Group '.mat']);

    %% Align leading neurons and following neurons from left to right
    for iPair = 1:size(conn_chain_go,1) % Go trials
        if conn_chain_go(iPair,3) < 0
            temp = conn_chain_go(iPair,1); conn_chain_go(iPair,1) = conn_chain_go(iPair,2); conn_chain_go(iPair,2) = temp;
            conn_chain_go(iPair,3) = -1*conn_chain_go(iPair,3);
            temp = pref_chain_go(iPair,1:7); pref_chain_go(iPair,1:7) = pref_chain_go(iPair,8:14); pref_chain_go(iPair,8:14) = temp;
            temp = reg_chain_go(iPair,1); reg_chain_go(iPair,1) = reg_chain_go(iPair,2); reg_chain_go(iPair,2) = temp;
            reg_chain_go(iPair,3) = -1*reg_chain_go(iPair,3);
        end
    end
    for iPair = 1:size(conn_chain_nogo,1) % NoGo trials
        if conn_chain_nogo(iPair,3) < 0
            temp = conn_chain_nogo(iPair,1); conn_chain_nogo(iPair,1) = conn_chain_nogo(iPair,2); conn_chain_nogo(iPair,2) = temp;
            conn_chain_nogo(iPair,3) = -1*conn_chain_nogo(iPair,3);
            temp = pref_chain_nogo(iPair,1:7); pref_chain_nogo(iPair,1:7) = pref_chain_nogo(iPair,8:14); pref_chain_nogo(iPair,8:14) = temp;
            temp = reg_chain_nogo(iPair,1); reg_chain_nogo(iPair,1) = reg_chain_nogo(iPair,2); reg_chain_nogo(iPair,2) = temp;
            reg_chain_nogo(iPair,3) = -1*reg_chain_nogo(iPair,3);
        end
    end
    
    %% Remove neuronal pairs (showing FC in Go trials) of non-memory neurons or switched-coding neurons, keep pairs with two neurons showing selectivity in the same bin
    RemainFCId = ~(all(pref_chain_go(:,3:6)==0,2) | all(pref_chain_go(:,10:13)==0,2) | ((any(pref_chain_go(:,3:6)==1,2) & any(pref_chain_go(:,3:6)==2,2)) |(any(pref_chain_go(:,10:13)==1,2) & any(pref_chain_go(:,10:13)==2,2)))) & pref_chain_go(:,2+iBin) == SamplePrefer & pref_chain_go(:,9+iBin) == SamplePrefer;
    conn_chain_go = conn_chain_go(RemainFCId,:);
    pref_chain_go = pref_chain_go(RemainFCId,:);
    reg_chain_go = reg_chain_go(RemainFCId,:);
    
    %% Region and neuron ID of neuronal pairs
    for iFC = 1:size(conn_chain_go,1)
        fprintf('Analyze %dth FC in total %d FC at bin %d of delay\n',iFC,size(conn_chain_go,1),iBin);
        tempUnit1 = conn_chain_go(iFC,1);
        tempUnit2 = conn_chain_go(iFC,2);
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
        if IsConsiderNoGoTrial == 1
            IsFCinNoGoTrial = find(conn_chain_nogo(:,1) == tempUnit1 & conn_chain_nogo(:,2) == tempUnit2);
        elseif IsConsiderNoGoTrial == 0
            IsFCinNoGoTrial = [];
        end
        if isempty(IsFCinNoGoTrial)
            tempTrialMark = SortedUnitTrialMark{tempUnit1,2};
            tempTrialID_Go = find(tempTrialMark(:,2)==1);
            tempTrialID_NoGo = find(tempTrialMark(:,2)==2);
            tempSpikeTime_Unit1 = SortedSpikeTimeIndividualTrial{tempUnit1,2};
            tempSpikeTime_Unit2 = SortedSpikeTimeIndividualTrial{tempUnit2,2};
            % spike time of two units in Go trials
            tempSpikeTime_Unit1_Go = tempSpikeTime_Unit1(:,tempTrialID_Go);
            tempSpikeTime_Unit2_Go = tempSpikeTime_Unit2(:,tempTrialID_Go);
            % spike time of two units in NoGo trials
            tempSpikeTime_Unit1_NoGo = tempSpikeTime_Unit1(:,tempTrialID_NoGo);
            tempSpikeTime_Unit2_NoGo = tempSpikeTime_Unit2(:,tempTrialID_NoGo);
            
            %% number of FCSP events of each Go trial
            FcspNum_Go = [];
            for itr = 1:length(tempSpikeTime_Unit1_Go)
                FcspStamp = [];
                for iSpike = 1:length(tempSpikeTime_Unit1_Go{itr})
                    if ~isempty(find(tempSpikeTime_Unit2_Go{itr}(:,1)-tempSpikeTime_Unit1_Go{itr}(iSpike)<=0.01 & tempSpikeTime_Unit2_Go{itr}(:,1)-tempSpikeTime_Unit1_Go{itr}(iSpike)>0.002))
                        FcspStamp = [FcspStamp tempSpikeTime_Unit1_Go{itr}(iSpike)];
                    end
                end
                FcspStamp(FcspStamp<BaseLen_SpikeRaster+SampOdorLen+iBin-1 | FcspStamp>=BaseLen_SpikeRaster+SampOdorLen+iBin) = [];
                FcspNum_Go = [FcspNum_Go; numel(FcspStamp)];
            end
            
            %% number of FCSP events of each NoGo trial
            FcspNum_NoGo = [];
            for itr = 1:length(tempSpikeTime_Unit1_NoGo)
                FcspStamp = [];
                for iSpike = 1:length(tempSpikeTime_Unit1_NoGo{itr})
                    if ~isempty(find(tempSpikeTime_Unit2_NoGo{itr}(:,1)-tempSpikeTime_Unit1_NoGo{itr}(iSpike)<=0.01 & tempSpikeTime_Unit2_NoGo{itr}(:,1)-tempSpikeTime_Unit1_NoGo{itr}(iSpike)>0.002))
                        FcspStamp = [FcspStamp tempSpikeTime_Unit1_NoGo{itr}(iSpike)];
                    end
                end
                FcspStamp(FcspStamp < BaseLen_SpikeRaster+SampOdorLen+iBin-1 | FcspStamp >= BaseLen_SpikeRaster+SampOdorLen+iBin) = [];
                FcspNum_NoGo = [FcspNum_NoGo; numel(FcspStamp)];
            end
            
            %% auROC of FCSP events
            [DV_Go,DV_Nogo] = DecisionVariableCalculation(FcspNum_Go,FcspNum_NoGo);
            [TPR,FPR] = TprFprCalculation(DV_Go,DV_Nogo);
            tempAUC = AucAnalysis(FPR,TPR);
            
            %% Rastergram plot
            FCperiod = BaseLen_SpikeRaster + SampOdorLen + iBin;
            if nnz(FcspNum_Go>=2) >= ShownTriNum_Go
                % ID of Go trials
                [~,sortid] = sortrows(FcspNum_Go,'descend');
                SelectedTrialID_Go = sortid(1:ShownTriNum_Go);
                SelectedFcspNum_Go = FcspNum_Go(SelectedTrialID_Go);
                % ID of NoGo trials
                id = randperm(numel(FcspNum_NoGo));
                SelectedTrialID_NoGo = id(1:ShownTriNum_NoGo);
                SelectedFcspNum_NoGo = FcspNum_NoGo(SelectedTrialID_NoGo);
                % spike raster
                figure('Color','w','Position',[500,200,500,600])
                SpikeData_LU = horzcat(tempSpikeTime_Unit1_Go(:,SelectedTrialID_Go),tempSpikeTime_Unit1_NoGo(:,SelectedTrialID_NoGo));
                SpikeData_FU = horzcat(tempSpikeTime_Unit2_Go(:,SelectedTrialID_Go),tempSpikeTime_Unit2_NoGo(:,SelectedTrialID_NoGo));
                SpikeData = vertcat(SpikeData_LU,SpikeData_FU);
                FcspNum = horzcat(SelectedFcspNum_Go(:)',SelectedFcspNum_NoGo(:)');
                for iTrial = 1:size(SpikeData,2)
                    for iSpike = 1:numel(SpikeData{1,iTrial}) % leading neuron
                        if ~isempty(find(SpikeData{2,iTrial} - SpikeData{1,iTrial}(iSpike) <= 0.01 & SpikeData{2,iTrial} - SpikeData{1,iTrial}(iSpike) > 0.002))
                            plot([SpikeData{1,iTrial}(iSpike) SpikeData{1,iTrial}(iSpike)], [1.5+3*(size(SpikeData,2)-iTrial) 2.5+3*(size(SpikeData,2)-iTrial)],'color',C2);
                        else
                            plot([SpikeData{1,iTrial}(iSpike) SpikeData{1,iTrial}(iSpike)], [1.5+3*(size(SpikeData,2)-iTrial) 2.5+3*(size(SpikeData,2)-iTrial)],'color',C1);
                        end
                        hold on
                    end
                    for iSpike = 1:numel(SpikeData{2,iTrial}) % following neuron
                        if ~isempty(find(SpikeData{2,iTrial}(iSpike) - SpikeData{1,iTrial} <= 0.01 & SpikeData{2,iTrial}(iSpike) - SpikeData{1,iTrial} > 0.002))
                            plot([SpikeData{2,iTrial}(iSpike) SpikeData{2,iTrial}(iSpike)], [0.5+3*(size(SpikeData,2)-iTrial) 1.5+3*(size(SpikeData,2)-iTrial)],'color',C3);
                        else
                            plot([SpikeData{2,iTrial}(iSpike) SpikeData{2,iTrial}(iSpike)], [0.5+3*(size(SpikeData,2)-iTrial) 1.5+3*(size(SpikeData,2)-iTrial)],'color',C1);
                        end
                        hold on
                    end
                    text(FCperiod-0.2,1.5+3*(size(SpikeData,2)-iTrial),num2str(FcspNum(iTrial)),'color','red'); hold on
                end
                box off;
                title(['auROC = ' num2str(tempAUC)]);
                if IsConsiderNoGoTrial == 1
                    set(gca,'XTick',BaseLen_SpikeRaster+SampOdorLen:1:BaseLen_SpikeRaster+SampOdorLen+DelayLen,'XTickLabel',{'','','','',''},'YTick',1.5:3:2.75+3*(numel(FcspNum)-1),'FontName','Arial','FontSize',16,'XLim',[BaseLen_SpikeRaster+SampOdorLen+iBin-1 BaseLen_SpikeRaster+SampOdorLen+iBin],'YLim',[0.25 2.75+3*(numel(FcspNum)-1)]);
                elseif IsConsiderNoGoTrial == 0
                    set(gca,'XTick',BaseLen_SpikeRaster+SampOdorLen:1:BaseLen_SpikeRaster+SampOdorLen+DelayLen,'XTickLabel',{'','','','',''},'YTick',1.5:3:2.75+3*(numel(FcspNum)-1),'FontName','Arial','FontSize',16,'XLim',[BaseLen_SpikeRaster+SampOdorLen+iBin-1 BaseLen_SpikeRaster+SampOdorLen+iBin],'YLim',[0.25+3*(numel(FcspNum)-ShownTriNum_Go) 2.75+3*(numel(FcspNum)-1)]);
                end
                set(gcf,'Render','Painter'); saveas(gcf,sprintf('FCauROC-%s-FCSP in %s#%d-%s#%d_at delay %dth s',Type,Reg_unit1,IDinReg_unit1,Reg_unit2,IDinReg_unit2,iBin),'fig'); close all;
            end
        end
    end
end









