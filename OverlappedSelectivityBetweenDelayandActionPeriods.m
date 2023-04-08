%% Proportion of neurons showing selectivity for hit vs. CR during action period and response period for STM-encoding neurons.

clear; clc; close all;

%% Assignment
BrainRegion = {'mPFC','aAIC'};
PreferType = {'GoPreferredUnits','NoGoPreferredUnits'};
BinsNuminWindow = 5;
BaseLen = 4;
SampOdorLen = 1;
DelayLen = 4;
TestOdorLen = 0.5;
RespLen = 0.5;
TimeGain = 10;
UnitsID_Memory_CongCodingHitCR_TestOdor = cell(1,numel(BrainRegion));
UnitsID_Memory_IncongCodingHitCR_TestOdor = cell(1,numel(BrainRegion));
UnitsID_Memory_NotCodingHitCR_TestOdor = cell(1,numel(BrainRegion));
UnitsID_Memory_CongCodingHitCR_Resp = cell(1,numel(BrainRegion));
UnitsID_Memory_IncongCodingHitCR_Resp = cell(1,numel(BrainRegion));
UnitsID_Memory_NotCodingHitCR_Resp = cell(1,numel(BrainRegion));

%% Enter path
CurrPath = uigetdir;
cd(CurrPath);

%% Load trialmark information for each neuron
load('NeuronOriginandSpikeRasterInformationforCtrlGroup.mat');

%% Whether neuron shows selectivity in action period
Prop = [];
for iReg = 1:numel(BrainRegion)
    UnitsID_Memory_CongCodingHitCR_TestOdor{iReg} = cell(numel(PreferType),1);
    UnitsID_Memory_IncongCodingHitCR_TestOdor{iReg} = cell(numel(PreferType),1);
    UnitsID_Memory_NotCodingHitCR_TestOdor{iReg} = cell(numel(PreferType),1);
    UnitsID_Memory_CongCodingHitCR_Resp{iReg} = cell(numel(PreferType),1);
    UnitsID_Memory_IncongCodingHitCR_Resp{iReg} = cell(numel(PreferType),1);
    UnitsID_Memory_NotCodingHitCR_Resp{iReg} = cell(numel(PreferType),1);
    load(sprintf('%sFRinS1S2_CtrlGroup.mat',BrainRegion{iReg}));
    if strcmp(BrainRegion{iReg},'mPFC')
        TrialMark = mPFCUnitTrialMark;
    elseif strcmp(BrainRegion{iReg},'aAIC')
        TrialMark = aAICUnitTrialMark;
    end
    for iNeuronType = 1:numel(PreferType)
        % ID of memory neurons (Go- and NoGo-preferred)
        if strcmp(BrainRegion{iReg},'mPFC') & strcmp(PreferType{iNeuronType},'GoPreferredUnits')
            tempfilename = dir('*ID-mPFC-GoPreferredUnits*.mat');
        elseif strcmp(BrainRegion{iReg},'mPFC') & strcmp(PreferType{iNeuronType},'NoGoPreferredUnits')
            tempfilename = dir('*ID-mPFC-NoGoPreferredUnits*.mat');
        elseif strcmp(BrainRegion{iReg},'aAIC') & strcmp(PreferType{iNeuronType},'GoPreferredUnits')
            tempfilename = dir('*ID-aAIC-GoPreferredUnits*.mat');
        elseif strcmp(BrainRegion{iReg},'aAIC') & strcmp(PreferType{iNeuronType},'NoGoPreferredUnits')
            tempfilename = dir('*ID-aAIC-NoGoPreferredUnits*.mat');
        end
        load(tempfilename.name);
        UnitsID_memory = tempNeuID;
        % ID of neurons with selectivity for hit vs. CR (not considering Go and No-Go preferred)
        if strcmp(BrainRegion{iReg},'mPFC')
            tempfilename = dir('*mPFCData_SelectivityforHitandCR*.mat');
        elseif strcmp(BrainRegion{iReg},'aAIC')
            tempfilename = dir('*aAICData_SelectivityforHitandCR*.mat');
        end
        load(tempfilename.name);
        UnitsID_codingHitCR = CodingUnitsID;
        UnitsIDinRange_codingHitCR = CodingUnitsID_NeuronsInRange;
        % ID of overlapped neurons between memory neurons and selective for hit vs. CR (action period) neurons
        TestOdorWindowID = (BaseLen+SampOdorLen+DelayLen+TestOdorLen)*10/BinsNuminWindow;
        RespWindowID = (BaseLen+SampOdorLen+DelayLen+TestOdorLen+RespLen)*10/BinsNuminWindow;
        for iUnit = UnitsID_memory(:)'
            % trialmark
            tempTrialMark = TrialMark{iUnit};
            tempTrialMark_Go = tempTrialMark(tempTrialMark(:,2)==1,:);
            tempHitIDinGo = find(tempTrialMark_Go(:,4)==1);
            tempTrialMark_NoGo = tempTrialMark(tempTrialMark(:,2)==2,:);
            tempCRIDinNoGo = find(tempTrialMark_NoGo(:,4)==4);
            % FR of single neuron % baseline: 1 s; Sample: 1 s; Delay: 4 s; Test: 0.5 s; Response window: 0.5 s
            tempUnitsFR_Hit = TargetBrainUnitsFRinS1{iUnit}(tempHitIDinGo,11:end);
            tempUnitsFR_CR = TargetBrainUnitsFRinS2{iUnit}(tempCRIDinNoGo,11:end);
            pos = find(UnitsID_codingHitCR == iUnit);
            if ~isempty(pos)
                tempUnitID = UnitsIDinRange_codingHitCR(pos);
                tempSigBinsID = SigBinsID{tempUnitID};
                % whether neuron showed selectivity for hit vs. CR during the test odor window
                pos_testodor = find(tempSigBinsID == TestOdorWindowID);
                if ~isempty(pos_testodor)
                    % difference in FR during test odor period of hit and CR trials
                    FRDiffBetwHitCR = mean(mean(tempUnitsFR_Hit(:,61:65))) - mean(mean(tempUnitsFR_CR(:,61:65)));
                    if (iNeuronType == 1 & FRDiffBetwHitCR > 0) | (iNeuronType == 2 & FRDiffBetwHitCR < 0)
                        UnitsID_Memory_CongCodingHitCR_TestOdor{iReg}{iNeuronType} = [UnitsID_Memory_CongCodingHitCR_TestOdor{iReg}{iNeuronType} iUnit];
                    else
                        UnitsID_Memory_IncongCodingHitCR_TestOdor{iReg}{iNeuronType} = [UnitsID_Memory_IncongCodingHitCR_TestOdor{iReg}{iNeuronType} iUnit];
                    end
                else
                    UnitsID_Memory_NotCodingHitCR_TestOdor{iReg}{iNeuronType} = [UnitsID_Memory_NotCodingHitCR_TestOdor{iReg}{iNeuronType} iUnit];
                end
                % whether neuron shows selectivity for hit vs. CR during response window
                pos_resp = find(tempSigBinsID == RespWindowID);
                if ~isempty(pos_resp)
                    % difference in FR during test odor period of hit and CR trials
                    FRDiffBetwHitCR = mean(mean(tempUnitsFR_Hit(:,66:70))) - mean(mean(tempUnitsFR_CR(:,66:70)));
                    if (iNeuronType == 1 & FRDiffBetwHitCR > 0) | (iNeuronType == 2 & FRDiffBetwHitCR < 0)
                        UnitsID_Memory_CongCodingHitCR_Resp{iReg}{iNeuronType} = [UnitsID_Memory_CongCodingHitCR_Resp{iReg}{iNeuronType} iUnit];
                    else
                        UnitsID_Memory_IncongCodingHitCR_Resp{iReg}{iNeuronType} = [UnitsID_Memory_IncongCodingHitCR_Resp{iReg}{iNeuronType} iUnit];
                    end
                else
                    UnitsID_Memory_NotCodingHitCR_Resp{iReg}{iNeuronType} = [UnitsID_Memory_NotCodingHitCR_Resp{iReg}{iNeuronType} iUnit];
                end
            else
                UnitsID_Memory_NotCodingHitCR_TestOdor{iReg}{iNeuronType} = [UnitsID_Memory_NotCodingHitCR_TestOdor{iReg}{iNeuronType} iUnit];
                UnitsID_Memory_NotCodingHitCR_Resp{iReg}{iNeuronType} = [UnitsID_Memory_NotCodingHitCR_Resp{iReg}{iNeuronType} iUnit];
            end
           %% PSTH plots
            figure('OuterPosition',[219 303 420 534]);
            plot((1:size(tempUnitsFR_Hit,2))/TimeGain, smooth(mean(tempUnitsFR_Hit,1),3)', 'color',[1 0 0],'linewidth',2);
            hold on
            plot((1:size(tempUnitsFR_CR,2))/TimeGain, smooth(mean(tempUnitsFR_CR,1),3)', 'color',[0 0 1],'linewidth',2);
            hold on
            patch([1 1 2 2],[0 5 5 0],'k','FaceAlpha',0.2,'edgecolor','none');
            hold on
            patch([6 6 6.5 6.5],[0 5 5 0],'k','FaceAlpha',0.2,'edgecolor','none');
            hold on
            plotshadow(tempUnitsFR_Hit,[1 0 0],2,3,0);
            hold on
            plotshadow(tempUnitsFR_CR,[0 0 1],2,3,0);
            % label time bins when selectivity for hit vs. CR is significant
            if ~isempty(pos)
                tempUnitID = UnitsIDinRange_codingHitCR(pos);
                tempSigBinsID = SigBinsID{tempUnitID}-3*TimeGain/BinsNuminWindow; % 3: duration of baseline period (4 s) - duration of baseline in FR result (1 s)
            else
                tempSigBinsID = [];
            end
            if ~isempty(tempSigBinsID)
                hold on
                for j=1:length(tempSigBinsID)
                    patch([tempSigBinsID(j)-1 tempSigBinsID(j)-1 tempSigBinsID(j) tempSigBinsID(j)]./(TimeGain/BinsNuminWindow),[1 2 2 1],[0 0 0],'edgecolor','none');
                    hold on
                end
            end
            set(gca,'XTick',1:1:9,'XTickLabel',{'0','','2','','4','','6','','8'},'FontName','Arial','FontSize',16,'xlim',[0.1 7]);
            % set(gca,'YTick',[0 TrialNum size(SingleUnitRG,1)],'YTickLabel',{'','',''},'FontName','Arial','FontSize',16,'ylim',[0 size(SingleUnitRG,1)+3]);
            box off;
            set(gcf, 'Renderer', 'Painter'); saveas(gcf,[BrainRegion{iReg} '-Unit ID-' num2str(iUnit) '-PSTH for hit and CR trials'],'fig'); close;
        end
        tempRatio_InconCoding_testodor = numel(UnitsID_Memory_CongCodingHitCR_TestOdor{iReg}{iNeuronType}) / numel(UnitsID_memory);
        tempRatio_NotCoding_testodor = numel(UnitsID_Memory_NotCodingHitCR_TestOdor{iReg}{iNeuronType}) / numel(UnitsID_memory);
        tempRatio_ConCoding_testodor = 1 - tempRatio_InconCoding_testodor - tempRatio_NotCoding_testodor;
        tempRatio_InconCoding_resp = numel(UnitsID_Memory_CongCodingHitCR_Resp{iReg}{iNeuronType}) / numel(UnitsID_memory);
        tempRatio_NotCoding_resp = numel(UnitsID_Memory_NotCodingHitCR_Resp{iReg}{iNeuronType}) / numel(UnitsID_memory);
        tempRatio_ConCoding_resp = 1 - tempRatio_InconCoding_resp - tempRatio_NotCoding_resp;
        figure('position',[200 200 600 400]);
        y = [tempRatio_ConCoding_testodor tempRatio_NotCoding_testodor tempRatio_InconCoding_testodor; tempRatio_ConCoding_resp tempRatio_NotCoding_resp tempRatio_InconCoding_resp];
        bar(y,'stacked');
        set(gca,'XTick',zeros(1,0),'XLim',[0.5 2.5],'YTick',0:0.2:1,'YTickLabel',{'0','20','40','60','80','100'},'YLim',[0 1]);
        box off
        set(gcf,'Render','Painter'); saveas(gcf,sprintf('Proportion of memory neurons losing and changing selectivity for Hit and CR after delay period in memory neurons-%s-%s',PreferType{iNeuronType},BrainRegion{iReg}),'fig'); close all;
    end
end
save('UnitIDKeepLoseSwitchSelectivityforHitCR_GoandNoGoPreferredUnits_mPFCandaAIC','UnitsID_Memory_CongCodingHitCR_TestOdor','UnitsID_Memory_IncongCodingHitCR_TestOdor','UnitsID_Memory_NotCodingHitCR_TestOdor','UnitsID_Memory_CongCodingHitCR_Resp','UnitsID_Memory_IncongCodingHitCR_Resp','UnitsID_Memory_NotCodingHitCR_Resp','-v7.3');



