%% Example movement-responsible neurons

clear; clc; close all;

%% Assignment
Reg = 'aAIC';
ShownTrialNum = 20;
BaseLen_ER = 9;
BaseLen_Movement = 0.5;
SampOdorLen = 1;
DelayLen = 4;
TestOdorLen = 0.5;
TimeGain = 60;

%% Load movement-related result and firing result
load(sprintf('Units responseness to movement-%s.mat',Reg));
load(sprintf('%sSelectivityData_CtrlGroup.mat',Reg));

% %% Example trunk movement-responsible neurons
% % activated
% SingleTrialSpikeRasterAndMovement(Reg,'Trunk movement activated',SortedActivatedUnitsID_Trunk,TargetBrainRGinS1,FinishedTrialID,TrialMark,TrunkMovement,TrunkMovementSig,ShownTrialNum,BaseLen_ER,BaseLen_Movement,SampOdorLen,DelayLen,TestOdorLen,TimeGain);
% % inhibited
% SingleTrialSpikeRasterAndMovement(Reg,'Trunk movement inhibited',SortedInhibitedUnitsID_Trunk,TargetBrainRGinS1,FinishedTrialID,TrialMark,TrunkMovement,TrunkMovementSig,ShownTrialNum,BaseLen_ER,BaseLen_Movement,SampOdorLen,DelayLen,TestOdorLen,TimeGain);
% 
% %% Example nose movement-responsible neurons
% % activated
% SingleTrialSpikeRasterAndMovement(Reg,'Nose movement activated',SortedActivatedUnitsID_Nose,TargetBrainRGinS1,FinishedTrialID,TrialMark,NoseMovement,NoseMovementSig,ShownTrialNum,BaseLen_ER,BaseLen_Movement,SampOdorLen,DelayLen,TestOdorLen,TimeGain);
% % inhibited
% SingleTrialSpikeRasterAndMovement(Reg,'Nose movement inhibited',SortedInhibitedUnitsID_Nose,TargetBrainRGinS1,FinishedTrialID,TrialMark,NoseMovement,NoseMovementSig,ShownTrialNum,BaseLen_ER,BaseLen_Movement,SampOdorLen,DelayLen,TestOdorLen,TimeGain);
% 
% %% Example pupil dilation-responsible neurons
% % activated
% SingleTrialSpikeRasterAndMovement(Reg,'Pupil dilation activated',SortedActivatedUnitsID_Pupil,TargetBrainRGinS1,FinishedTrialID,TrialMark,PupilDilation,PupilDilationSig,ShownTrialNum,BaseLen_ER,BaseLen_Movement,SampOdorLen,DelayLen,TestOdorLen,TimeGain);
% inhibited
SingleTrialSpikeRasterAndMovement(Reg,'Pupil dilation inhibited',SortedInhibitedUnitsID_Pupil,TargetBrainRGinS1,FinishedTrialID,TrialMark,PupilDilation,PupilDilationSig,ShownTrialNum,BaseLen_ER,BaseLen_Movement,SampOdorLen,DelayLen,TestOdorLen,TimeGain);




function SingleTrialSpikeRasterAndMovement(Region,MovementName,MovementModulatedUnitsID,AllUnitsSpikeRasterInGoTrials,CompletedTrialID,TrialMark,Movement,MovementSig,ShownTrialNum,BaseLen_spike,BaseLen_move,SampOdorLen,DelayLen,TestOdorLen,TimeGain)
for iUnit = 1:numel(MovementModulatedUnitsID)
    fprintf('Processing %dth unit\n',iUnit);
    UnitID = MovementModulatedUnitsID(iUnit);
    % spike raster in Go trials
    SpikeRaster_Go = AllUnitsSpikeRasterInGoTrials{UnitID};
    % trunk movement in Go trials
    UnitMovement = Movement{UnitID};
    UnitMovement = UnitMovement(CompletedTrialID{UnitID},:);
    UnitMovement = UnitMovement(TrialMark{UnitID}(:,3)==1,:);
    % trunk movement significance in Go trials
    UnitMovementSig = MovementSig{UnitID};
    UnitMovementSig = UnitMovementSig(CompletedTrialID{UnitID},:);
    UnitMovementSig = UnitMovementSig(TrialMark{UnitID}(:,3)==1,:);
    % randomly select example trials
    ExampleTrialID = randperm(size(UnitMovement,1));
    ExampleTrialID = ExampleTrialID(1:ShownTrialNum);
    SpikeRasterInExampleTrial = SpikeRaster_Go(ExampleTrialID,:);
    UnitTrunkMovementInExampleTrial = UnitMovement(ExampleTrialID,:);
    UnitTrunkMovementSigInExampleTrial = UnitMovementSig(ExampleTrialID,:);
    % plot movement and spike raster in example Go trials
    figure('position',[150 250 2000 600]);
    for iTrial = 1:ShownTrialNum
        % movement
        subplot(6,ShownTrialNum/2,iTrial+ShownTrialNum*(ceil(2*iTrial/ShownTrialNum)-1));
        tempBM = UnitTrunkMovementInExampleTrial(iTrial,:);
        tempBM = tempBM/mean(tempBM(1:BaseLen_move*TimeGain));
        tempBaseLevel = prctile(tempBM,5);
        tempSD = std(tempBM);
        tempThreshold = tempBaseLevel + 1.5*tempSD;
        plot([1:1:numel(tempBM)]/TimeGain,tempBM,'k'); hold on
        plot([1 numel(tempBM)]/TimeGain,[tempThreshold tempThreshold],'--k'); hold on
        patch([BaseLen_move BaseLen_move BaseLen_move+SampOdorLen BaseLen_move+SampOdorLen],[0 5 5 0],'k','FaceAlpha',0.2,'edgecolor','none'); hold on
        patch([BaseLen_move+SampOdorLen+DelayLen BaseLen_move+SampOdorLen+DelayLen BaseLen_move+SampOdorLen+DelayLen+TestOdorLen BaseLen_move+SampOdorLen+DelayLen+TestOdorLen],[0 5 5 0],'k','FaceAlpha',0.2,'edgecolor','none'); hold on
        set(gca,'XTick',0:1:numel(tempBM)/TimeGain,'XTickLabel',num2cell(0:1:numel(tempBM)/TimeGain),'XLim',[0 BaseLen_move+SampOdorLen+DelayLen+TestOdorLen],'YLim',[0.8 1.2]);
        box off;
        % movement significance
        subplot(6,ShownTrialNum/2,iTrial+ShownTrialNum*(ceil(2*iTrial/ShownTrialNum)-1)+ShownTrialNum/2);
        plot(1:1:numel(tempBM),UnitTrunkMovementSigInExampleTrial(iTrial,:),'k'); hold on
        patch([BaseLen_move*TimeGain BaseLen_move*TimeGain (BaseLen_move+SampOdorLen)*TimeGain (BaseLen_move+SampOdorLen)*TimeGain],[-1 2 2 -1],'k','FaceAlpha',0.2,'edgecolor','none'); hold on
        patch([(BaseLen_move+SampOdorLen+DelayLen)*TimeGain (BaseLen_move+SampOdorLen+DelayLen)*TimeGain (BaseLen_move+SampOdorLen+DelayLen+TestOdorLen)*TimeGain (BaseLen_move+SampOdorLen+DelayLen+TestOdorLen)*TimeGain],[-1 2 2 -1],'k','FaceAlpha',0.2,'edgecolor','none'); hold on
        set(gca,'XTick',0:TimeGain:numel(tempBM),'XTickLabel',num2cell(0:1:numel(tempBM)/TimeGain),'XLim',[0 (BaseLen_move+SampOdorLen+DelayLen+TestOdorLen)*TimeGain],'YLim',[-1 2]);
        box off;
        % spike raster
        subplot(6,ShownTrialNum/2,iTrial+ShownTrialNum*(ceil(2*iTrial/ShownTrialNum)-1)+ShownTrialNum);
        for iSpike = 1:numel(SpikeRasterInExampleTrial{iTrial})
            plot([SpikeRasterInExampleTrial{iTrial}(iSpike) SpikeRasterInExampleTrial{iTrial}(iSpike)], [0.5 1.5],'k'); hold on
        end
        patch([BaseLen_spike BaseLen_spike BaseLen_spike+SampOdorLen BaseLen_spike+SampOdorLen],[0 2 2 0],'k','FaceAlpha',0.2,'edgecolor','none'); hold on
        patch([BaseLen_spike+SampOdorLen+DelayLen BaseLen_spike+SampOdorLen+DelayLen BaseLen_spike+SampOdorLen+DelayLen+TestOdorLen BaseLen_spike+SampOdorLen+DelayLen+TestOdorLen],[0 2 2 0],'k','FaceAlpha',0.2,'edgecolor','none'); hold on
        set(gca,'XTick',0:1:BaseLen_spike+SampOdorLen+DelayLen+TestOdorLen,'XTickLabel',num2cell(0:1:BaseLen_spike+SampOdorLen+DelayLen+TestOdorLen),'XLim',[BaseLen_spike-BaseLen_move BaseLen_spike+SampOdorLen+DelayLen+TestOdorLen],'YLim',[0 2]);
        box off
    end
    set(gcf,'Renderer','Painter'); saveas(gcf,sprintf('%s-%s-Rank-%d-Unit ID-%d in Go trials',Region,MovementName,iUnit,UnitID),'fig'); close all;
end
end

