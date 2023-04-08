%% Fano factor and distribution of FR across trials for sustained neurons.

clear; clc; close all;

%% Assignment
TargetBrain = 'mPFC';
BinNum = 20;

%% Load results of spike raster, FR, and selectivity
load('NeuronOriginandSpikeRasterInformationforCtrlGroup.mat');
if strcmp(TargetBrain,'mPFC')
    load('mPFCSelectivityData_CtrlGroup.mat'); 
    TargetBrainSustainedUnitID = SustainedUnitID; 
    load('mPFCFRinS1S2_CtrlGroup.mat');  % FR
else
    load('aAICSelectivityData_CtrlGroup.mat');
    TargetBrainSustainedUnitID = SustainedUnitID; 
    load('aAICFRinS1S2_CtrlGroup.mat'); % FR
end
SustainedUnitFRinS1Trials = TargetBrainUnitsFRinS1(:,TargetBrainSustainedUnitID);
SustainedUnitFRinS1Trials = cellfun(@(x) x(:,11:end),SustainedUnitFRinS1Trials,'UniformOutput', 0);
SustainedUnitFRinS2Trials = TargetBrainUnitsFRinS2(:,TargetBrainSustainedUnitID);
SustainedUnitFRinS2Trials = cellfun(@(x) x(:,11:end),SustainedUnitFRinS2Trials,'UniformOutput', 0);
SustainedUnitsNum = size(SustainedUnitFRinS1Trials,2);
% ID of neurons in SustainedUnitSigSelectivity
UnitID = [];
for iUnit = 1:SustainedUnitsNum
    UnitID = [UnitID find(SortedID==TargetBrainSustainedUnitID(iUnit))];
end
SustainedUnitSigSelectivity = TargetBrainSigSelectivity(UnitID,:);
clearvars -except TargetBrain group TargetBrainSustainedUnitID BinNum SustainedUnitFRinS1Trials SustainedUnitFRinS2Trials SustainedUnitSigSelectivity SustainedUnitsNum
% /// These variables are used to compute fano factor, plot firing rate in each trial and get proportion of trials with different FR value in S1 and S2 trias during early and late delay period ///
AllSustainedUnitsFiringVariation=struct('FanoFactor',[],'FRinTwoTrials',[],'RatioOfTrialsWithVariousFR',[],'IsSigTrialBasedDifference',zeros(SustainedUnitsNum,2),'ID',TargetBrainSustainedUnitID);
AllColor = [];
for iSustainedUnit = 1:SustainedUnitsNum
    
    %% FR of individual trial
    tempFRinS1 = SustainedUnitFRinS1Trials{1,iSustainedUnit}; 
    tempFRinS2 = SustainedUnitFRinS2Trials{1,iSustainedUnit};
    for iTrial = 1:size(tempFRinS1,1)
        tempFRinS1(iTrial,:) = smooth(tempFRinS1(iTrial,:),5)';
    end
    for iTrial = 1:size(tempFRinS2,1)
        tempFRinS2(iTrial,:) = smooth(tempFRinS2(iTrial,:),5)';
    end
    AllSustainedUnitsFiringVariation.FRinTwoTrials{1,iSustainedUnit} = tempFRinS1;
    AllSustainedUnitsFiringVariation.FRinTwoTrials{2,iSustainedUnit} = tempFRinS2;
    
    %% Fano factor
    tempVariance = [];
    tempMeanValue = [];
    if SustainedUnitSigSelectivity(iSustainedUnit,40) > 0
        tempMeanValue = mean(tempFRinS1);
        tempVariance = var(tempFRinS1);
        tempColor = {[1 0 0]};
    else
        tempMeanValue = mean(tempFRinS2);
        tempVariance = var(tempFRinS2);
        tempColor = {[0 0 1]};
    end
    AllColor = [AllColor; tempColor];
    AllSustainedUnitsFiringVariation.FanoFactor = [AllSustainedUnitsFiringVariation.FanoFactor; tempVariance./tempMeanValue];
    
    %% Proportion of trials with different FR during early and late delay
    % S1 trials
    meanFRDuringDelayinS1 = []; 
    for iTrial = 1:size(SustainedUnitFRinS1Trials{1,iSustainedUnit},1)
        meanFRDuringDelayinS1(iTrial,1) = mean(SustainedUnitFRinS1Trials{1,iSustainedUnit}(iTrial,21:40),2); % early delay
        meanFRDuringDelayinS1(iTrial,2) = mean(SustainedUnitFRinS1Trials{1,iSustainedUnit}(iTrial,41:60),2); % late delay
    end
    % S2 trials
    meanFRDuringDelayinS2 = [];
    for iTrial = 1:size(SustainedUnitFRinS2Trials{1,iSustainedUnit},1)
        meanFRDuringDelayinS2(iTrial,1) = mean(SustainedUnitFRinS2Trials{1,iSustainedUnit}(iTrial,21:40),2); % early delay
        meanFRDuringDelayinS2(iTrial,2) = mean(SustainedUnitFRinS2Trials{1,iSustainedUnit}(iTrial,41:60),2); % late delay
    end
    AllSustainedUnitsFiringVariation.RatioOfTrialsWithVariousFR{1,iSustainedUnit} = meanFRDuringDelayinS1;
    AllSustainedUnitsFiringVariation.RatioOfTrialsWithVariousFR{2,iSustainedUnit} = meanFRDuringDelayinS2;
end
AllSustainedUnitsFiringVariation.FanoFactor(isnan(AllSustainedUnitsFiringVariation.FanoFactor)) = 0;

%% Plot fano factor for each sustained neuron
X = -1:0.1:6;
X = X + 0.05;
X = X(1:end-1);
for iSustainedUnit = 1:SustainedUnitsNum
    figure;
    plot(X,AllSustainedUnitsFiringVariation.FanoFactor(iSustainedUnit,:),'color',AllColor{iSustainedUnit},'linewidth',3);
    hold on 
    plot([0 0],[min(AllSustainedUnitsFiringVariation.FanoFactor(iSustainedUnit,:))-0.2 max(AllSustainedUnitsFiringVariation.FanoFactor(iSustainedUnit,:))+0.2],'k--');
    plot([1 1],[min(AllSustainedUnitsFiringVariation.FanoFactor(iSustainedUnit,:))-0.2 max(AllSustainedUnitsFiringVariation.FanoFactor(iSustainedUnit,:))+0.2],'k--');
    plot([5 5],[min(AllSustainedUnitsFiringVariation.FanoFactor(iSustainedUnit,:))-0.2 max(AllSustainedUnitsFiringVariation.FanoFactor(iSustainedUnit,:))+0.2],'k--');
    plot([5.5 5.5],[min(AllSustainedUnitsFiringVariation.FanoFactor(iSustainedUnit,:))-0.2 max(AllSustainedUnitsFiringVariation.FanoFactor(iSustainedUnit,:))+0.2],'k--');
    axis([-0.5 6 min(AllSustainedUnitsFiringVariation.FanoFactor(iSustainedUnit,:))-0.2 max(AllSustainedUnitsFiringVariation.FanoFactor(iSustainedUnit,:))+0.2]);
    set(gca,'XTick',0:1:6,'XTickLabel',{'0','1','2','3','4','5','6'},'FontName','Arial','FontSize',16);
    set(gcf,'Renderer','Painter'); saveas(gcf,[TargetBrain '_FanoFactor_SustainedUnit#' num2str(TargetBrainSustainedUnitID(iSustainedUnit)) ],'fig'); close;
end

%% Plot FR in individual trial for each sustained neuron
for iSustainedUnit = 1:SustainedUnitsNum
    figure;
    for iTrial = 1:size(AllSustainedUnitsFiringVariation.FRinTwoTrials{1,iSustainedUnit},1) % S1 trials
        plot(X,AllSustainedUnitsFiringVariation.FRinTwoTrials{1,iSustainedUnit}(iTrial,:),'color',[230 143 144]/255,'linewidth',1);
        hold on
    end
    for iTrial = 1:size(AllSustainedUnitsFiringVariation.FRinTwoTrials{2,iSustainedUnit}) % S2 trials
        plot(X,AllSustainedUnitsFiringVariation.FRinTwoTrials{2,iSustainedUnit}(iTrial,:),'color',[142 144 192]/255,'linewidth',1);
        hold on
    end
    plot(X,mean(AllSustainedUnitsFiringVariation.FRinTwoTrials{1,iSustainedUnit}),'color',[1 0 0],'linewidth',3);
    plot(X,mean(AllSustainedUnitsFiringVariation.FRinTwoTrials{2,iSustainedUnit}),'color',[0 0 1],'linewidth',3);
    Ymin = floor(min(min(min(AllSustainedUnitsFiringVariation.FRinTwoTrials{1,iSustainedUnit})),min(min(AllSustainedUnitsFiringVariation.FRinTwoTrials{2,iSustainedUnit}))))-2;
    Ymax = ceil(max(max(max(AllSustainedUnitsFiringVariation.FRinTwoTrials{1,iSustainedUnit})),max(max(AllSustainedUnitsFiringVariation.FRinTwoTrials{2,iSustainedUnit}))))+2;
    axis([-0.5 6  Ymin Ymax]);
    hold on 
    YStepTick = ceil((Ymax-Ymin+4)/5);
    plot([0 0],[Ymin-0.2 Ymax+0.2],'k--');
    plot([1 1],[Ymin-0.2 Ymax+0.2],'k--');
    plot([5 5],[Ymin-0.2 Ymax+0.2],'k--');
    plot([5.5 5.5],[Ymin-0.2 Ymax+0.2],'k--');
    set(gca,'XTick',0:1:6,'XTickLabel',{'0','1','2','3','4','5','6'},'FontName','Arial','FontSize',16);
    set(gca,'YTick',Ymin:YStepTick:Ymax,'YTickLabel',{num2str(Ymin),num2str(Ymin+YStepTick),num2str(Ymin+YStepTick*2),num2str(Ymin+YStepTick*3),num2str(Ymin+YStepTick*4)},'FontName','Arial','FontSize',16);
    set(gcf,'Renderer','Painter'); saveas(gcf,[TargetBrain '_SingleTrialFiringRate_SustainedUnit#' num2str(TargetBrainSustainedUnitID(iSustainedUnit)) ],'fig'); close;
end

%% Histogram plot for distribution of trial with differernt FR
p = [];
startcenter = 0.5;
for iSustainedUnit = 1:SustainedUnitsNum
    S1TrialFRinEarlyDelay = AllSustainedUnitsFiringVariation.RatioOfTrialsWithVariousFR{1,iSustainedUnit}(:,1);
    S2TrialFRinEarlyDelay = AllSustainedUnitsFiringVariation.RatioOfTrialsWithVariousFR{2,iSustainedUnit}(:,1);
    S1TrialFRinLateDelay = AllSustainedUnitsFiringVariation.RatioOfTrialsWithVariousFR{1,iSustainedUnit}(:,2); 
    S2TrialFRinLateDelay = AllSustainedUnitsFiringVariation.RatioOfTrialsWithVariousFR{2,iSustainedUnit}(:,2);
    MaxFRinEarlyDelay = ceil(max(max(S1TrialFRinEarlyDelay),max(S2TrialFRinEarlyDelay)));
    MinFRinEarlyDelay = floor(min(min(S1TrialFRinEarlyDelay),min(S2TrialFRinEarlyDelay)));
    MaxFRinLateDelay = ceil(max(max(S1TrialFRinLateDelay),max(S2TrialFRinLateDelay)));
    MinFRinLateDelay = floor(min(min(S1TrialFRinLateDelay),min(S2TrialFRinLateDelay)));
    figure;
    [Counts,Centers] = hist(S1TrialFRinEarlyDelay, startcenter:(MaxFRinEarlyDelay-startcenter)/(BinNum-1):MaxFRinEarlyDelay);
    bar(Centers,Counts/length(S1TrialFRinEarlyDelay),'FaceColor',[1 0 0],'FaceAlpha',0.5,'EdgeColor','none');
    hold on
    [Counts,Centers] = hist(S2TrialFRinEarlyDelay, startcenter:(MaxFRinEarlyDelay-startcenter)/(BinNum-1):MaxFRinEarlyDelay);
    bar(Centers,Counts/length(S2TrialFRinEarlyDelay),'FaceColor',[0 0 1],'FaceAlpha',0.5,'EdgeColor','none');
    view(90,-90);set(gca,'XDir','reverse');title('Early delay');
    set(gcf,'Renderer','Painter'); saveas(gcf,[TargetBrain 'TrialsDistributionofFR_SustainedUnit#' num2str(TargetBrainSustainedUnitID(iSustainedUnit)) '_EarlyDelay'],'fig'); close;
    figure;
    [Counts,Center] = hist(S1TrialFRinLateDelay, startcenter:(MaxFRinLateDelay-startcenter)/(BinNum-1):MaxFRinLateDelay);
    bar(Centers,Counts/length(S1TrialFRinLateDelay),'FaceColor',[1 0 0],'FaceAlpha',0.5,'EdgeColor','none');
    hold on
    [Counts,Center] = hist(S2TrialFRinLateDelay, startcenter:(MaxFRinLateDelay-startcenter)/(BinNum-1):MaxFRinLateDelay);
    bar(Centers,Counts/length(S2TrialFRinLateDelay),'FaceColor',[0 0 1],'FaceAlpha',0.5,'EdgeColor','none');
    view(90,-90);set(gca,'XDir','reverse');title('Late delay');
    set(gcf,'Renderer','Painter'); saveas(gcf,[TargetBrain 'TrialsDistributionofFR_SustainedUnit#' num2str(TargetBrainSustainedUnitID(iSustainedUnit)) '_LateDelay'],'fig'); close;
    % whether FR between different sample trials is different
    LowBoundaryForBiggerFR = []; 
    MaxPercentileforBiggerFRTrials = []; 
    BelowHighBoundaryTrialID =[];
    BeyondLowBoundaryTrialID = []; 
    for iDelay = 1:2  % early and late delay
        if AllColor{iSustainedUnit} == [1 0 0]
            BiggerFiring = AllSustainedUnitsFiringVariation.RatioOfTrialsWithVariousFR{1,iSustainedUnit};
            SmallerFiring = AllSustainedUnitsFiringVariation.RatioOfTrialsWithVariousFR{2,iSustainedUnit};
        else
            BiggerFiring = AllSustainedUnitsFiringVariation.RatioOfTrialsWithVariousFR{2,iSustainedUnit};
            SmallerFiring = AllSustainedUnitsFiringVariation.RatioOfTrialsWithVariousFR{1,iSustainedUnit};
        end
        LowBoundaryForBiggerFR(iDelay) = prctile(BiggerFiring(:,iDelay),2.5);
        HighBoundaryForSmallerFR(iDelay) = prctile(SmallerFiring(:,iDelay),97.5);
        BelowHighBoundaryTrialID{iDelay} =  find(BiggerFiring(:,iDelay)<HighBoundaryForSmallerFR(iDelay));
        BeyondLowBoundaryTrialID{iDelay} =  find(SmallerFiring(:,iDelay)>LowBoundaryForBiggerFR(iDelay));
        p(iSustainedUnit,iDelay) = mean([(length(BelowHighBoundaryTrialID{iDelay})+1)/(size(BiggerFiring,1)+1) (length(BeyondLowBoundaryTrialID{iDelay})+1)/(size(SmallerFiring,1)+1)]);
        if p(iSustainedUnit,iDelay) <= 0.05
            AllSustainedUnitsFiringVariation.IsSigTrialBasedDifference(iSustainedUnit,iDelay) = 1;
        end
    end
end
save(['AllSustained' TargetBrain 'CrossTrialVariation' group],'AllSustainedUnitsFiringVariation');





