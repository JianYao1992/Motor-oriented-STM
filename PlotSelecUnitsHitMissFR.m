function [Samp1CorrErrorTriDiffFRIsSig,Samp1CorrErrorTriDiffFRIsSig_BelowChance,CorrErrTriActivity] = PlotSelecUnitsHitMissFR(UnitsID...
    ,Samp1CorrTriFR,Samp1ErrorTriFR,X,SampOdorLen,DelayLen,TestOdorLen,RespLen,TriLen,TimeGain...
    ,BaselinePer,RaworNormFR,BootsStrapNum,Title,IsClusterBasedPermuTest,WorkerNum)    
    
%% Establish parallel pool
if IsClusterBasedPermuTest==1 && WorkerNum>1
    poolobj = gcp('nocreate');
    if isempty(poolobj)
        myCluster = parcluster('local'); myCluster.NumWorkers = WorkerNum; parpool(myCluster,WorkerNum);
    end   
end

%% Averaged FR of Hit and Miss trials for real and shuffled result
MeanFR_Hit = zeros(length(UnitsID),TriLen-1); % real 
MeanFR_Miss = zeros(length(UnitsID),TriLen-1);
ShuffleFR_Hit = zeros(length(UnitsID),1000,TriLen-1); % shuffle
ShuffleFR_Miss = zeros(length(UnitsID),1000,TriLen-1);
for iUnit = 1:length(UnitsID)
    if mod(iUnit,10)==0
      disp(['Current neuron number=' num2str(iUnit)])
    end
    tempUnitID = UnitsID(iUnit);
    tempUnit_Hit = Samp1CorrTriFR{tempUnitID};
    tempUnit_Miss = Samp1ErrorTriFR{tempUnitID};    
    % averaged FR in Hit and Miss trials
    [BSSamp1CorrTriFR,BSSamp1ErrorTriFR,ShuffleSamp1CorrTriFR,ShuffleSamp1ErrorTriFR] = ReSampCorrErrorTrialsForComparison(tempUnit_Hit,tempUnit_Miss,...
        BootsStrapNum,IsClusterBasedPermuTest,TriLen);
    MeanFR_Hit(iUnit,:) = BSSamp1CorrTriFR;
    MeanFR_Miss(iUnit,:) = BSSamp1ErrorTriFR;
    ShuffleFR_Hit(iUnit,:,:) = ShuffleSamp1CorrTriFR;
    ShuffleFR_Miss(iUnit,:,:) = ShuffleSamp1ErrorTriFR;
end

%% Raw or normalized FR
if RaworNormFR == 1
    % real
    Samp1CorrTriMeanFR_New = MeanFR_Hit;
    Samp1ErrorTriMeanFR_New = MeanFR_Miss;
    % shuffle
    ShuffleFR_Hit_New = mean(ShuffleFR_Hit,1);
    ShuffleFR_Hit_New = reshape(ShuffleFR_Hit_New,size(ShuffleFR_Hit_New,2),size(ShuffleFR_Hit_New,3));
    ShuffleFR_Miss_New = mean(ShuffleFR_Miss,1);
    ShuffleFR_Miss_New = reshape(ShuffleFR_Miss_New,size(ShuffleFR_Miss_New,2),size(ShuffleFR_Miss_New,3));
elseif RaworNormFR == 2
    % real
    Samp1CorrTriMeanFR_New = NormFrCalculation(MeanFR_Hit,BaselinePer);
    Samp1ErrorTriMeanFR_New = NormFrCalculation(MeanFR_Miss,BaselinePer);
    % shuffle
    ShuffleFR_Hit_New = zeros(size(ShuffleFR_Hit));
    ShuffleFR_Miss_New = zeros(size(ShuffleFR_Miss));
    for iUnit = 1:length(UnitsID) 
        [NormShuffleSamp1CorrTriFR,NormShuffleSamp1ErrorTriFR] = NormFrCalculationForShuffledData(ShuffleFR_Hit(iUnit,:,:),ShuffleFR_Miss(iUnit,:,:),TriLen,BaselinePer);
        ShuffleFR_Hit_New(iUnit,:,:) = NormShuffleSamp1CorrTriFR;
        ShuffleFR_Miss_New(iUnit,:,:) = NormShuffleSamp1ErrorTriFR;
    end
    ShuffleFR_Hit_New = mean(ShuffleFR_Hit_New,1); % mean shuffled FR averaging across units
    ShuffleFR_Hit_New = reshape(ShuffleFR_Hit_New,size(ShuffleFR_Hit_New,2),size(ShuffleFR_Hit_New,3));
    ShuffleFR_Miss_New = mean(ShuffleFR_Miss_New,1);
    ShuffleFR_Miss_New = reshape(ShuffleFR_Miss_New,size(ShuffleFR_Miss_New,2),size(ShuffleFR_Miss_New,3));
end

%% Difference in FR of Hit and Miss trials
% real
Sam1CorrErrorTriDiffFR = mean(Samp1ErrorTriMeanFR_New)-mean(Samp1CorrTriMeanFR_New);
% shuffle
Sam1ShuffleCorrErrorTriDiffFR = ShuffleFR_Miss_New-ShuffleFR_Hit_New;

%% Cluster-based permutation test
[Samp1CorrErrorTriDiffFRIsSig,Samp1CorrErrorTriDiffFRIsSig_BelowChance] = PermutationTest_ClusterBasedOrNot(Sam1CorrErrorTriDiffFR,Sam1ShuffleCorrErrorTriDiffFR,1,1000);

%% Figure
figure;
% plot real result
[YMin1,YMax1] = PlotMeanAndSEM(Samp1CorrTriMeanFR_New,X,[1 0 0],[-3 6.4],'-',1);
[YMin2,YMax2] = PlotMeanAndSEM(Samp1ErrorTriMeanFR_New,X,[26 72 157]/255,[-3 6.4],'-',1);
% maximal and minimal FR
YMax = max([YMax1 YMax2]);
YMin = min([YMin1 YMin2]);
hold on
% label ID of bins with significant difference of FR between FA and CR trials
for i = 1:length(Samp1CorrErrorTriDiffFRIsSig)
    if Samp1CorrErrorTriDiffFRIsSig(i) == 1
        plot([X(i)-1/TimeGain/2 X(i)+1/TimeGain/2],YMax*1.05*ones(1,2),'-','color',[0 0 0],'linewidth',2)
    end
end
hold on
for i = 1:length(Samp1CorrErrorTriDiffFRIsSig_BelowChance)
    if Samp1CorrErrorTriDiffFRIsSig_BelowChance(i) == 1
        plot([X(i)-1/TimeGain/2 X(i)+1/TimeGain/2],YMax*1.05*ones(1,2),'-','color',[0 0 0],'linewidth',2)
    end
end
hold on
AddEventCurve(SampOdorLen,DelayLen,TestOdorLen,RespLen,YMax*1.3,YMin-0.2,6);
hold on
text(SampOdorLen+0.3,YMax*1.1,['NeuNumber-' num2str(length(UnitsID))],'color',[1 0 0])
axis([-1 6 YMin-0.2 YMax*1.3])
xlabel('Time from sample onset £¨s£©','fontsize',12)
if RaworNormFR == 1
    ylabel('Firing rate (Hz)','fontsize',12);
else
    ylabel('Normalized FR','fontsize',12);
end
title(Title,'fontsize',14);

%% FR of neurons in hit, FA, and CR trials
CorrErrTriActivity = [{Samp1CorrTriMeanFR_New};{Samp1ErrorTriMeanFR_New}];
