function [Samp2CorrErrorTriDiffFRIsSig,Samp2CorrErrorTriDiffFRIsSig_BelowChance,CorrErrTriActivity] = PlotSelecUnitsPerfDependFR(UnitsID...
    ,Samp1CorrTriFR,Samp2CorrTriFR,Samp2ErrorTriFR,X,SampOdorLen,DelayLen,TestOdorLen,RespLen,TriLen,TimeGain...
    ,BaselinePer,RaworNormFR,BootsStrapNum,Title,IsClusterBasedPermuTest,WorkerNum)    
    
%% Establish parallel pool
if IsClusterBasedPermuTest == 1 && WorkerNum > 1
    poolobj = gcp('nocreate');
    if isempty(poolobj)
        myCluster = parcluster('local'); myCluster.NumWorkers = WorkerNum; parpool(myCluster,WorkerNum);
    end   
end

%% Averaged FR of hit, FA, and CR trials for real and shuffled result
MeanFR_Hit = zeros(length(UnitsID),TriLen-1); % real 
MeanFR_CR = zeros(length(UnitsID),TriLen-1);
MeanFR_FA = zeros(length(UnitsID),TriLen-1);
ShuffleFR_CR = zeros(length(UnitsID),1000,TriLen-1); % shuffle
ShuffleFR_FA = zeros(length(UnitsID),1000,TriLen-1);
for iUnit = 1:length(UnitsID)
    if mod(iUnit,10)==0
      disp(['Current neuron number=' num2str(iUnit)])
    end
    tempUnitID = UnitsID(iUnit);
    tempUnit_Hit = Samp1CorrTriFR{tempUnitID};
    tempUnit_CR = Samp2CorrTriFR{tempUnitID};
    tempUnit_FA = Samp2ErrorTriFR{tempUnitID};    
    % averaged FR in hit trials
    MeanFR_Hit(iUnit,:) = mean(tempUnit_Hit(:,1:end-1));
    % averaged FR in FA and CR trials
    [BSSamp2CorrTriFR,BSSamp2ErrorTriFR,ShuffleSamp2CorrTriFR,ShuffleSamp2ErrorTriFR] = ReSampCorrErrorTrialsForComparison(tempUnit_CR,tempUnit_FA,...
        BootsStrapNum,IsClusterBasedPermuTest,TriLen);
    MeanFR_CR(iUnit,:) = BSSamp2CorrTriFR;
    MeanFR_FA(iUnit,:) = BSSamp2ErrorTriFR;
    ShuffleFR_CR(iUnit,:,:) = ShuffleSamp2CorrTriFR;
    ShuffleFR_FA(iUnit,:,:) = ShuffleSamp2ErrorTriFR;
end

%% Raw or normalized FR
if RaworNormFR == 1
    % real
    Samp1CorrTriMeanFR_New = MeanFR_Hit;
    Samp2CorrTriMeanFR_New = MeanFR_CR;
    Samp2ErrorTriMeanFR_New = MeanFR_FA;
    % shuffle
    ShuffleFR_CR_New = mean(ShuffleFR_CR,1);
    ShuffleFR_CR_New = reshape(ShuffleFR_CR_New,size(ShuffleFR_CR_New,2),size(ShuffleFR_CR_New,3));
    ShuffleFR_FA_New = mean(ShuffleFR_FA,1);
    ShuffleFR_FA_New = reshape(ShuffleFR_FA_New,size(ShuffleFR_FA_New,2),size(ShuffleFR_FA_New,3));
elseif RaworNormFR == 2
    % real
    Samp1CorrTriMeanFR_New = NormFrCalculation(MeanFR_Hit,BaselinePer);
    Samp2CorrTriMeanFR_New = NormFrCalculation(MeanFR_CR,BaselinePer);
    Samp2ErrorTriMeanFR_New = NormFrCalculation(MeanFR_FA,BaselinePer);
    % shuffle
    ShuffleFR_CR_New = zeros(size(ShuffleFR_CR));
    ShuffleFR_FA_New = zeros(size(ShuffleFR_FA));
    for iUnit = 1:length(UnitsID) 
        [NormShuffleSamp2CorrTriFR,NormShuffleSamp2ErrorTriFR] = NormFrCalculationForShuffledData(ShuffleFR_CR(iUnit,:,:),ShuffleFR_FA(iUnit,:,:),TriLen,BaselinePer);
        ShuffleFR_CR_New(iUnit,:,:) = NormShuffleSamp2CorrTriFR;
        ShuffleFR_FA_New(iUnit,:,:) = NormShuffleSamp2ErrorTriFR;
    end
    ShuffleFR_CR_New = mean(ShuffleFR_CR_New,1); % mean shuffled FR averaging across units
    ShuffleFR_CR_New = reshape(ShuffleFR_CR_New,size(ShuffleFR_CR_New,2),size(ShuffleFR_CR_New,3));
    ShuffleFR_FA_New = mean(ShuffleFR_FA_New,1);
    ShuffleFR_FA_New = reshape(ShuffleFR_FA_New,size(ShuffleFR_FA_New,2),size(ShuffleFR_FA_New,3));
end

%% Difference in FR of correct and error trials
% real
Sam2CorrErrorTriDiffFR = mean(Samp2ErrorTriMeanFR_New)-mean(Samp2CorrTriMeanFR_New);
% shuffle
Sam2ShuffleCorrErrorTriDiffFR = ShuffleFR_FA_New-ShuffleFR_CR_New;

%% Cluster-based permutation test
[Samp2CorrErrorTriDiffFRIsSig,Samp2CorrErrorTriDiffFRIsSig_BelowChance] = PermutationTest_ClusterBasedOrNot(Sam2CorrErrorTriDiffFR,Sam2ShuffleCorrErrorTriDiffFR,1,1000);

%% Figure
figure;
% plot real result
[YMin1,YMax1] = PlotMeanAndSEM(Samp1CorrTriMeanFR_New,X,[228 0 127]/255,[-3 6.4],'-',1);
[YMin3,YMax3] = PlotMeanAndSEM(Samp2ErrorTriMeanFR_New,X,[243 151 0]/255,[-3 6.4],'-',1);
[YMin2,YMax2] = PlotMeanAndSEM(Samp2CorrTriMeanFR_New,X,[26 72 157]/255,[-3 6.4],'-',1);
% maximal and minimal FR
YMax = max([YMax1 YMax2 YMax3]);
YMin = min([YMin1 YMin2 YMin3]);
hold on
% label ID of bins with significant difference of FR between FA and CR trials
for i = 1:length(Samp2CorrErrorTriDiffFRIsSig)
    if Samp2CorrErrorTriDiffFRIsSig(i) == 1
        plot([X(i)-1/TimeGain/2 X(i)+1/TimeGain/2],YMax*1.05*ones(1,2),'-','color',[0 0 0],'linewidth',2)
    end
end
hold on
for i = 1:length(Samp2CorrErrorTriDiffFRIsSig_BelowChance)
    if Samp2CorrErrorTriDiffFRIsSig_BelowChance(i) == 1
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
CorrErrTriActivity = [{Samp1CorrTriMeanFR_New};{Samp2CorrTriMeanFR_New};{Samp2ErrorTriMeanFR_New}];
