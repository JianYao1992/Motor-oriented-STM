% Receiver operating characteristic(ROC) analysis for sustained and transient neurons

clear; clc; close all;

%% Assignment
TargetBrain = 'aAIC';
UnitsType = 1; % 1: only sustained; 2: only transient; 3: sustained and transient
IsDVorFRROC = 1; % 1: decision variable; 2: raw FR 
IsBetweenInhibitionAndCtrl = 0; 
TimeGain = 10;
switch IsBetweenInhibitionAndCtrl
    case 1
        Group = 'NpHRGroup';
        C = [0 125 0]/255;
    case 0
        Group = 'ChR2Group';
        C = [67 106 178]/255;
end

%% Load result of control group
load([TargetBrain 'FRinS1S2_CtrlGroup']); % FR
load([TargetBrain 'SelectivityData_CtrlGroup']); % ID of sustained and transient neurons
period = 1:size(TargetBrainUnitsFRinS1{1,1},2); 

%% ID of target neurons of control group
switch UnitsType
    case 1
        TarUnitsID = SustainedUnitID(:);
        savefilename = 'Sustained units ROC';
    case 2
        TarUnitsID = [TransientCodingNeuronID(:); SwitchedSustainedUnitID(:)];
        savefilename = 'Transient units ROC';
    case 3
        TarUnitsID = [SustainedUnitID(:); TransientCodingNeuronID(:); SwitchedSustainedUnitID(:)];
        savefilename = 'Memory units ROC';  
end

%% auROC of neurons of control group
TarUnitsAUC_Ctrl = zeros(length(TarUnitsID),floor(length(period)/TimeGain));
for iUnit = 1:length(TarUnitsID) 
    if mod(iUnit,100) == 0 || iUnit == 1 || iUnit == length(TarUnitsID)
        disp(['Current Neuron number=' num2str(iUnit)]);
    end
    for iSecBin = 1:floor(length(period)/TimeGain)
        % step 1: binned FR
        tempEventPeriod = (iSecBin-1)*TimeGain+1:iSecBin*TimeGain;
        AverFR_S1 = mean(TargetBrainUnitsFRinS1{TarUnitsID(iUnit)}(:,tempEventPeriod),2);
        AverFR_S2 = mean(TargetBrainUnitsFRinS2{TarUnitsID(iUnit)}(:,tempEventPeriod),2);
        if IsDVorFRROC == 1
            % step 2: DV of S1 and S2 trials
            [DV_S1,DV_S2] = DecisionVariableCalculation(AverFR_S1,AverFR_S2);
            % step 3: TPR and FPR
            [TPR,FPR] = TprFprCalculation(DV_S1,DV_S2);
        else
            % step 3: TPR and FPR
            [TPR,FPR] = TprFprCalculation(AverFR_S1,AverFR_S2);
        end
        AUC = AucAnalysis(FPR,TPR);
        TarUnitsAUC_Ctrl(iUnit,iSecBin) = AUC;
    end
end

%% Load results of experimental group
load([TargetBrain 'FRinS1S2_' Group '.mat']);
load([TargetBrain 'SelectivityData_' Group]);

%% ID of neurons of experimental group
switch UnitsType
    case 1
        TarUnitsID = SustainedUnitID(:);
    case 2
        TarUnitsID = [TransientCodingNeuronID(:); SwitchedSustainedUnitID(:)];
    case 3
        TarUnitsID = [SustainedUnitID(:); TransientCodingNeuronID(:); SwitchedSustainedUnitID(:)];
end

%% auROC of neurons of experimental group
TarUnitsAUC_Exp = zeros(length(TarUnitsID),floor(length(period)/TimeGain));
for iUnit = 1:length(TarUnitsID)
    if mod(iUnit,100) == 0 || iUnit == 1 || iUnit == length(TarUnitsID)
        disp(['Current Neuron number=' num2str(iUnit)])
    end
    for iSecBin = 1:floor(length(period)/TimeGain)
        % step 1: binned FR
        tempEventPeriod = (iSecBin-1)*TimeGain + 1:iSecBin*TimeGain;
        AverFR_S1 = mean(TargetBrainUnitsFRinS1{TarUnitsID(iUnit)}(:,tempEventPeriod),2);
        AverFR_S2 = mean(TargetBrainUnitsFRinS2{TarUnitsID(iUnit)}(:,tempEventPeriod),2);
        if IsDVorFRROC==1
            % step 2: DV of S1 and S2 trials
            [DV_S1,DV_S2] = DecisionVariableCalculation(AverFR_S1,AverFR_S2);
            % step 3: TPR and FPR
            [TPR,FPR] = TprFprCalculation(DV_S1,DV_S2);
        else
            % step 3: TPR and FPR
            [TPR,FPR] = TprFprCalculation(AverFR_S1,AverFR_S2);
        end
        % step 4: auROC
        AUC = AucAnalysis(FPR,TPR);
        TarUnitsAUC_Exp(iUnit,iSecBin) = AUC;
    end
end
TarUnitsAUC_Ctrl(:,1) = []; 
TarUnitsAUC_Ctrl(:,end) = [];
TarUnitsAUC_Exp(:,1) = [];
TarUnitsAUC_Exp(:,end) = [];

%% Plot auROC of neurons, following projection-specific manipulation
figure;
errorbar(-0.5:1:size(TarUnitsAUC_Ctrl,2)-1.5,mean(TarUnitsAUC_Ctrl,1),std(TarUnitsAUC_Ctrl,0,1)/sqrt(size(TarUnitsAUC_Ctrl,1)),'color','k','marker','o','markerfacecolor',[0 0 0],'linewidth',3);
hold on
errorbar(-0.5:size(TarUnitsAUC_Exp,2)-1.5,mean(TarUnitsAUC_Exp,1),std(TarUnitsAUC_Exp,0,1)/sqrt(size(TarUnitsAUC_Exp,1)),'color',C,'marker','o','markerfacecolor',C,'linewidth',3);
hold on
plot([0 0],[0 1],'k','linestyle','- -','linewidth',3); hold on
plot([1 1],[0 1],'k','linestyle','- -','linewidth',3); hold on
plot([5 5],[0 1],'k','linestyle','- -','linewidth',3); hold on
plot([5.5 5.5],[0 1],'k','linestyle','- -','linewidth',3); hold on
set(gca,'XTick',0:2:6,'YTickLabel',{'0','2','4','6'},'xlim',[-1 6],'FontName','Arial','FontSize',16);
set(gca,'YTick',0.5:0.1:0.9,'YTickLabel',{'0.5','0.6','0.7','0.8','0.9'},'FontName','Arial','FontSize',16,'ylim',[0.5 0.8]);
xlabel('Time from sample onset ( s )','FontSize',18,'FontName','Arial');
ylabel('Averaged delay auROC','FontSize',18,'FontName','Arial');
box off;

%% Tw-ANOVA-md
between_group = [];
Unit_ID = [];
within_group = [];
performance = [];
TarUnitsAUC_Exp = TarUnitsAUC_Exp(:,3:6); % delay period
TarUnitsAUC_Ctrl = TarUnitsAUC_Ctrl(:,3:6);
% experimental group
for i = 1:size(TarUnitsAUC_Exp,1) % neuron
    BinNum = size(TarUnitsAUC_Exp,2);
    between_group = [between_group; 1*ones(BinNum,1)];
    Unit_ID = [Unit_ID; i*ones(BinNum,1)];
    for j = 1:BinNum % second
        within_group = [within_group; j];
        performance = [performance;TarUnitsAUC_Exp(i,j)];
    end
end
% control group
for i = 1:size(TarUnitsAUC_Ctrl,1) % neuron
    between_group = [between_group; 2*ones(BinNum,1)];
    Unit_ID = [Unit_ID; (i+size(TarUnitsAUC_Exp,1))*ones(BinNum,1)];
    for j = 1:BinNum % second
        within_group = [within_group; j];
        performance = [performance; TarUnitsAUC_Ctrl(i,j)];
    end
end
Frame = [performance between_group within_group Unit_ID];
[SSQs, DFs, MSQs, Fs, Ps] = mixed_between_within_anova(Frame,0,strcat(savefilename,Group));
title(sprintf('p=%d',Ps{1,1}));
set(gcf,'Renderer','Painter'); saveas(gcf,[savefilename TargetBrain 'BetweenCtrland' Group],'fig'); close;







