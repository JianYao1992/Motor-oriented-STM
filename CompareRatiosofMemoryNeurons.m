%% Proportion of memory, transient, and sustained neurons, following projection-specific manipulation.

clear; clc; close all;

%% Assignment
TargetBrain = 'aAIC'; 
IsBetweenInhibitionAndCtrl = 0;
switch IsBetweenInhibitionAndCtrl
    case 1
        chara = 'NpHRGroup';
        C = [0 0.5 0];
    case 0
        chara = 'ChR2Group';
        C = [67 106 178]/255;
end

%% Result of control and experimental group
% control group
load([TargetBrain 'SelectivityData_CtrlGroup.mat']);
SortedID_Ctrl = SortedID;
SustUnitID_Ctrl = SustainedUnitID; SwitSustUnitID_Ctrl = SwitchedSustainedUnitID; TranUnitID_Ctrl = TransientCodingNeuronID; 
PropOfMemUnits_Ctrl = (length(SustUnitID_Ctrl)+length(SwitSustUnitID_Ctrl)+length(TranUnitID_Ctrl))/length(SortedID_Ctrl);
PropOfSustUnits_Ctrl = length(SustUnitID_Ctrl)/length(SortedID_Ctrl);
PropOfTranUnits_Ctrl = (length(SwitSustUnitID_Ctrl)+length(TranUnitID_Ctrl))/length(SortedID_Ctrl);
UnitIDWithCoding_1s_Ctrl = UnitIDWith1sCoding; UnitIDWithCoding_2s_Ctrl = UnitIDWith2sCoding; UnitIDWithCoding_3s_Ctrl = UnitIDWith3sCoding; UnitIDWithSwitCoding_Ctrl = [SwitchedSustainedUnitID' Switched3sCodingUnitID Switched2sCodingUnitID'];
PropOfUnitsWithCoding_1s_Ctrl = length(UnitIDWithCoding_1s_Ctrl)/length(SortedID_Ctrl);
PropOfUnitsWithCoding_2s_Ctrl = length(UnitIDWithCoding_2s_Ctrl)/length(SortedID_Ctrl);
PropOfUnitsWithCoding_3s_Ctrl = length(UnitIDWithCoding_3s_Ctrl)/length(SortedID_Ctrl);
PropOfUnitsWithSwitCoding_Ctrl = length(UnitIDWithSwitCoding_Ctrl)/length(SortedID_Ctrl);
DelayBinCodingUnitID_Ctrl = DelayBinCodingUnitID; 
% Experimental group
load([TargetBrain 'SelectivityData_' chara '.mat']);
SortedID_Exp = SortedID;
SustUnitID_Exp = SustainedUnitID; SwitSustUnitID_Exp = SwitchedSustainedUnitID; TranUnitID_Exp = TransientCodingNeuronID;
PropOfMemUnits_Exp = (length(SustUnitID_Exp)+length(SwitSustUnitID_Exp)+length(TranUnitID_Exp))/length(SortedID_Exp);
PropOfSustUnits_Exp = length(SustUnitID_Exp)/length(SortedID_Exp);
PropOfTranUnits_Exp = (length(SwitSustUnitID_Exp)+length(TranUnitID_Exp))/length(SortedID_Exp);
UnitIDWithCoding_1s_Exp = UnitIDWith1sCoding; UnitIDWithCoding_2s_Exp = UnitIDWith2sCoding; UnitIDWithCoding_3s_Exp = UnitIDWith3sCoding;  UnitIDWithSwitCoding_Exp = [SwitchedSustainedUnitID' Switched3sCodingUnitID Switched2sCodingUnitID'];
PropOfUnitsWithCoding_1s_Exp = length(UnitIDWithCoding_1s_Exp)/length(SortedID_Exp);
PropOfUnitsWithCoding_2s_Exp = length(UnitIDWithCoding_2s_Exp)/length(SortedID_Exp);
PropOfUnitsWithCoding_3s_Exp = length(UnitIDWithCoding_3s_Exp)/length(SortedID_Exp);
PropOfUnitsWithSwitCoding_Exp = length(UnitIDWithSwitCoding_Exp)/length(SortedID_Exp);
DelayBinCodingUnitID_Exp = DelayBinCodingUnitID;  

%% Plot proportion of memory neurons
figure('OuterPosition',[219 303 420 534]);
bar(1.1,PropOfMemUnits_Ctrl,0.6,'k','edgecolor','none');
hold on
bar(1.9,PropOfMemUnits_Exp,0.6,'facecolor',C,'edgecolor','none');
set(gca,'XTick',zeros(1,0),'xlim',[0.6 2.4]);
set(gca,'YTick',0:0.2:0.8,'YTickLabel',{'0','20','40','60','80'},'FontName','Arial','FontSize',16,'ylim',[0 0.8]);
ylabel('Proportion (%)','FontSize',18,'FontName','Arial');
box off;
[~,p,~,~] = prop_test([length(SustUnitID_Ctrl)+length(SwitSustUnitID_Ctrl)+length(TranUnitID_Ctrl) length(SustUnitID_Exp)+length(SwitSustUnitID_Exp)+length(TranUnitID_Exp)],[length(SortedID_Ctrl) length(SortedID_Exp)],false);
title(sprintf('p=%d',p));
set(gcf,'Renderer','Painter'); saveas(gcf,['CompareMemory' TargetBrain 'NeuronsBetweenCtrland' chara],'fig'); close;

%% Plot proportion of transient and sustained neurons
figure('OuterPosition',[219 303 420 534]);
bar(1.1,PropOfSustUnits_Ctrl,0.6,'k','edgecolor','none');
hold on
bar(1.9,PropOfSustUnits_Exp,0.6,'facecolor',C,'edgecolor','none');
hold on
bar(2.7,PropOfTranUnits_Ctrl,0.6,'k','edgecolor','none');
hold on
bar(3.5,PropOfTranUnits_Exp,0.6,'facecolor',C,'edgecolor','none');
hold on
set(gca,'XTick',zeros(1,0),'xlim',[0.6 4]);
set(gca,'YTick',0:0.2:0.8,'YTickLabel',{'0','20','40','60','80'},'FontName','Arial','FontSize',16,'ylim',[0 0.8]);
ylabel('Proportion ( % )','FontSize',18,'FontName','Arial');
box off;
[~,p1,~,~] = prop_test([length(SwitSustUnitID_Ctrl)+length(TranUnitID_Ctrl) length(SwitSustUnitID_Exp)+length(TranUnitID_Exp)],[length(SortedID_Ctrl) length(SortedID_Exp)],false);
[~,p2,~,~] = prop_test([length(SustUnitID_Ctrl) length(SustUnitID_Exp)],[length(SortedID_Ctrl) length(SortedID_Exp)],false);
title(sprintf('Ptran=%d, Psust=%d',p1,p2));
set(gcf,'Renderer','Painter'); saveas(gcf,['CompareTranAndSust' TargetBrain 'NeuronsBetweenCtrland' chara],'fig'); close;

%% Plot proportion of neurons with 1-second STM-encoding ability 
figure('OuterPosition',[219 303 420 534]);
bar(1.1,PropOfUnitsWithCoding_1s_Ctrl,0.6,'k','edgecolor','none');
hold on
bar(1.9,PropOfUnitsWithCoding_1s_Exp,0.6,'facecolor',C,'edgecolor','none');
set(gca,'XTick',zeros(1,0),'xlim',[0.6 2.4]);
set(gca,'YTick',0:0.05:0.25,'YTickLabel',{'0','5','10','15','20','25'},'FontName','Arial','FontSize',16,'ylim',[0 0.25]);
ylabel('Proportion (%)','FontSize',18,'FontName','Arial');
box off;
[~,p,~,~] = prop_test([length(UnitIDWithCoding_1s_Ctrl) length(UnitIDWithCoding_1s_Exp)],[length(SortedID_Ctrl) length(SortedID_Exp)],false);
title(sprintf('p=%d',p));
set(gcf,'Renderer','Painter'); saveas(gcf,['Compare' TargetBrain 'NeuronsWith1sCodingBetweenCtrland' chara],'fig'); close all;

%% Plot proportion of neurons with 2-second STM-encoding ability
figure('OuterPosition',[219 303 420 534]);
bar(1.1,PropOfUnitsWithCoding_2s_Ctrl,0.6,'k','edgecolor','none');
hold on
bar(1.9,PropOfUnitsWithCoding_2s_Exp,0.6,'facecolor',C,'edgecolor','none');
set(gca,'XTick',zeros(1,0),'xlim',[0.6 2.4]);
set(gca,'YTick',0:0.05:0.2,'YTickLabel',{'0','5','10','15','20'},'FontName','Arial','FontSize',16,'ylim',[0 0.2]);
ylabel('Proportion (%)','FontSize',18,'FontName','Arial');
box off;
[~,p,~,~] = prop_test([length(UnitIDWithCoding_2s_Ctrl) length(UnitIDWithCoding_2s_Exp)],[length(SortedID_Ctrl) length(SortedID_Exp)],false);
title(sprintf('p=%d',p));
set(gcf,'Renderer','Painter'); saveas(gcf,['Compare' TargetBrain 'NeuronsWith2sCodingBetweenCtrland' chara],'fig'); close;

%% Plot proportion of neurons with 3-second STM-encoding ability
figure('OuterPosition',[219 303 420 534]);
bar(1.1,PropOfUnitsWithCoding_3s_Ctrl,0.6,'k','edgecolor','none');
hold on
bar(1.9,PropOfUnitsWithCoding_3s_Exp,0.6,'facecolor',C,'edgecolor','none');
set(gca,'XTick',zeros(1,0),'xlim',[0.6 2.4]);
set(gca,'YTick',0:0.02:0.10,'YTickLabel',{'0','2','4','6','8','10'},'FontName','Arial','FontSize',16,'ylim',[0 0.1]);
ylabel('Proportion (%)','FontSize',18,'FontName','Arial');
box off;
[~,p,~,~] = prop_test([length(UnitIDWithCoding_3s_Ctrl) length(UnitIDWithCoding_3s_Exp)],[length(SortedID_Ctrl) length(SortedID_Exp)],false);
title(sprintf('p=%d',p));
set(gcf,'Renderer','Painter'); saveas(gcf,['Compare' TargetBrain 'NeuronsWith3sCodingBetweenCtrland' chara],'fig'); close;

%% Plot proportion of neurons with switched coding
figure('OuterPosition',[219 303 420 534]);
bar(1.1,PropOfUnitsWithSwitCoding_Ctrl,0.6,'k','edgecolor','none');
hold on
bar(1.9,PropOfUnitsWithSwitCoding_Exp,0.6,'facecolor',C,'edgecolor','none');
set(gca,'XTick',zeros(1,0),'xlim',[0.6 2.4]);
set(gca,'YTick',0:0.05:0.15,'YTickLabel',{'0','5','10','15'},'FontName','Arial','FontSize',16,'ylim',[0 0.15]);
ylabel('Proportion selective neurons ( % )','FontSize',18,'FontName','Arial');
box off;
[~,p,~,~] = prop_test([length(UnitIDWithSwitCoding_Ctrl) length(UnitIDWithSwitCoding_Exp)],[length(SortedID_Ctrl) length(SortedID_Exp)],false);
title(sprintf('p=%d',p));
set(gcf,'Renderer','Painter'); saveas(gcf,['Compare' TargetBrain 'NeuronsWithSwitchedCodingBetweenCtrland' chara],'fig'); close;

%% Plot proportion of neurons over the time course of delay
RatioOfCodingUnitsInSpecificBin_Ctrl = [];
RatioOfCodingUnitsInSpecificBin_Exp = [];
p = [];
for iBin = 1:length(DelayBinCodingUnitID_Ctrl)
    RatioOfCodingUnitsInSpecificBin_Ctrl = [RatioOfCodingUnitsInSpecificBin_Ctrl length(DelayBinCodingUnitID_Ctrl{iBin})/length(SortedID_Ctrl)];
    RatioOfCodingUnitsInSpecificBin_Exp = [RatioOfCodingUnitsInSpecificBin_Exp length(DelayBinCodingUnitID_Exp{iBin})/length(SortedID_Exp)];
    [~,p(iBin), ~,~] = prop_test([length(DelayBinCodingUnitID_Ctrl{iBin}) length(DelayBinCodingUnitID_Exp{iBin})],[length(SortedID_Ctrl) length(SortedID_Exp)],false);
end
p = p*length(DelayBinCodingUnitID_Ctrl);
figure;
plot([1.5 2.5 3.5 4.5],RatioOfCodingUnitsInSpecificBin_Ctrl,'k','linewidth',3,'marker','o','markersize',10,'markerfacecolor','k',...
    'markeredgecolor','none');
hold on
plot([1.5 2.5 3.5 4.5],RatioOfCodingUnitsInSpecificBin_Exp,'color',C,'linewidth',3,'marker','o','markersize',10,'markerfacecolor',C,...
    'markeredgecolor','none');
hold on
plot([0 0],[0 0.5],'k','linestyle','- -','linewidth',3);
plot([1 1],[0 0.5],'k','linestyle','- -','linewidth',3);
plot([5 5],[0 0.5],'k','linestyle','- -','linewidth',3);
plot([5.5 5.5],[0 0.5],'k','linestyle','- -','linewidth',3);
set(gca,'XTick',0:2:6,'XTickLabel',{'0','2','4','6'},'xlim',[-0.5 6],'FontName','Arial','FontSize',16);
set(gca,'YTick',0:0.1:0.5,'YTickLabel',{'0','10','20','30','40','50'},'FontName','Arial','FontSize',16,'ylim',[0 0.5]);
xlabel('Time from sample onset ( s )','FontSize',18,'FontName','Arial');
ylabel('Proportion (%)','FontSize',18,'FontName','Arial');
box off;
title(sprintf('p(1)=%d p(2)=%d p(3)=%d p(4)=%d',p(1),p(2),p(3),p(4)));
set(gcf,'Renderer','Painter'); saveas(gcf,['Compare' TargetBrain 'NeuronsWithCodinginSpecificBinBetweenCtrland' chara],'fig'); close all;
