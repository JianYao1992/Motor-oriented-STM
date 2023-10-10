%% FR difference between different samples of specific type of FC pairs, following projection-specific manipulation

clear; clc; close all;

%% Assignment
Target = 'mPFC-aAIC';
FcPairsID = 4;
PairsType = {'Con Act','Con Inact','Incon Act','Incon Inact','Non-memory'};

%% Load result of FR difference
% optogenetic suppression group
load(sprintf('Functional Coupling Firing Rate Difference_5pairs_%s_NpHRGroup',Target));
data_inhibition = DiffinFR{FcPairsID};
% control group
load(sprintf('Functional Coupling Firing Rate Difference_5pairs_%s_CtrlGroup',Target));
data_Ctrl = DiffinFR{FcPairsID};
% optogenetic activation group
load(sprintf('Functional Coupling Firing Rate Difference_5pairs_%s_ChR2Group',Target));
data_activation = DiffinFR{FcPairsID};

%% Error bar plot
figure('OuterPosition',[219 150 420 350]);
AverFrDiff_inhibition = cellfun(@(x) mean(x),data_inhibition,'UniformOutput',true);
StdFRDiff_inhibition = cellfun(@(x) std(x),data_inhibition,'UniformOutput',true);
n_inhibition = cellfun(@(x) size(x,1),data_inhibition,'UniformOutput',true);
AverFrDiff_Ctrl = cellfun(@(x) mean(x),data_Ctrl,'UniformOutput',true);
StdFrDiff_Ctrl = cellfun(@(x) std(x),data_Ctrl,'UniformOutput',true); 
n_Ctrl = cellfun(@(x) size(x,1),data_Ctrl,'UniformOutput',true);
AverFrDiff_activation = cellfun(@(x) mean(x),data_activation,'UniformOutput',true);
StdFRDiff_activation = cellfun(@(x) std(x),data_activation,'UniformOutput',true);
n_activation = cellfun(@(x) size(x,1),data_activation,'UniformOutput',true);
errorbar(1.5:1:4.5,AverFrDiff_inhibition,StdFRDiff_inhibition./sqrt(n_inhibition),'color',[0 125 0]/255,'linewidth',3,'marker','o','markerfacecolor',[0 125 0]/255,'markeredgecolor','none'); hold on
errorbar(1.5:1:4.5,AverFrDiff_Ctrl,StdFrDiff_Ctrl./sqrt(n_Ctrl),'color','k','linewidth',3,'marker','o','markerfacecolor','k','markeredgecolor','none'); hold on
errorbar(1.5:1:4.5,AverFrDiff_activation,StdFRDiff_activation./sqrt(n_activation),'color',[67 106 178]/255,'linewidth',3,'marker','o','markerfacecolor',[67 106 178]/255,'markeredgecolor','none'); hold on
MaxY = max([AverFrDiff_inhibition+StdFRDiff_inhibition./sqrt(n_inhibition) AverFrDiff_Ctrl+StdFrDiff_Ctrl./sqrt(n_Ctrl) AverFrDiff_activation+StdFRDiff_activation./sqrt(n_activation)])+1;
MinY = max([AverFrDiff_inhibition-StdFRDiff_inhibition./sqrt(n_inhibition) AverFrDiff_Ctrl-StdFrDiff_Ctrl./sqrt(n_Ctrl) AverFrDiff_activation-StdFRDiff_activation./sqrt(n_activation)])-1;
plot([0 0],[MinY MaxY],'k--'); hold on
plot([1 1],[MinY MaxY],'k--'); hold on
plot([5 5],[MinY MaxY],'k--'); hold on
plot([5.5 5.5],[MinY MaxY],'k--'); hold on

%% Two-way ANOVA
[p,tb1] = unbalanced_anova_test(data_inhibition,data_Ctrl,data_activation);
title(['F(' num2str(tb1{2,3}) ')=' num2str(tb1{2,6}) 'P(1) = ' num2str(tb1{2,7}) 'F(' num2str(tb1{3,3}) ')=' num2str(tb1{3,6}) 'P(2) = ' num2str(tb1{3,7})]);
set(gca,'XTick',0:1:6,'YTickLabel',{'0','1','2','3','4','5','6'},'xlim',[-0.5 6],'FontName','Arial','FontSize',16);
set(gca,'YTick',0:2:10,'YTickLabel',{'0','2','4','6','8','10'},'FontName','Arial','FontSize',16,'ylim',[MinY MaxY]);
xlabel('Time from sample onset (s)','FontSize',18,'FontName','Arial');
ylabel('FR difference (Hz)','FontSize',18,'FontName','Arial');
box off;
set(gcf,'Renderer','Painter'); saveas(gcf,['CompareFcFrDiff_' Target '_' PairsType{FcPairsID} '_BetweenControlandOptoGroup'],'fig'); close;

