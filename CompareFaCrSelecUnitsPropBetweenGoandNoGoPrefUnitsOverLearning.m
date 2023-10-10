%% Compare proportion of neurons showing selectivity to FA and CR between Go- and NoGo-preferred neurons.

clear; clc; close all;

%% Assignment
LearningDayNum = 6;
C1 = [1 0 0];
C2 = [0 0 1];

%% load result
load('Proportion of performance dependent neurons in memory neurons over learning_RanksumTest-CtrlGroup.mat')

%% Proportion of neurons showing selectivity to FA and CR in memory neurons
% mPFC Go-preferred neurons, mPFC NoGo-preferred neurons, aAIC Go-preferred neurons, and aAIC NoGo-preferred neurons
PropPerfDependUnits_GoPrefer_mPFC = cell(1,LearningDayNum);
PropPerfDependUnits_NoGoPrefer_mPFC = cell(1,LearningDayNum);
PropPerfDependUnits_GoPrefer_aAIC = cell(1,LearningDayNum);
PropPerfDependUnits_NoGoPrefer_aAIC = cell(1,LearningDayNum);
for iDay = 1:LearningDayNum
    if ~isempty(IsPerfRelatedUnit.mPFC.GoPreferredUnits{iDay}{1})
        PropPerfDependUnits_GoPrefer_mPFC{1,iDay} = mean(IsPerfRelatedUnit.mPFC.GoPreferredUnits{iDay}{1}(2,:));
    else
        PropPerfDependUnits_GoPrefer_mPFC{1,iDay} = NaN;
    end
    if ~isempty(IsPerfRelatedUnit.mPFC.NoGoPreferredUnits{iDay}{1})
        PropPerfDependUnits_NoGoPrefer_mPFC{1,iDay} = mean(IsPerfRelatedUnit.mPFC.NoGoPreferredUnits{iDay}{1}(2,:));
    else
        PropPerfDependUnits_NoGoPrefer_mPFC{1,iDay} = NaN;
    end
    if ~isempty(IsPerfRelatedUnit.aAIC.GoPreferredUnits{iDay}{1})
        PropPerfDependUnits_GoPrefer_aAIC{1,iDay} = mean(IsPerfRelatedUnit.aAIC.GoPreferredUnits{iDay}{1}(2,:));
    else
        PropPerfDependUnits_GoPrefer_aAIC{1,iDay} = NaN;
    end
    if ~isempty(IsPerfRelatedUnit.aAIC.NoGoPreferredUnits{iDay}{1})
        PropPerfDependUnits_NoGoPrefer_aAIC{1,iDay} = mean(IsPerfRelatedUnit.aAIC.NoGoPreferredUnits{iDay}{1}(2,:));
    else
        PropPerfDependUnits_NoGoPrefer_aAIC{1,iDay} = NaN;
    end
end

%% Compare proportion between Go- and NoGo-preferred neurons
% mPFC
figure('position',[200 200 400 600]);
plot(1:1:LearningDayNum,cell2mat(PropPerfDependUnits_GoPrefer_mPFC),'color',C1,'marker','o','markerfacecolor',C1,'markeredgecolor','none'); hold on
plot(1:1:LearningDayNum,cell2mat(PropPerfDependUnits_NoGoPrefer_mPFC),'color',C2,'marker','o','markerfacecolor',C2,'markeredgecolor','none'); hold on
set(gca,'XTick',1:1:LearningDayNum,'XTickLabel',num2cell(1:1:LearningDayNum),'XLim',[0.5 LearningDayNum+0.5],'YTIck',0:0.2:1,'YTickLabel',num2cell(0:20:100),'YLim',[0 1]);
xlabel('Learning day');
ylabel('Porportion of performance-dependent neurons (%)');
box off
set(gcf,'Render','Painter'); saveas(gcf,'Proportion of performance-dependent neurons in memory neurons over learning-mPFC','fig'); close all;
[p,tb1,stats] = unbalanced_anova_test(PropPerfDependUnits_GoPrefer_mPFC,PropPerfDependUnits_NoGoPrefer_mPFC); % Two-way ANOVA
save('Statistics for proportion of performance-dependent neurons in memory neurons over learning-mPFC','p','tb1','stats','-v7.3');
% aAIC
figure('position',[200 200 400 600]);
plot(1:1:LearningDayNum,cell2mat(PropPerfDependUnits_GoPrefer_aAIC),'color',C1,'marker','o','markerfacecolor',C1,'markeredgecolor','none'); hold on
plot(1:1:LearningDayNum,cell2mat(PropPerfDependUnits_NoGoPrefer_aAIC),'color',C2,'marker','o','markerfacecolor',C2,'markeredgecolor','none'); hold on
set(gca,'XTick',1:1:LearningDayNum,'XTickLabel',num2cell(1:1:LearningDayNum),'XLim',[0.5 LearningDayNum+0.5],'YTIck',0:0.2:1,'YTickLabel',num2cell(0:20:100),'YLim',[0 1]);
xlabel('Learning day');
ylabel('Porportion of performance-dependent neurons (%)');
box off
set(gcf,'Render','Painter'); saveas(gcf,'Proportion of performance-dependent neurons in memory neurons over learning-aAIC','fig'); close all;
[p,tb1,stats] = unbalanced_anova_test(PropPerfDependUnits_GoPrefer_aAIC,PropPerfDependUnits_NoGoPrefer_aAIC); % Two-way ANOVA
save('Statistics for proportion of performance-dependent neurons in memory neurons over learning-aAIC','p','tb1','stats','-v7.3');



