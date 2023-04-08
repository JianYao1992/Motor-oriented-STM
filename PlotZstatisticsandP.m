%% Plot Z-statistics and P for difference in FR across bins between S1 and S2 trials for example neuron

clear; clc; close all;

%% Assignment
TargetBrain = 'aAIC'; 
UnitID = 150; 
X = [0.5 4.5];

%% Load selectivity results
load([TargetBrain 'SelectivityData_CtrlGroup']);

%% Plot Z-statistics and P
figure1 = figure('Color','w','Position',[20,20,500,500]);
% Z-statistics
axes1 = axes('Parent',figure1,'position',[0.1 0.1 0.8 0.8],'FontSize',14,'FontName','Times New Roman');
box(axes1,'on');
hold(axes1,'all');
colormap('jet');
imagesc(ZStatistics{1,UnitID},[-2 2]);
hold on
plot([0.5 0.5],[0.5 4.5],'k','LineStyle','--','linewidth',3);
plot([4.5 4.5],[0.5 4.5],'k','LineStyle','--','linewidth',3);
plot([0.5 4.5],[0.5 0.5],'k','LineStyle','--','linewidth',3);
plot([0.5 4.5],[4.5 4.5],'k','LineStyle','--','linewidth',3);
set(gca,'YDir','normal');
set(gca,'XTick',[1.5 3.5],'XTickLabel',{'2','4'},'YTick',[1.5 3.5],'YTickLabel',{'2','4'},...
    'xlim',[0.3 4.7],'ylim',[0.3 4.7],'FontName','Arial','FontSize',16);
% P
for i = 1:size(P{1,UnitID},1)
    for j = 1:size(P{1,UnitID},1)
        text(j,i,num2str(P{1,UnitID}(i,j)),'fontsize',5); hold on;
    end
end
xlabel('Reference bin ( s )','FontSize',18,'FontName','Arial');
ylabel('Target bin ( s )','FontSize',18,'FontName','Arial');
set(gcf,'Renderer','Painter'); saveas(gcf,[TargetBrain '_UnitID ' num2str(UnitID) '_ZstatisticandP'],'fig'); close all;

