%% Plot FC density of five types of pairs.

clear; clc; close all;

%% Load result of FC density
file = dir('*SessionBased FC density within transient or sustained neurons*.mat');
for i = 1:size(file,1)
    filename = file(i,1).name; 
    if ~isempty(regexp(filename,'Within mPFC'))
        load(filename);
        FcDensity_mPFC_Tran = FcDensity_Tran_hit;
        FcDensity_mPFC_Sust = FcDensity_Sust_hit;
    elseif ~isempty(regexp(filename,'Within aAIC'))
        load(filename);
        FcDensity_aAIC_Tran = FcDensity_Tran_hit;
        FcDensity_aAIC_Sust = FcDensity_Sust_hit;
    elseif ~isempty(regexp(filename,'mPFC-aAIC'))
        load(filename);
        FcDensity_mPFCaAIC_Tran = FcDensity_Tran_hit;
        FcDensity_mPFCaAIC_Sust = FcDensity_Sust_hit;
    elseif ~isempty(regexp(filename,'aAIC-mPFC'))
        load(filename);
        FcDensity_aAICmPFC_Tran = FcDensity_Tran_hit;
        FcDensity_aAICmPFC_Sust = FcDensity_Sust_hit;
    end
end

%% FC density within transient neurons
% within mPFC
AverFcDensity_mPFC_Tran = cellfun(@mean,FcDensity_mPFC_Tran,'UniformOutput',1); 
SemFcDensity_mPFC_Tran = cellfun(@(x) std(x)/sqrt(size(x,1)),FcDensity_mPFC_Tran,'UniformOutput',1);
% within aAIC
AverFcDensity_aAIC_Tran = cellfun(@mean,FcDensity_aAIC_Tran,'UniformOutput',1); 
SemFcDensity_aAIC_Tran = cellfun(@(x) std(x)/sqrt(size(x,1)),FcDensity_aAIC_Tran,'UniformOutput',1);
% mPFC-aAIC
AverFcDensity_mPFCaAIC_Tran = cellfun(@mean,FcDensity_mPFCaAIC_Tran,'UniformOutput',1); 
SemFcDensity_mPFCaAIC_Tran = cellfun(@(x) std(x)/sqrt(size(x,1)),FcDensity_mPFCaAIC_Tran,'UniformOutput',1);
% aAIC-mPFC
AverFcDensity_aAICmPFC_Tran = cellfun(@mean,FcDensity_aAICmPFC_Tran,'UniformOutput',1); 
SemFcDensity_aAICmPFC_Tran = cellfun(@(x) std(x)/sqrt(size(x,1)),FcDensity_aAICmPFC_Tran,'UniformOutput',1);

%% Compare FC density within transient neurons between mPFC, aAIC, and mPFC-aAIC
% Two-way ANOVA
data = vertcat(FcDensity_mPFC_Tran,FcDensity_aAIC_Tran,FcDensity_mPFCaAIC_Tran);
num = cellfun(@numel,data,'UniformOutput',1);
Data = [];
RegionID = [];
PairID = [];
for i = 1:size(num,1)
    for j = 1:size(num,2)
        Data = [Data; data{i,j}];
        RegionID = [RegionID i*ones(1,num(i,j))];
        PairID = [PairID j*ones(1,num(i,j))]
    end
end
[p,tbl,stats] = anovan(Data,{RegionID PairID},'model','interaction','varnames',{'RegionID','PairID'});
save('FC density commparison_within transient neurons_between regions_mPFC_aAIC_mPFCaAIC','p','tbl','stats','-v7.3'); close;
% bar plot
figure('OuterPosition',[219 303 420 534]);
bar(1.1:3.5:1.1+3*3.5,AverFcDensity_mPFC_Tran,0.25,'facecolor','k','edgecolor','k');
hold on
errorbar(1.1:3.5:1.1+3*3.5,AverFcDensity_mPFC_Tran,SemFcDensity_mPFC_Tran,'k','linestyle','none','marker','none','capsize',10);
hold on
bar(2.1:3.5:2.1+3*3.5,AverFcDensity_aAIC_Tran,0.25,'facecolor','w','edgecolor','k');
hold on
errorbar(2.1:3.5:2.1+3*3.5,AverFcDensity_aAIC_Tran,SemFcDensity_aAIC_Tran,'k','linestyle','none','marker','none','capsize',10);
hold on
bar(3.1:3.5:3.1+3*3.5,AverFcDensity_mPFCaAIC_Tran,0.25,'facecolor','r','edgecolor','r');
hold on
errorbar(3.1:3.5:3.1+3*3.5,AverFcDensity_mPFCaAIC_Tran,SemFcDensity_mPFCaAIC_Tran,'r','linestyle','none','marker','none','capsize',10);
hold on
set(gca,'XTick',zeros(1,0),'xlim',[0 14.7]);
set(gca,'YTick',0:0.05:0.2,'YTickLabel',num2cell(0:5:20),'FontName','Arial','FontSize',16,'ylim',[0 0.2]);
ylabel('FC density (%)','FontSize',18,'FontName','Arial'); box off;
set(gcf,'Renderer','Painter'); saveas(gcf,'FC density commparison_within transient neurons_between regions_mPFC_aAIC_mPFCaAIC','fig'); close;

%% FC density within sustained neurons
% within mPFC
AverFcDensity_mPFC_Sust = cellfun(@mean,FcDensity_mPFC_Sust,'UniformOutput',1); 
SemFcDensity_mPFC_Sust = cellfun(@(x) std(x)/sqrt(size(x,1)),FcDensity_mPFC_Sust,'UniformOutput',1);
% within aAIC
AverFcDensity_aAIC_Sust = cellfun(@mean,FcDensity_aAIC_Sust,'UniformOutput',1); 
SemFcDensity_aAIC_Sust = cellfun(@(x) std(x)/sqrt(size(x,1)),FcDensity_aAIC_Sust,'UniformOutput',1);
% mPFC-aAIC
AverFcDensity_mPFCaAIC_Sust = cellfun(@mean,FcDensity_mPFCaAIC_Sust,'UniformOutput',1); 
SemFcDensity_mPFCaAIC_Sust = cellfun(@(x) std(x)/sqrt(size(x,1)),FcDensity_mPFCaAIC_Sust,'UniformOutput',1);
% aAIC-mPFC
AverFcDensity_aAICmPFC_Sust = cellfun(@mean,FcDensity_aAICmPFC_Sust,'UniformOutput',1); 
SemFcDensity_aAICmPFC_Sust = cellfun(@(x) std(x)/sqrt(size(x,1)),FcDensity_aAICmPFC_Sust,'UniformOutput',1);

%% Compare FC density within sustained neurons between mPFC, aAIC, and mPFC-aAIC
% Two-way ANOVA
data = vertcat(FcDensity_mPFC_Sust,FcDensity_aAIC_Sust,FcDensity_mPFCaAIC_Sust);
num = cellfun(@numel,data,'UniformOutput',1);
Data = [];
RegionID = [];
PairID = [];
for i = 1:size(num,1)
    for j = 1:size(num,2)
        Data = [Data; data{i,j}];
        RegionID = [RegionID i*ones(1,num(i,j))];
        PairID = [PairID j*ones(1,num(i,j))]
    end
end
[p,tbl,stats] = anovan(Data,{RegionID PairID},'model','interaction','varnames',{'RegionID','PairID'});
save('FC density commparison_within sustained neurons_between regions_mPFC_aAIC_mPFCaAIC','p','tbl','stats','-v7.3'); close;
% bar plot
figure('OuterPosition',[219 303 420 534]);
bar(1.1:3.5:1.1+3.5,AverFcDensity_mPFC_Sust,0.25,'facecolor','k','edgecolor','k');
hold on
errorbar(1.1:3.5:1.1+3.5,AverFcDensity_mPFC_Sust,SemFcDensity_mPFC_Sust,'k','linestyle','none','marker','none','capsize',10);
hold on
bar(2.1:3.5:2.1+3.5,AverFcDensity_aAIC_Sust,0.25,'facecolor','w','edgecolor','k');
hold on
errorbar(2.1:3.5:2.1+3.5,AverFcDensity_aAIC_Sust,SemFcDensity_aAIC_Sust,'k','linestyle','none','marker','none','capsize',10);
hold on
bar(3.1:3.5:3.1+3.5,AverFcDensity_mPFCaAIC_Sust,0.25,'facecolor','r','edgecolor','r');
hold on
errorbar(3.1:3.5:3.1+3.5,AverFcDensity_mPFCaAIC_Sust,SemFcDensity_mPFCaAIC_Sust,'r','linestyle','none','marker','none','capsize',10);
hold on
set(gca,'XTick',zeros(1,0),'xlim',[0 7.7]);
set(gca,'YTick',0:0.1:0.5,'YTickLabel',num2cell(0:10:50),'FontName','Arial','FontSize',16,'ylim',[0 0.5]);
ylabel('FC density (%)','FontSize',18,'FontName','Arial'); box off;
set(gcf,'Renderer','Painter'); saveas(gcf,'FC density commparison_within sustained neurons_between regions_mPFC_aAIC_mPFCaAIC','fig'); close;




