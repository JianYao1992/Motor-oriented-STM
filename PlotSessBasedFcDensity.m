%% Plot FC density of five types of pairs.

clear; clc; close all;

%% Load result of FC density
file = dir('*SessionBased FC density_1msbin_5types_hit*.mat');
for i = 1:size(file,1)
    filename = file(i,1).name; 
    if ~isempty(regexp(filename,'Within mPFC'))
        data_mPFC = load(filename);
    elseif ~isempty(regexp(filename,'Within aAIC'))
        data_aAIC = load(filename);
    elseif ~isempty(regexp(filename,'mPFC-aAIC'))
        data_mPFCaAIC = load(filename);
    elseif ~isempty(regexp(filename,'aAIC-mPFC'))
        data_aAICmPFC = load(filename);
    end
end

%% FC density
% within mPFC
AverFcDensity_mPFC = cellfun(@mean,data_mPFC.FcDensity_hit,'UniformOutput',1); 
SemFcDensity_mPFC = cellfun(@(x) std(x)/sqrt(size(x,1)),data_mPFC.FcDensity_hit,'UniformOutput',1);
% within aAIC
AverFcDensity_aAIC = cellfun(@mean,data_aAIC.FcDensity_hit,'UniformOutput',1); 
SemFcDensity_aAIC = cellfun(@(x) std(x)/sqrt(size(x,1)),data_aAIC.FcDensity_hit,'UniformOutput',1);
% mPFC-aAIC
AverFcDensity_mPFCaAIC = cellfun(@mean,data_mPFCaAIC.FcDensity_hit,'UniformOutput',1); 
SemFcDensity_mPFCaAIC = cellfun(@(x) std(x)/sqrt(size(x,1)),data_mPFCaAIC.FcDensity_hit,'UniformOutput',1);
% % aAIC-mPFC
% AverFcDensity_aAICmPFC = cellfun(@mean,data_aAICmPFC.FcDensity_hit,'UniformOutput',1); 
% SemFcDensity_aAICmPFC = cellfun(@(x) std(x)/sqrt(size(x,1)),data_aAICmPFC.FcDensity_hit,'UniformOutput',1);

%% Compare FC density between mPFC, aAIC, and mPFC-aAIC
% Two-way ANOVA
data = vertcat(data_mPFC.FcDensity_hit,data_aAIC.FcDensity_hit,data_mPFCaAIC.FcDensity_hit);
num = cellfun(@numel,data,'UniformOutput',1);
Data = [];
RegionID = [];
PairID = [];
for i = 1:size(num,1)
    for j = 1:size(num,2)
        Data = [Data; data{i,j}];
        RegionID = [RegionID i*ones(1,num(i,j))];
        PairID = [PairID j*ones(1,num(i,j))];
    end
end
[p,tbl,stats] = anovan(Data,{RegionID PairID},'model','interaction','varnames',{'RegionID','PairID'});
save('FC density commparison_1msbin_between types_between regions_mPFC_aAIC_mPFCaAIC','p','tbl','stats','-v7.3'); close;
% bar plot
figure('OuterPosition',[219 303 420 534]);
bar(1.1:3.5:1.1+4*3.5,AverFcDensity_mPFC,0.25,'facecolor','k','edgecolor','k');
hold on
errorbar(1.1:3.5:1.1+4*3.5,AverFcDensity_mPFC,SemFcDensity_mPFC,'k','linestyle','none','marker','none','capsize',10);
hold on
bar(2.1:3.5:2.1+4*3.5,AverFcDensity_aAIC,0.25,'facecolor','w','edgecolor','k');
hold on
errorbar(2.1:3.5:2.1+4*3.5,AverFcDensity_aAIC,SemFcDensity_aAIC,'b','linestyle','none','marker','none','capsize',10);
hold on
bar(3.1:3.5:3.1+4*3.5,AverFcDensity_mPFCaAIC,0.25,'facecolor','r','edgecolor','r');
hold on
errorbar(3.1:3.5:3.1+4*3.5,AverFcDensity_mPFCaAIC,SemFcDensity_mPFCaAIC,'r','linestyle','none','marker','none','capsize',10);
hold on
set(gca,'XTick',zeros(1,0),'xlim',[0 18.2]);
set(gca,'YTick',0:0.05:0.2,'YTickLabel',num2cell(0:5:20),'FontName','Arial','FontSize',16,'ylim',[0 0.2]);
ylabel('FC density (%)','FontSize',18,'FontName','Arial'); box off;
set(gcf,'Renderer','Painter'); saveas(gcf,'FC density commparison_1msbin_between types_between regions_mPFC_aAIC_mPFCaAIC','fig'); close;

%% Compare FC density between mPFC-aAIC and aAIC-mPFC
% Two-way ANOVA
data = vertcat(data_mPFCaAIC.FcDensity_hit,data_aAICmPFC.FcDensity_hit);
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
save('FC density commparison_between types_between regions_mPFCaAIC_aAICmPFC','p','tbl','stats','-v7.3'); close;
% bar plot
figure('OuterPosition',[219 303 420 534]);
bar(1.1:3:1.1+4*3,AverFcDensity_mPFCaAIC,0.25,'facecolor','r','edgecolor','r');
hold on
errorbar(1.1:3:1.1+4*3,AverFcDensity_mPFCaAIC,SemFcDensity_mPFCaAIC,'r','linestyle','none','marker','none','capsize',10);
hold on
bar(2.1:3:2.1+4*3,AverFcDensity_aAICmPFC,0.25,'facecolor','k','edgecolor','k');
hold on
errorbar(2.1:3:2.1+4*3,AverFcDensity_aAICmPFC,SemFcDensity_aAICmPFC,'b','linestyle','none','marker','none','capsize',10);
set(gca,'XTick',zeros(1,0),'xlim',[0 15.2]);
set(gca,'YTick',0:0.02:0.1,'YTickLabel',num2cell(0:2:10),'FontName','Arial','FontSize',16,'ylim',[0 0.1]);
ylabel('FC density (%)','FontSize',18,'FontName','Arial'); box off;
set(gcf,'Renderer','Painter'); saveas(gcf,'FC density commparison_between types_between regions_mPFCaAIC_aAICmPFC','fig'); close;





