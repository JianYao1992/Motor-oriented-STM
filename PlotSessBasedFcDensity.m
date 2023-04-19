%% Plot FC density of five types of pairs.

clear; clc; close all;

%% Load result of FC density
file = dir('*SessionBased FC density_5types_hit*.mat');
for i = 1:size(file,1)
    filename = file(i,1).name; 
    if ~isempty(regexp(filename,'Within mPFC'))
        data_mPFC = load(filename);
    elseif ~isempty(regexp(filename,'Within aAIC'))
        data_aAIC = load(filename);
    elseif ~isempty(regexp(filename,'mPFC-aAIC'))
        data_Proj = load(filename);
    end
end

%% Error bar plot
AverFcDensity_mPFC = cellfun(@mean,data_mPFC.FcDensity_hit,'UniformOutput',1); % mPFC
SemFcDensity_mPFC = cellfun(@(x) std(x)/sqrt(size(x,1)),data_mPFC.FcDensity_hit,'UniformOutput',1);
AverFcDensity_aAIC = cellfun(@mean,data_aAIC.FcDensity_hit,'UniformOutput',1); % aAIC
SemFcDensity_aAIC = cellfun(@(x) std(x)/sqrt(size(x,1)),data_aAIC.FcDensity_hit,'UniformOutput',1);
AverFcDensity_Proj = cellfun(@mean,data_Proj.FcDensity_hit,'UniformOutput',1); % mPFC-aAIC
SemFcDensity_Proj = cellfun(@(x) std(x)/sqrt(size(x,1)),data_Proj.FcDensity_hit,'UniformOutput',1);
figure('OuterPosition',[219 303 420 534]);
bar(1.1:3.5:1.1+4*3.5,AverFcDensity_mPFC,0.25,'facecolor','k','edgecolor','k');
hold on
errorbar(1.1:3.5:1.1+4*3.5,AverFcDensity_mPFC,SemFcDensity_mPFC,'k','linestyle','none','marker','none','capsize',10);
hold on
bar(2.1:3.5:2.1+4*3.5,AverFcDensity_aAIC,0.25,'facecolor','w','edgecolor','k');
hold on
errorbar(2.1:3.5:2.1+4*3.5,AverFcDensity_aAIC,SemFcDensity_aAIC,'b','linestyle','none','marker','none','capsize',10);
hold on
bar(3.1:3.5:3.1+4*3.5,AverFcDensity_Proj,0.25,'facecolor','r','edgecolor','r');
hold on
errorbar(3.1:3.5:3.1+4*3.5,AverFcDensity_Proj,SemFcDensity_Proj,'r','linestyle','none','marker','none','capsize',10);
hold on
set(gca,'XTick',zeros(1,0),'xlim',[0 18.2]);
set(gca,'YTick',0:0.05:0.2,'YTickLabel',num2cell(0:5:20),'FontName','Arial','FontSize',16,'ylim',[0 0.2]);
ylabel('FC density (%)','FontSize',18,'FontName','Arial'); box off;
set(gcf,'Renderer','Painter'); saveas(gcf,'FC density commparison_between types_between regions','fig'); close;

%% Two-way ANOVA
data = vertcat(data_mPFC.CD_Hit_alltypes,data_aAIC.CD_Hit_alltypes,data_Proj.CD_Hit_alltypes);
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
save('FC density commparison_between types_between regions','p','tbl','stats','-v7.3'); close;





