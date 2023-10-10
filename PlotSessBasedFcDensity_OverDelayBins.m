%% Plot FC density of five types of pairs in each delay bin.

clear; clc; close all;

%% Load result of FC density
file = dir('*SessionBased FC density_OverDelayBins_5types_hit*.mat');
for i = 1:size(file,1)
    filename = file(i,1).name; 
    if ~isempty(regexp(filename,'Within mPFC'))
        load(filename);
        FcDensity_mPFC = FcDensity_hit;
    elseif ~isempty(regexp(filename,'Within aAIC'))
        load(filename);
        FcDensity_aAIC = FcDensity_hit;
    elseif ~isempty(regexp(filename,'mPFC-aAIC'))
        load(filename);
        FcDensity_mPFCaAIC = FcDensity_hit;
    elseif ~isempty(regexp(filename,'aAIC-mPFC'))
        load(filename);
        FcDensity_aAICmPFC = FcDensity_hit;
    end
end

%% FC density of each delay bin
for iBin = 1:4
    % within mPFC
    TempFcDensity_mPFC = cell(0);
    for iType = 1:size(FcDensity_mPFC,2)
        TempFcDensity_mPFC{1,end+1} = FcDensity_mPFC{iType}{iBin};
    end
    AverFcDensity_mPFC = cellfun(@mean,TempFcDensity_mPFC,'UniformOutput',1);
    SemFcDensity_mPFC = cellfun(@(x) std(x)/sqrt(numel(x)),TempFcDensity_mPFC,'UniformOutput',1);
    % within aAIC
    TempFcDensity_aAIC = cell(0);
    for iType = 1:size(FcDensity_aAIC,2)
        TempFcDensity_aAIC{1,end+1} = FcDensity_aAIC{iType}{iBin};
    end
    AverFcDensity_aAIC = cellfun(@mean,TempFcDensity_aAIC,'UniformOutput',1);
    SemFcDensity_aAIC = cellfun(@(x) std(x)/sqrt(numel(x)),TempFcDensity_aAIC,'UniformOutput',1);
    % mPFC-aAIC
    TempFcDensity_mPFCaAIC = cell(0);
    for iType = 1:size(FcDensity_mPFCaAIC,2)
        TempFcDensity_mPFCaAIC{1,end+1} = FcDensity_mPFCaAIC{iType}{iBin};
    end
    AverFcDensity_mPFCaAIC = cellfun(@mean,TempFcDensity_mPFCaAIC,'UniformOutput',1);
    SemFcDensity_mPFCaAIC = cellfun(@(x) std(x)/sqrt(numel(x)),TempFcDensity_mPFCaAIC,'UniformOutput',1);
    % aAIC-mPFC
    TempFcDensity_aAICmPFC = cell(0);
    for iType = 1:size(FcDensity_aAICmPFC,2)
        TempFcDensity_aAICmPFC{1,end+1} = FcDensity_mPFCaAIC{iType}{iBin};
    end
    AverFcDensity_aAICmPFC = cellfun(@mean,TempFcDensity_aAICmPFC,'UniformOutput',1);
    SemFcDensity_aAICmPFC = cellfun(@(x) std(x)/sqrt(numel(x)),TempFcDensity_aAICmPFC,'UniformOutput',1);
    
    %% Compare FC density between mPFC, aAIC, and mPFC-aAIC
    % Two-way ANOVA
    data = vertcat(TempFcDensity_mPFC,TempFcDensity_aAIC,TempFcDensity_mPFCaAIC);
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
    save(sprintf('FC density comparison_at delay bin %d_between regions_mPFC_aAIC_mPFCaAIC',iBin),'p','tbl','stats','-v7.3'); close;
    % bar plot
    figure('OuterPosition',[219 303 420 534]);
    bar(1.1:3.5:1.1+4*3.5,AverFcDensity_mPFC,0.25,'facecolor','k','edgecolor','k');
    hold on
    errorbar(1.1:3.5:1.1+4*3.5,AverFcDensity_mPFC,SemFcDensity_mPFC,'k','linestyle','none','marker','none','capsize',10);
    hold on
    bar(2.1:3.5:2.1+4*3.5,AverFcDensity_aAIC,0.25,'facecolor','w','edgecolor','k');
    hold on
    errorbar(2.1:3.5:2.1+4*3.5,AverFcDensity_aAIC,SemFcDensity_aAIC,'k','linestyle','none','marker','none','capsize',10);
    hold on
    bar(3.1:3.5:3.1+4*3.5,AverFcDensity_mPFCaAIC,0.25,'facecolor','r','edgecolor','r');
    hold on
    errorbar(3.1:3.5:3.1+4*3.5,AverFcDensity_mPFCaAIC,SemFcDensity_mPFCaAIC,'r','linestyle','none','marker','none','capsize',10);
    hold on
    set(gca,'XTick',zeros(1,0),'xlim',[0 18.2]);
    set(gca,'YTick',0:0.1:1,'YTickLabel',num2cell(0:10:100),'FontName','Arial','FontSize',16,'ylim',[0 0.301]);
    ylabel('FC density (%)','FontSize',18,'FontName','Arial'); box off;
    set(gcf,'Renderer','Painter'); saveas(gcf,sprintf('FC density comparison_at delay bin %d_between regions_mPFC_aAIC_mPFCaAIC',iBin),'fig'); close;
end







