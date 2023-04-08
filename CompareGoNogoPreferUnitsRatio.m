%% Compare proportion of Go-preferred and NoGo-preferred neurons, following optogenetic manipulation

clear; clc; close all;

%% Enter path
currpath = uigetdir;
cd(currpath);

%% Assignment
ToCompGroup = 'ChR2Group';
TarReg = 'aAIC';
if strcmp(ToCompGroup,'NpHRGroup')
    C = [0 125 0]/255;
elseif strcmp(ToCompGroup,'ChR2Group')
    C = [67 106 178]/255;
end

%% Load result of Go- and NoGo-preferred neurons
% control group
load('Go-NoGo-PreferredUnitsID-CtrlGroup.mat');
if strcmp(TarReg,'mPFC')
    Data_Ctrl = PreferCodingNeuronID.mPFC;
elseif strcmp(TarReg,'aAIC')
    Data_Ctrl = PreferCodingNeuronID.aAIC;
end
% experimental group
load(sprintf('Go-NoGo-PreferredUnitsID-%s.mat',ToCompGroup));
if strcmp(TarReg,'mPFC')
    Data_Exp = PreferCodingNeuronID.mPFC;
elseif strcmp(TarReg,'aAIC')
    Data_Exp = PreferCodingNeuronID.aAIC;
end

%% Compare Go-preferred and NoGo-perferred neurons
% Go-preferred neurons
GoPrefUnitsNum_Ctrl = numel(Data_Ctrl.GoPrefUnitsID);
TotalUnitsNum_Ctrl = Data_Ctrl.totalNum;
GoPrefUnitsNum_Exp = numel(Data_Exp.GoPrefUnitsID);
TotalUnitsNum_Exp = Data_Exp.totalNum;
[h,p_go,chi2stat,df] = prop_test([GoPrefUnitsNum_Ctrl GoPrefUnitsNum_Exp],[TotalUnitsNum_Ctrl TotalUnitsNum_Exp], false);
% NoGo-preferred neurons
NoGoPrefUnitsNum_Ctrl = numel(Data_Ctrl.NoGoPrefUnitsID);
NoGoPrefUnitsNum_Exp = numel(Data_Exp.NoGoPrefUnitsID);
[h,p_nogo,chi2stat,df] = prop_test([NoGoPrefUnitsNum_Ctrl NoGoPrefUnitsNum_Exp],[TotalUnitsNum_Ctrl TotalUnitsNum_Exp], false);
% bar plot
maxYvalue = 1.1*max(horzcat(GoPrefUnitsNum_Ctrl/TotalUnitsNum_Ctrl,GoPrefUnitsNum_Exp/TotalUnitsNum_Exp,NoGoPrefUnitsNum_Ctrl/TotalUnitsNum_Ctrl,NoGoPrefUnitsNum_Exp/TotalUnitsNum_Exp));
figure('OuterPosition',[219 303 420 534]);
bar(1.1,GoPrefUnitsNum_Ctrl/TotalUnitsNum_Ctrl,0.6,'k','edgecolor','none');
hold on
bar(1.9,GoPrefUnitsNum_Exp/TotalUnitsNum_Exp,0.6,'facecolor',C,'edgecolor','none');
hold on
bar(2.7,NoGoPrefUnitsNum_Ctrl/TotalUnitsNum_Ctrl,0.6,'k','edgecolor','none');
hold on
bar(3.5,NoGoPrefUnitsNum_Exp/TotalUnitsNum_Exp,0.6,'facecolor',C,'edgecolor','none');
hold on
title(['p_go = ' num2str(p_go) '; p_nogo = ' num2str(p_nogo)]);
set(gca,'XTick',zeros(1,0),'xlim',[0.6 4]);
set(gca,'YTick',0:0.1:1,'YTickLabel',num2cell(0:0.1:1),'FontName','Arial','FontSize',16,'ylim',[0 maxYvalue]);
ylabel('Proportion ( % )','FontSize',18,'FontName','Arial');
box off;
set(gcf,'Renderer','Painter'); saveas(gcf,['Compare-GoandNoGo-Preferred' TarReg 'NeuronsBetweenCtrland' ToCompGroup],'fig'); close;




