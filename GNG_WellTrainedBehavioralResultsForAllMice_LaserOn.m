%% Plot OGNG behavioral result in well-trained phase ////// in the condition of laser in all trials //////

clear; clc; close all;

%% Enter path
CurrPath = uigetdir;
cd(CurrPath);


%% defination
Task = 'OGNG';
Mani = 'activate mPFC-aAIC';
C1 = [0 0 0];
if strcmp(Mani,'suppress mPFC-aAIC')
    C2 = [0 125 0]/255;
else
    C2 = [67 106 178]/255;
end

%% Hit rate and CR rate in well-trained phase
TrialNumInUse = 40;
HitRate_Ctrl = [];
CRRate_Ctrl = [];
HitRate_Exp = [];
CRRate_Exp = [];
% load result
[C_Filename Pathway] = uigetfile({'*.mat','Matlab files(*.mat)';},'Pick some file','MultiSelect','on');
for iMouse = 1:size(C_Filename,2)
    load(C_Filename{1,iMouse});
    tempOutcome = GNGTrialMarker{1}(end-(TrialNumInUse-1):end);
    tempHitNum = nnz(tempOutcome==1);
    tempMissNum = nnz(tempOutcome==2);
    tempFalseNum = nnz(tempOutcome==3);
    tempCRNum = nnz(tempOutcome==4)
    HitRate_Ctrl = [HitRate_Ctrl; tempHitNum/(tempHitNum+tempMissNum)];
    CRRate_Ctrl = [CRRate_Ctrl; tempCRNum/(tempFalseNum+tempCRNum)];
end
[E_Filename Pathway] = uigetfile({'*.mat','Matlab files(*.mat)';},'Pick some file','MultiSelect','on');
for iMouse = 1:size(E_Filename,2)
    load(E_Filename{1,iMouse});
    tempOutcome = GNGTrialMarker{1}(end-(TrialNumInUse-1):end);
    tempHitNum = nnz(tempOutcome==1);
    tempMissNum = nnz(tempOutcome==2);
    tempFalseNum = nnz(tempOutcome==3);
    tempCRNum = nnz(tempOutcome==4)
    HitRate_Exp = [HitRate_Exp; tempHitNum/(tempHitNum+tempMissNum)];
    CRRate_Exp = [CRRate_Exp; tempCRNum/(tempFalseNum+tempCRNum)];
end

%% Mann-Whitney U test
p_hit = ranksum(HitRate_Ctrl,HitRate_Exp);
p_cr = ranksum(CRRate_Ctrl,CRRate_Exp);

%% Plot behavioral performance
figure('position',[200 200 500 300]);
errorbar(1,mean(HitRate_Ctrl),std(HitRate_Ctrl)/size(HitRate_Ctrl,1),'color',C1,'marker','o','markerfacecolor',[1 1 1],'markeredgecolor',C1);  
hold on
errorbar(1,mean(HitRate_Exp),std(HitRate_Exp)/size(HitRate_Exp,1),'color',C2,'marker','o','markerfacecolor',[1 1 1],'markeredgecolor',C2); 
hold on
errorbar(2,mean(CRRate_Ctrl),std(CRRate_Ctrl)/size(CRRate_Ctrl,1),'color',C1,'marker','o','markerfacecolor',C1,'markeredgecolor',C1);  
hold on
errorbar(2,mean(CRRate_Exp),std(CRRate_Exp)/size(CRRate_Exp,1),'color',C2,'marker','o','markerfacecolor',C2,'markeredgecolor',C2); 
hold on
set(gca,'XTick',[1 2],'XLim',[0 3],'YLim',[0 1]);
box off
title(['ranksum p_hit = ' num2str(p_hit) '-p_cr = ' num2str(p_cr)]);
set(gcf,'Render','Painter'); saveas(gcf,['Hit and CR rates-' Mani '-WellTrained-' Task],'fig'); close all;
