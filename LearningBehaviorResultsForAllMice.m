%% Behavioral performance for population of mice in learning phase.

clear; clc; close all;

%% Load files
[ C_Filename Pathway ] = uigetfile({'*.mat','Matlab files(*.mat)';},...
    'Pick some file','MultiSelect','on');
for iMouse = 1:size(C_Filename,2)
    C_Dataofmice{iMouse} = load(C_Filename{1,iMouse});
end
[ E_Filename Pathway ] = uigetfile({'*.mat','Matlab files(*.mat)';},...
    'Pick some file','MultiSelect','on');
for iMouse = 1:size(E_Filename,2)
    E_Dataofmice{iMouse} = load(E_Filename{1,iMouse});
end

%% Assignment
Beforetrial = 2; DRTSampleDura = 1; Delay = 4; DRTResponseDura = 0.5; RewardDura = 1; Aftertrial = 3;
Mani = 'Suppress mPFC';
C1 = [0 0 0];
CD = cd; OptoGroup = 'VGAT-ChR2';
if strcmp(OptoGroup,'NpHR')
    C2 = [0 125 0]/255;
else
    C2 = [67 106 178]/255;
end
LearningDay = 5; WellTrainedDay = 2;
TimeGain = 10;

%% Merge behavioral measurements
% performance correct rate, hit rate, CR rate, d', early lick
C_Perf = cell(1,LearningDay); C_HitRate = cell(1,LearningDay); C_CRRate = cell(1,LearningDay);  C_EarlyLickRatio = cell(1,LearningDay);
C_LickPSTH_Go = cell(1,LearningDay); C_LickPSTH_NoGo = cell(1,LearningDay);
C_LickRasterinTrial = cell(1,LearningDay); C_DRTTrialMarker = cell(1,LearningDay); C_dprime = cell(1,LearningDay);
for iDay = 1:LearningDay  % session
    for j = 1:size(C_Dataofmice,2) % mice
        if iDay <= size(C_Dataofmice{1,j}.HitRate,2) % task days
            C_Perf{iDay} = [C_Perf{iDay}; C_Dataofmice{1,j}.Perf(iDay)];
            C_HitRate{iDay} = [C_HitRate{iDay}; C_Dataofmice{1,j}.HitRate(iDay)];
            C_CRRate{iDay} = [C_CRRate{iDay}; C_Dataofmice{1,j}.CRRate(iDay)];
            C_EarlyLickRatio{iDay} = [C_EarlyLickRatio{iDay}; C_Dataofmice{1,j}.AbortedTrialsRatio(iDay)];
            tempTrialMarker = C_Dataofmice{1,j}.DRTTrialMarker{iDay};
            tempLickPSTH = C_Dataofmice{1,j}.LickPSTH{iDay};
            tempLickPSTH_Go = tempLickPSTH(tempTrialMarker==1 | tempTrialMarker==2,:);
            tempLickPSTH_NoGo = tempLickPSTH(tempTrialMarker==3 | tempTrialMarker==4,:);
            C_LickPSTH_Go{iDay} = [C_LickPSTH_Go{iDay}; TimeGain*mean(tempLickPSTH_Go)];
            C_LickPSTH_NoGo{iDay} = [C_LickPSTH_NoGo{iDay}; TimeGain*mean(tempLickPSTH_NoGo)];
            C_LickRasterinTrial{iDay} = [C_LickRasterinTrial{iDay}; {C_Dataofmice{1,j}.LickinTrial{iDay}}];
            C_DRTTrialMarker{iDay} = [C_DRTTrialMarker{iDay}; {tempTrialMarker}];
            tempHitRate = C_Dataofmice{1,j}.HitRate(iDay);
            tempFalseRate = 1-C_Dataofmice{1,j}.CRRate(iDay);
            tempTrialMarker = C_Dataofmice{1,j}.DRTTrialMarker{iDay};
            tempGoNum = nnz(tempTrialMarker==1 | tempTrialMarker==2);
            tempNogoNum = nnz(tempTrialMarker==3 | tempTrialMarker==4);
            if tempHitRate == 1
                a = 1-1/(2*tempGoNum);
            elseif tempHitRate == 0
                a = 1/(2*tempGoNum);
            else
                a = norminv(tempHitRate);
            end
            if tempFalseRate == 1
                b = 1-1/(2*tempNogoNum);
            elseif tempFalseRate == 0
                b = 1/(2*tempNogoNum);
            else
                b = norminv(tempFalseRate);
            end
            C_dprime{iDay} = [C_dprime{iDay}; a-b];
        end
    end
end
E_Perf = cell(1,LearningDay); E_HitRate = cell(1,LearningDay); E_CRRate = cell(1,LearningDay);  E_EarlyLickRatio = cell(1,LearningDay);
E_LickPSTH_Go = cell(1,LearningDay); E_LickPSTH_NoGo = cell(1,LearningDay);
E_LickRasterinTrial = cell(1,LearningDay); E_DRTTrialMarker = cell(1,LearningDay); E_dprime = cell(1,LearningDay);
for iDay = 1:LearningDay  % session
    for j = 1:size(E_Dataofmice,2) % mice
        if iDay <= size(E_Dataofmice{1,j}.HitRate,2) % task days
            E_Perf{iDay} = [E_Perf{iDay}; E_Dataofmice{1,j}.Perf(iDay)];
            E_HitRate{iDay} = [E_HitRate{iDay}; E_Dataofmice{1,j}.HitRate(iDay)];
            E_CRRate{iDay} = [E_CRRate{iDay}; E_Dataofmice{1,j}.CRRate(iDay)];
            E_EarlyLickRatio{iDay} = [E_EarlyLickRatio{iDay}; E_Dataofmice{1,j}.AbortedTrialsRatio(iDay)];
            tempTrialMarker = E_Dataofmice{1,j}.DRTTrialMarker{iDay};
            tempLickPSTH = E_Dataofmice{1,j}.LickPSTH{iDay};
            tempLickPSTH_Go = tempLickPSTH(tempTrialMarker==1 | tempTrialMarker==2,:);
            tempLickPSTH_NoGo = tempLickPSTH(tempTrialMarker==3 | tempTrialMarker==4,:);
            E_LickPSTH_Go{iDay} = [E_LickPSTH_Go{iDay}; TimeGain*mean(tempLickPSTH_Go)];
            E_LickPSTH_NoGo{iDay} = [E_LickPSTH_NoGo{iDay}; TimeGain*mean(tempLickPSTH_NoGo)];
            E_LickRasterinTrial{iDay} = [E_LickRasterinTrial{iDay}; {E_Dataofmice{1,j}.LickinTrial{iDay}}];
            E_DRTTrialMarker{iDay} = [E_DRTTrialMarker{iDay}; {E_Dataofmice{1,j}.DRTTrialMarker{iDay}}];
            tempHitRate = E_Dataofmice{1,j}.HitRate(iDay);
            tempFalseRate = 1-E_Dataofmice{1,j}.CRRate(iDay);
            tempTrialMarker = E_Dataofmice{1,j}.DRTTrialMarker{iDay};
            tempGoNum = nnz(tempTrialMarker==1 | tempTrialMarker==2);
            tempNogoNum = nnz(tempTrialMarker==3 | tempTrialMarker==4);
            if tempHitRate == 1
                a = 1-1/(2*tempGoNum);
            elseif tempHitRate == 0
                a = 1/(2*tempGoNum);
            else
                a = norminv(tempHitRate);
            end
            if tempFalseRate == 1
                b = 1-1/(2*tempNogoNum);
            elseif tempFalseRate == 0
                b = 1/(2*tempNogoNum);
            else
                b = norminv(tempFalseRate);
            end
            E_dprime{iDay} = [E_dprime{iDay}; a-b];
        end
    end
end

%% Plot session-based performance correct rate
figure('OuterPosition',[219 203 750 600]);
C_MeanPerf = cellfun(@mean, C_Perf);
C_StdPerf = cellfun(@std,C_Perf);
[nrows,ncols] = cellfun(@size,C_Perf);
errorbar(1:LearningDay,C_MeanPerf,C_StdPerf./sqrt(nrows),'color',C1,'LineWidth',2,'marker','o','markerfacecolor',C1,'markeredgecolor','none','markersize',10);
hold on
E_MeanPerf = cellfun(@mean, E_Perf);
E_StdPerf = cellfun(@std,E_Perf);
[nrows,ncols] = cellfun(@size,E_Perf);
errorbar(1:LearningDay,E_MeanPerf,E_StdPerf./sqrt(nrows),'color',C2,'LineWidth',2,'marker','o','markerfacecolor',C2,'markeredgecolor','none','markersize',10);
legend(['Vgat-ChR2 -, n = ' num2str(size(C_Dataofmice,2))],['Vgat-ChR2 +, n = ' num2str(size(E_Dataofmice,2))],'Location','southeast');
legend('boxoff');
set(gca,'XTick',0:1:LearningDay,'XTickLabel',{'0','1','2','3','4','5','6'},'FontName','Arial','FontSize',16,'xlim',[1-0.2 LearningDay+0.2]);
set(gca,'YTick',0:0.2:1,'YTickLabel',{'0','20','40','60','80','100'},'FontName','Arial','FontSize',16,'ylim',[0.4 1]);
xlabel('Training Day','FontSize',18,'FontName','Arial');
ylabel('Correct Rate ( % )','FontSize',18,'FontName','Arial');
box off;
set(gcf, 'Renderer', 'Painter'); saveas(gcf,['Correct Rate-' Mani],'fig'); close;
pValue = GetDataMatrixForMixedRepeatedAnova(C_Perf, E_Perf);
p_perf = pValue{1,1};

%% Plot session-based hit and CR rate
figure('OuterPosition',[219 303 750 600]);
C_MeanHitRate = cellfun(@mean, C_HitRate); % hit Rate
C_StdHitRate = cellfun(@std,C_HitRate);
[nrows,ncols] = cellfun(@size,C_HitRate);
errorbar(1:LearningDay,C_MeanHitRate,C_StdHitRate./sqrt(nrows),'color',C1,'linestyle','--','LineWidth',2,'marker','o','markerfacecolor','w','markeredgecolor',C1,'markersize',10);
hold on
E_MeanHitRate = cellfun(@mean, E_HitRate);
E_StdHitRate = cellfun(@std,E_HitRate);
[nrows,ncols] = cellfun(@size,E_HitRate);
errorbar(1:LearningDay,E_MeanHitRate,E_StdHitRate./sqrt(nrows),'color',C2,'linestyle','--','LineWidth',2,'marker','o','markerfacecolor','w','markeredgecolor',C2,'markersize',10);
pValue = GetDataMatrixForMixedRepeatedAnova(C_HitRate, E_HitRate);
p_HitRate = pValue{1,1};
C_MeanCRRate = cellfun(@mean, C_CRRate); % CR Rate
C_StdCRRate = cellfun(@std,C_CRRate);
[nrows,ncols] = cellfun(@size,C_CRRate);
errorbar(1:LearningDay,C_MeanCRRate,C_StdCRRate./sqrt(nrows),'color',C1,'LineWidth',2,'marker','o','markerfacecolor',C1,'markeredgecolor','none','markersize',10);
hold on
E_MeanCRRate = cellfun(@mean, E_CRRate);
E_StdCRRate = cellfun(@std,E_CRRate);
[nrows,ncols] = cellfun(@size,E_CRRate);
errorbar(1:LearningDay,E_MeanCRRate,E_StdCRRate./sqrt(nrows),'color',C2,'LineWidth',2,'marker','o','markerfacecolor',C2,'markeredgecolor','none','markersize',10);
pValue = GetDataMatrixForMixedRepeatedAnova(C_CRRate, E_CRRate);
p_CRRate = pValue{1,1};
legend(['Vgat-ChR2 -, n = ' num2str(size(C_Dataofmice,2))],['Vgat-ChR2 +, n = ' num2str(size(E_Dataofmice,2))],'Location','southeast');
legend('boxoff');
set(gca,'XTick',0:1:LearningDay,'XTickLabel',{'0','1','2','3','4','5','6'},'FontName','Arial','FontSize',16,'xlim',[1-0.2 LearningDay+0.2]);
set(gca,'YTick',0:0.2:1,'YTickLabel',{'0','20','40','60','80','100'},'FontName','Arial','FontSize',16,'ylim',[0 1]);
xlabel('Training Day','FontSize',18,'FontName','Arial');
ylabel('Rate ( % )','FontSize',18,'FontName','Arial');
box off;
set(gcf, 'Renderer', 'Painter'); saveas(gcf,['Hit and CR Rates-' Mani],'fig'); close;

%% Plot session-based d' (discriminability)
figure('OuterPosition',[219 200 400 350]);
C_Meand = cellfun(@mean, C_dprime);
C_Stdd = cellfun(@std,C_dprime);
[nrows,ncols] = cellfun(@size,C_dprime);
errorbar(1:LearningDay,C_Meand,C_Stdd./sqrt(nrows),'color',C1,'LineWidth',2,'marker','o','markerfacecolor',C1,'markeredgecolor','none','markersize',10);
hold on
E_Meand = cellfun(@mean, E_dprime);
E_Stdd = cellfun(@std,E_dprime);
[nrows,ncols] = cellfun(@size,E_dprime);
errorbar(1:LearningDay,E_Meand,E_Stdd./sqrt(nrows),'color',C2,'LineWidth',2,'marker','o','markerfacecolor',C2,'markeredgecolor','none','markersize',10);
legend(['eYFP, n = ' num2str(size(C_Dataofmice,2))],['ChR2, n = ' num2str(size(E_Dataofmice,2))],'Location','southeast');
legend('boxoff');
set(gca,'XTick',0:1:LearningDay,'XTickLabel',{'0','1','2','3','4','5','6'},'FontName','Arial','FontSize',16,'xlim',[1-0.2 LearningDay+0.2]);
set(gca,'YTick',-1:1:3,'YTickLabel',{'-1','0','1','2','3'},'FontName','Arial','FontSize',16,'ylim',[-1.2 3]);
xlabel('Learning day','FontSize',18,'FontName','Arial');
ylabel('Performance ( d )','FontSize',18,'FontName','Arial');
box off;
set(gcf, 'Renderer', 'Painter'); saveas(gcf,['dprime-' Mani],'fig');
close;
pValue = GetDataMatrixForMixedRepeatedAnova(C_dprime, E_dprime);
p_dprime = pValue{1,1};

%% Plot Day-based delay lick ratios for control and experimental group mice.
figure('OuterPosition',[219 103 750 600]);
C_MeanEarlyLickRatio = cellfun(@mean, C_EarlyLickRatio);
C_StdEarlyLickRatio = cellfun(@std,C_EarlyLickRatio);
[nrows,ncols] = cellfun(@size,C_EarlyLickRatio);
errorbar(1:LearningDay,C_MeanEarlyLickRatio,C_StdEarlyLickRatio./sqrt(nrows),'color',C1,'LineWidth',2,'marker','o','markerfacecolor',C1,'markeredgecolor','none','markersize',10);
hold on
E_MeanEarlyLickRatio = cellfun(@mean, E_EarlyLickRatio);
E_StdEarlyLickRatio = cellfun(@std,E_EarlyLickRatio);
[nrows,ncols] = cellfun(@size,E_EarlyLickRatio);
errorbar(1:LearningDay,E_MeanEarlyLickRatio,E_StdEarlyLickRatio./sqrt(nrows),'color',C2,'LineWidth',2,'marker','o','markerfacecolor',C2,'markeredgecolor','none','markersize',10);
legend(['mcherry, n = ' num2str(size(C_Dataofmice,2))],['ChR2, n = ' num2str(size(E_Dataofmice,2))],'Location','southeast');
legend('boxoff');
set(gca,'XTick',0:1:LearningDay,'XTickLabel',{'0','1','2','3','4','5','6'},'FontName','Arial','FontSize',16,'xlim',[1-0.2 LearningDay+0.2]);
set(gca,'YTick',0:0.2:1,'YTickLabel',{'0','20','40','60','80','100'},'FontName','Arial','FontSize',16,'ylim',[0 1]);
xlabel('Training Day','FontSize',18,'FontName','Arial');
ylabel('Proportion of abortion ( % )','FontSize',18,'FontName','Arial');
box off;
pValue = GetDataMatrixForMixedRepeatedAnova(C_EarlyLickRatio, E_EarlyLickRatio);
p_AbortionRatio = pValue{1,1};
title(sprintf('p=%d',p_AbortionRatio));
set(gcf, 'Renderer', 'Painter');
saveas(gcf,sprintf('Ratio of abortion-%s',Mani),'fig'); close all;
save('TwoGroupsMiceIDandPValue','C_Filename','E_Filename','p_perf','p_HitRate','p_CRRate','p_dprime','p_AbortionRatio');
