%% Lick efficiency, following optogenetic manipulation.

clear; clc; close all;

%% Load files
% control group
[C_Filename Pathway] = uigetfile({'*.mat','Matlab files(*.mat)';},'Pick some file','MultiSelect','on');
for iSess = 1:size(C_Filename,2)
    C_Dataofmice{iSess} = load(C_Filename{1,iSess});
end
% experimental group
[E_Filename Pathway] = uigetfile({'*.mat','Matlab files(*.mat)';},'Pick some file','MultiSelect','on');
for iSess = 1:size(E_Filename,2)
    E_Dataofmice{iSess} = load(E_Filename{1,iSess});
end

%% Assignment
Mani = 'suppress mPFC-to-aAIC';
C1 = [0 0 0];
OptoGroup = 'NpHR';
if ~isempty(regexp(OptoGroup,'NpHR'))
    C2 = [0 125 0]/255;
else
    C2 = [67 106 178]/255;
end
Learningday = 6; WellTrainedDay = 2; BeforeTrial = 2; SampleDuration = 1; DelayDuration = 10; RespOdorDuration = 0.5;
LearningOrWellTrain = 1; % learning:1; well-trained:2
RespWindow = [BeforeTrial + SampleDuration + DelayDuration + RespOdorDuration BeforeTrial + SampleDuration + DelayDuration + RespOdorDuration + 0.4];

%% Lick number in the response window
if LearningOrWellTrain == 1 % learning phase
    C_LickNuminWindow = cell(1,Learningday);
    for iSess = 1:Learningday  % session
        for j = 1:size(C_Dataofmice,2) % Mice
            if iSess <= size(C_Dataofmice{1,j}.HitRate,2) % //choosing only task days
                tempTrialMarker = C_Dataofmice{1,j}.DRTTrialMarker{iSess};
                tempLickinGoTrials = find(tempTrialMarker==1 | tempTrialMarker==2); % Go trials
                tempLickinGoTrials = C_Dataofmice{1,j}.LickinTrial{iSess}(:,tempLickinGoTrials);
                tempLickNuminGoTrials = cellfun(@(x) nnz(x>RespWindow(1) & x<=RespWindow(2)),tempLickinGoTrials,'UniformOutput',1);
                tempLickinNogoTrials = find(tempTrialMarker==3 | tempTrialMarker==4); % No-go trials
                tempLickinNogoTrials = C_Dataofmice{1,j}.LickinTrial{iSess}(:,tempLickinNogoTrials);
                tempLickNuminNogoTrials = cellfun(@(x) nnz(x>RespWindow(1) & x<=RespWindow(2)),tempLickinNogoTrials,'UniformOutput',1);
                C_LickNuminWindow{iSess} = [C_LickNuminWindow{iSess}; horzcat(sum(tempLickNuminGoTrials,2),sum(tempLickNuminNogoTrials,2))];
            end
        end
    end
    E_LickNuminWindow = cell(1,Learningday);
    for iSess = 1:Learningday  % Day
        for j = 1:size(E_Dataofmice,2) % Mice
            if iSess <= size(E_Dataofmice{1,j}.HitRate,2) % //choosing only task days
                tempTrialMarker = E_Dataofmice{1,j}.DRTTrialMarker{iSess};
                tempLickinGoTrials = find(tempTrialMarker==1 | tempTrialMarker==2); % Go trials
                tempLickinGoTrials = E_Dataofmice{1,j}.LickinTrial{iSess}(:,tempLickinGoTrials);
                tempLickNuminGoTrials = cellfun(@(x) nnz(x>RespWindow(1) & x<=RespWindow(2)),tempLickinGoTrials,'UniformOutput',1);
                tempLickinNogoTrials = find(tempTrialMarker==3 | tempTrialMarker==4); % No-go trials
                tempLickinNogoTrials = E_Dataofmice{1,j}.LickinTrial{iSess}(:,tempLickinNogoTrials);
                tempLickNuminNogoTrials = cellfun(@(x) nnz(x>RespWindow(1) & x<=RespWindow(2)),tempLickinNogoTrials,'UniformOutput',1);
                E_LickNuminWindow{iSess} = [E_LickNuminWindow{iSess}; horzcat(sum(tempLickNuminGoTrials,2),sum(tempLickNuminNogoTrials,2))];
            end
        end
    end
    LickingEfficiency_Ctrl = cellfun(@(x) x(:,1)./sum(x,2),C_LickNuminWindow,'UniformOutput',0);
    LickingEfficiency_Opto = cellfun(@(x) x(:,1)./sum(x,2),E_LickNuminWindow,'UniformOutput',0);
    
    %% Plot session-based lick efficiency
    figure('OuterPosition',[219 303 350 400]);
    C_MeanPerf = cellfun(@mean, LickingEfficiency_Ctrl);
    C_StdPerf = cellfun(@std,LickingEfficiency_Ctrl);
    [nrows,ncols] = cellfun(@size,LickingEfficiency_Ctrl);
    errorbar(1:Learningday,C_MeanPerf,C_StdPerf./sqrt(nrows),'color',C1,'LineWidth',2,'marker','o','markerfacecolor',C1,'markeredgecolor','none','markersize',10);
    hold on
    E_MeanPerf = cellfun(@mean, LickingEfficiency_Opto);
    E_StdPerf = cellfun(@std,LickingEfficiency_Opto);
    [nrows,ncols] = cellfun(@size,LickingEfficiency_Opto);
    errorbar(1:Learningday,E_MeanPerf,E_StdPerf./sqrt(nrows),'color',C2,'LineWidth',2,'marker','o','markerfacecolor',C2,'markeredgecolor','none','markersize',10);
    legend(['Vgat-ChR2 -, n = ' num2str(size(C_Dataofmice,2))],['Vgat-ChR2 +, n = ' num2str(size(E_Dataofmice,2))],'Location','southeast');
    legend('boxoff');
    set(gca,'XTick',0:1:5,'XTickLabel',{'0','1','2','3','4','5'},'FontName','Arial','FontSize',16,'xlim',[1-0.2 Learningday+0.2]);
    set(gca,'YTick',0.5:0.1:1,'YTickLabel',{'50','60','70','80','90','100'},'FontName','Arial','FontSize',16,'ylim',[0.47 1]);
    xlabel('Learning Day','FontSize',18,'FontName','Arial');
    ylabel('Lick efficiency ( % )','FontSize',18,'FontName','Arial');
    box off;
    set(gcf, 'Renderer', 'Painter'); saveas(gcf,['Lick efficiency-' Mani],'fig'); close;
    pValue = GetDataMatrixForMixedRepeatedAnova(LickingEfficiency_Ctrl,LickingEfficiency_Opto,strcat('Lick efficiency-',Mani)); % Tw-ANOVA-md
    p_lickefficiency = pValue{1,1};
else % well-trained phase
    ControlMiceNum = length(C_Dataofmice); OptoMiceNum = length(E_Dataofmice); interval = 0.1; popinterval = 0.3;
    XPosttionforControlLaserOff = (1-interval*floor(ControlMiceNum/2)):interval:(1-interval*floor(ControlMiceNum/2)+interval*(ControlMiceNum-1));
    XPosttionforControlLaserOn = (2-interval*floor(ControlMiceNum/2)):interval:(2-interval*floor(ControlMiceNum/2)+interval*(ControlMiceNum-1));
    XPosttionforOptoLaserOff = (4-interval*floor(OptoMiceNum/2)):interval:(4-interval*floor(OptoMiceNum/2)+interval*(OptoMiceNum-1));
    XPosttionforOptoLaserOn = (5-interval*floor(OptoMiceNum/2)):interval:(5-interval*floor(OptoMiceNum/2)+interval*(OptoMiceNum-1));
    TrialNumPerBlock = 20;
    % control group
    C_LickNuminWindow = [];
    for j = 1:size(C_Dataofmice,2) % mice
        if iscell(C_Dataofmice{1,j}.DRTTrialMarker)
            tempTrialMarker = C_Dataofmice{1,j}.DRTTrialMarker{1,WellTrainedDay};
        else
            tempTrialMarker = C_Dataofmice{1,j}.DRTTrialMarker;
        end
        if length(C_Dataofmice{1,j}.LickinTrial) <= 5
            tempLick = C_Dataofmice{1,j}.LickinTrial{1,WellTrainedDay};
        else
            tempLick = C_Dataofmice{1,j}.LickinTrial;
        end
        tempBlockNum = length(tempTrialMarker)/TrialNumPerBlock;
        tempTrialMarker_laseroff = []; tempTrialMarker_laseron = [];
        tempLicklaseroff = []; tempLicklaseron = [];
        for k = 1:2:tempBlockNum
            tempTrialMarker_laseroff = [tempTrialMarker_laseroff; tempTrialMarker(TrialNumPerBlock*(k-1)+1:TrialNumPerBlock*k)];
            tempLicklaseroff = [tempLicklaseroff tempLick(:,TrialNumPerBlock*(k-1)+1:TrialNumPerBlock*k)];
        end
        for k = 2:2:tempBlockNum
            tempTrialMarker_laseron = [tempTrialMarker_laseron; tempTrialMarker(TrialNumPerBlock*(k-1)+1:TrialNumPerBlock*k)];
            tempLicklaseron = [tempLicklaseron tempLick(:,TrialNumPerBlock*(k-1)+1:TrialNumPerBlock*k)];
        end
        % laser off condition
        tempGoTrialMarker = find(tempTrialMarker_laseroff==1 | tempTrialMarker_laseroff==2);
        tempNogoTrialMarker = find(tempTrialMarker_laseroff==3 | tempTrialMarker_laseroff==4);
        LickinGoTrials_laseroff = tempLicklaseroff(:,tempGoTrialMarker);
        LickinNogoTrials_laseroff = tempLicklaseroff(tempNogoTrialMarker);
        % laser on condition
        tempGoTrialMarker = find(tempTrialMarker_laseron==1 | tempTrialMarker_laseron==2);
        tempNogoTrialMarker = find(tempTrialMarker_laseron==3 | tempTrialMarker_laseron==4);
        LickinGoTrials_laseron = tempLicklaseron(:,tempGoTrialMarker);
        LickinNogoTrials_laseron = tempLicklaseron(tempNogoTrialMarker);
        LickinGoTrials_laseroff = cellfun(@(x) nnz(x>RespWindow(1) & x<=RespWindow(2)),LickinGoTrials_laseroff,'UniformOutput',1);
        LickinNogoTrials_laseroff = cellfun(@(x) nnz(x>RespWindow(1) & x<=RespWindow(2)),LickinNogoTrials_laseroff,'UniformOutput',1);
        LickinGoTrials_laseron = cellfun(@(x) nnz(x>RespWindow(1) & x<=RespWindow(2)),LickinGoTrials_laseron,'UniformOutput',1);
        LickinNogoTrials_laseron = cellfun(@(x) nnz(x>RespWindow(1) & x<=RespWindow(2)),LickinNogoTrials_laseron,'UniformOutput',1);
        C_LickNuminWindow = [C_LickNuminWindow; horzcat(sum(LickinGoTrials_laseroff,2)/(sum(LickinGoTrials_laseroff,2)+sum(LickinNogoTrials_laseroff,2)),sum(LickinGoTrials_laseron,2)/(sum(LickinGoTrials_laseron,2)+sum(LickinNogoTrials_laseron,2)))];
    end
    C_LickNuminWindow(any(isnan(C_LickNuminWindow),2) == 1,:) = [];
    % experimental group
    E_LickNuminWindow = [];
    for j = 1:size(E_Dataofmice,2) % mice
        if iscell(E_Dataofmice{1,j}.DRTTrialMarker)
            tempTrialMarker = E_Dataofmice{1,j}.DRTTrialMarker{1,WellTrainedDay};
        else
            tempTrialMarker = E_Dataofmice{1,j}.DRTTrialMarker;
        end
        if length(E_Dataofmice{1,j}.LickinTrial)<=5
            tempLick = E_Dataofmice{1,j}.LickinTrial{1,WellTrainedDay};
        else
            tempLick = E_Dataofmice{1,j}.LickinTrial;
        end
        tempBlockNum = length(tempTrialMarker)/TrialNumPerBlock;
        tempTrialMarker_laseroff = []; tempTrialMarker_laseron = [];
        tempLicklaseroff = []; tempLicklaseron = [];
        for k = 1:2:tempBlockNum
            tempTrialMarker_laseroff = [tempTrialMarker_laseroff; tempTrialMarker(TrialNumPerBlock*(k-1)+1:TrialNumPerBlock*k)];
            tempLicklaseroff = [tempLicklaseroff tempLick(:,TrialNumPerBlock*(k-1)+1:TrialNumPerBlock*k)];
        end
        for k = 2:2:tempBlockNum
            tempTrialMarker_laseron = [tempTrialMarker_laseron; tempTrialMarker(TrialNumPerBlock*(k-1)+1:TrialNumPerBlock*k)];
            tempLicklaseron = [tempLicklaseron tempLick(:,TrialNumPerBlock*(k-1)+1:TrialNumPerBlock*k)];
        end
        % laser off condition
        tempGoTrialMarker = find(tempTrialMarker_laseroff==1 | tempTrialMarker_laseroff==2);
        tempNogoTrialMarker = find(tempTrialMarker_laseroff==3 | tempTrialMarker_laseroff==4);
        LickinGoTrials_laseroff = tempLicklaseroff(:,tempGoTrialMarker);
        LickinNogoTrials_laseroff = tempLicklaseroff(tempNogoTrialMarker);
        % laser on condition
        tempGoTrialMarker = find(tempTrialMarker_laseron==1 | tempTrialMarker_laseron==2);
        tempNogoTrialMarker = find(tempTrialMarker_laseron==3 | tempTrialMarker_laseron==4);
        LickinGoTrials_laseron = tempLicklaseron(:,tempGoTrialMarker);
        LickinNogoTrials_laseron = tempLicklaseron(tempNogoTrialMarker);
        LickinGoTrials_laseroff = cellfun(@(x) nnz(x>RespWindow(1) & x<=RespWindow(2)),LickinGoTrials_laseroff,'UniformOutput',1);
        LickinNogoTrials_laseroff = cellfun(@(x) nnz(x>RespWindow(1) & x<=RespWindow(2)),LickinNogoTrials_laseroff,'UniformOutput',1);
        LickinGoTrials_laseron = cellfun(@(x) nnz(x>RespWindow(1) & x<=RespWindow(2)),LickinGoTrials_laseron,'UniformOutput',1);
        LickinNogoTrials_laseron = cellfun(@(x) nnz(x>RespWindow(1) & x<=RespWindow(2)),LickinNogoTrials_laseron,'UniformOutput',1);
        E_LickNuminWindow = [E_LickNuminWindow; horzcat(sum(LickinGoTrials_laseroff,2)/(sum(LickinGoTrials_laseroff,2)+sum(LickinNogoTrials_laseroff,2)),sum(LickinGoTrials_laseron,2)/(sum(LickinGoTrials_laseron,2)+sum(LickinNogoTrials_laseron,2)))];
    end
    E_LickNuminWindow(any(isnan(E_LickNuminWindow),2) == 1,:) = [];
    PlotWellTrainedBehaviourData(C_LickNuminWindow,E_LickNuminWindow,0.1,0.3,OptoGroup,0,0,0);
    ylabel('Lick efficiency','FontSize',18,'FontName','Arial');
    box off;
    
    %% Wilcoxon signed-rank test and Mann-Whitney U test
    p1=signrank(C_LickNuminWindow(:,1),C_LickNuminWindow(:,2));
    p2=signrank(E_LickNuminWindow(:,1),E_LickNuminWindow(:,2));
    diffinCtrl = C_LickNuminWindow(:,2)-C_LickNuminWindow(:,1);
    diffinInhibition = E_LickNuminWindow(:,2)-E_LickNuminWindow(:,1);
    p3 = ranksum(diffinCtrl,diffinInhibition);
    title(sprintf('p1=%d p2=%d p3=%d',p1,p2,p3));
    set(gca,'YTick',0.5:0.1:1,'YTickLabel',{'50','60','70','80','90','100'},'FontName','Arial','FontSize',16,'ylim',[0.47 1]);
    set(gcf, 'Renderer', 'Painter'); saveas(gcf,['Welltrained lick efficiency-' Mani],'fig'); close;    
end



