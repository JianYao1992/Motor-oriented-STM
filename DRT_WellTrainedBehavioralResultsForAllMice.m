%% Behavioral results of population of mice in the well-trained phase

clear; clc; close all;

%% Load files
% control group
[ C_Filename Pathway ] = uigetfile({'*.mat','Matlab files(*.mat)';},...
    'Pick some file','MultiSelect','on');
for iDay = 1:size(C_Filename,2)
    C_Dataofmice{iDay} = load(C_Filename{1,iDay});
end
% experimental group
[ E_Filename Pathway ] = uigetfile({'*.mat','Matlab files(*.mat)';},...
    'Pick some file','MultiSelect','on');
for iDay = 1:size(E_Filename,2)
    E_Dataofmice{iDay} = load(E_Filename{1,iDay});
end

%% Assignment
Mani = 'Suppress mPFC or aAIC';
OptoGroup = 'VGAT-ChR2';
C1 = [0 0 0];
if strcmp(OptoGroup,'NpHR')
    C2 = [0 125 0]/255;
else
    C2 = [67 106 178]/255;
end
WelltrainedDayID = 1; TrialNumPerBlock = 20;

%% Performance correct rate
c_perf = []; e_perf = [];
for j = 1:size(C_Dataofmice,2) % control group
    if strcmp(class(C_Dataofmice{1,1}.Perf),'cell')
        TotalBlockNum = size(C_Dataofmice{1,j}.Perf{1,WelltrainedDayID},2);
        c_perf(j,1) =  mean(C_Dataofmice{1,j}.Perf{1,WelltrainedDayID}(:,1:2:TotalBlockNum-1)); % laser off
        c_perf(j,2) =  mean(C_Dataofmice{1,j}.Perf{1,WelltrainedDayID}(:,2:2:TotalBlockNum)); % laser on
    else
        TotalBlockNum = size(C_Dataofmice{1,j}.Perf,2);
        c_perf(j,1) =  mean(C_Dataofmice{1,j}.Perf(:,1:2:TotalBlockNum-1)); % laser off
        c_perf(j,2) =  mean(C_Dataofmice{1,j}.Perf(:,2:2:TotalBlockNum)); % laser on
    end
end
for j = 1:size(E_Dataofmice,2) % experimental group
    if strcmp(class(C_Dataofmice{1,1}.Perf),'cell')
        TotalBlockNum = size(E_Dataofmice{1,j}.Perf{1,WelltrainedDayID},2);
        e_perf(j,1) =  mean(E_Dataofmice{1,j}.Perf{1,WelltrainedDayID}(:,1:2:TotalBlockNum-1)); % laser off
        e_perf(j,2) =  mean(E_Dataofmice{1,j}.Perf{1,WelltrainedDayID}(:,2:2:TotalBlockNum)); % laser on
    else
        TotalBlockNum = size(E_Dataofmice{1,j}.Perf,2);
        e_perf(j,1) =  mean(E_Dataofmice{1,j}.Perf(:,1:2:TotalBlockNum-1)); % laser off
        e_perf(j,2) =  mean(E_Dataofmice{1,j}.Perf(:,2:2:TotalBlockNum)); % laser on
    end
end
PlotWellTrainedResults(c_perf,e_perf,0.1,0.3,OptoGroup,0,0,0);
p1_perf=signrank(c_perf(:,1),c_perf(:,2));
p2_perf=signrank(e_perf(:,1),e_perf(:,2));
diffinCtrl = c_perf(:,2)-c_perf(:,1);
diffinExp = e_perf(:,2)-e_perf(:,1);
p3_perf = ranksum(diffinCtrl,diffinExp);
title(sprintf('p1=%d p2=%d p3=%d',p1_perf,p2_perf,p3_perf));
ylabel('Correct Rate ( % )','FontSize',18,'FontName','Arial');
ylim([0.49 1]);
box off;
set(gcf, 'Renderer', 'Painter'); saveas(gcf,['Welltrained Correct Rate-' Mani],'fig'); close;

%% Hit rate
c_hit = []; e_hit = [];
for j = 1:size(C_Dataofmice,2) % control
    if strcmp(class(C_Dataofmice{1,1}.HitRate),'cell')
        TotalBlockNum = size(C_Dataofmice{1,j}.HitRate{1,WelltrainedDayID},2);
        c_hit(j,1) =  mean(C_Dataofmice{1,j}.HitRate{1,WelltrainedDayID}(:,1:2:TotalBlockNum-1)); % laser off
        c_hit(j,2) =  mean(C_Dataofmice{1,j}.HitRate{1,WelltrainedDayID}(:,2:2:TotalBlockNum)); % laser on
    else
        TotalBlockNum = size(C_Dataofmice{1,j}.HitRate,2);
        c_hit(j,1) =  mean(C_Dataofmice{1,j}.HitRate(:,1:2:TotalBlockNum-1)); % laser off
        c_hit(j,2) =  mean(C_Dataofmice{1,j}.HitRate(:,2:2:TotalBlockNum)); % laser on
    end
end
for j = 1:size(E_Dataofmice,2) % experimental
    if strcmp(class(E_Dataofmice{1,1}.HitRate),'cell')
        TotalBlockNum = size(E_Dataofmice{1,j}.HitRate{1,WelltrainedDayID},2);
        e_hit(j,1) =  mean(E_Dataofmice{1,j}.HitRate{1,WelltrainedDayID}(:,1:2:TotalBlockNum-1)); % laser off
        e_hit(j,2) =  mean(E_Dataofmice{1,j}.HitRate{1,WelltrainedDayID}(:,2:2:TotalBlockNum)); % laser on
    else
        TotalBlockNum = size(E_Dataofmice{1,j}.HitRate,2);
        e_hit(j,1) =  mean(E_Dataofmice{1,j}.HitRate(:,1:2:TotalBlockNum-1)); % laser off
        e_hit(j,2) =  mean(E_Dataofmice{1,j}.HitRate(:,2:2:TotalBlockNum)); % laser on
    end
end
p1_Hit=signrank(c_hit(:,1),c_hit(:,2));
p2_Hit=signrank(e_hit(:,1),e_hit(:,2));
diffinCtrl = c_hit(:,2)-c_hit(:,1);
diffinExp = e_hit(:,2)-e_hit(:,1);
p3_Hit = ranksum(diffinCtrl,diffinExp);

%% CR rate
c_rej = []; e_false = [];
for j = 1:size(C_Dataofmice,2) % control
    if strcmp(class(C_Dataofmice{1,1}.CRRate),'cell')
        TotalBlockNum = size(C_Dataofmice{1,j}.CRRate{1,WelltrainedDayID},2);
        c_rej(j,1) =  mean(C_Dataofmice{1,j}.CRRate{1,WelltrainedDayID}(:,1:2:TotalBlockNum-1)); % laser off
        c_rej(j,2) =  mean(C_Dataofmice{1,j}.CRRate{1,WelltrainedDayID}(:,2:2:TotalBlockNum)); % laser on
    else
        TotalBlockNum = size(C_Dataofmice{1,j}.CRRate,2);
        c_rej(j,1) =  mean(C_Dataofmice{1,j}.CRRate(:,1:2:TotalBlockNum-1)); % laser off
        c_rej(j,2) =  mean(C_Dataofmice{1,j}.CRRate(:,2:2:TotalBlockNum)); % laser on
    end
end
for j = 1:size(E_Dataofmice,2) % experimental
    if strcmp(class(E_Dataofmice{1,1}.CRRate),'cell')
        TotalBlockNum = size(E_Dataofmice{1,j}.CRRate{1,WelltrainedDayID},2);
        e_rej(j,1) =  mean(E_Dataofmice{1,j}.CRRate{1,WelltrainedDayID}(:,1:2:TotalBlockNum-1)); % OFF
        e_rej(j,2) =  mean(E_Dataofmice{1,j}.CRRate{1,WelltrainedDayID}(:,2:2:TotalBlockNum)); % ON
    else
        TotalBlockNum = size(E_Dataofmice{1,j}.CRRate,2);
        e_rej(j,1) =  mean(E_Dataofmice{1,j}.CRRate(:,1:2:TotalBlockNum-1)); % OFF
        e_rej(j,2) =  mean(E_Dataofmice{1,j}.CRRate(:,2:2:TotalBlockNum)); % ON
    end
end
p1_CR=signrank(c_rej(:,1),c_rej(:,2));
p2_CR=signrank(e_rej(:,1),e_rej(:,2));
diffinCtrl = c_rej(:,2)-c_rej(:,1);
diffinExp = e_rej(:,2)-e_rej(:,1);
p3_CR = ranksum(diffinCtrl,diffinExp);
figure
PlotWellTrainedResults(c_hit,e_hit,0.1,0.3,OptoGroup,0,0,1); hold on
hold on
PlotWellTrainedResults(c_rej,e_rej,0.1,0.3,OptoGroup,0,0,0);
ylim([0 1]);
title(sprintf('p1_Hit=%d p1_CR=%d p2_Hit=%d p2_CR=%d p3_Hit=%d p3_CR=%d n1=%d n2=%d',p1_Hit,p1_CR,p2_Hit,p2_CR,p3_Hit,p3_CR,length(C_Dataofmice),length(E_Dataofmice)));
box off;
set(gcf, 'Renderer', 'Painter'); saveas(gcf,['Welltrained Hit and corr rej rates-' Mani],'fig'); close;

%% d'(discriminality)
c_perf = []; e_perf = []; C_dprime = []; E_dprime = []; 
% control group
for j = 1:size(C_Dataofmice,2) % mice
    if iscell(C_Dataofmice{1,j}.DRTTrialMarker)
        tempTrialMarker = C_Dataofmice{1,j}.DRTTrialMarker{1,WelltrainedDayID};
    else
        tempTrialMarker = C_Dataofmice{1,j}.DRTTrialMarker;
    end
    tempBlockNum = length(tempTrialMarker)/TrialNumPerBlock;
    tempTrialMarker_laseroff = []; tempTrialMarker_laseron = [];
    for k = 1:2:tempBlockNum
        tempTrialMarker_laseroff = [tempTrialMarker_laseroff tempTrialMarker(TrialNumPerBlock*(k-1)+1:TrialNumPerBlock*k)];
    end
    for k = 2:2:tempBlockNum
        tempTrialMarker_laseron = [tempTrialMarker_laseron tempTrialMarker(TrialNumPerBlock*(k-1)+1:TrialNumPerBlock*k)];
    end
    tempGoNum_laseroff = nnz(tempTrialMarker_laseroff==1 | tempTrialMarker_laseroff==2);
    tempNogoNum_laseroff = nnz(tempTrialMarker_laseroff==3 | tempTrialMarker_laseroff==4);
    tempGoNum_laseron = nnz(tempTrialMarker_laseron==1 | tempTrialMarker_laseron==2);
    tempNogoNum_laseron = nnz(tempTrialMarker_laseron==3 | tempTrialMarker_laseron==4);
    tempHitRate_laseroff = nnz(tempTrialMarker_laseroff==1)/tempGoNum_laseroff;
    tempFalseRate_laseroff = nnz(tempTrialMarker_laseroff==3)/tempNogoNum_laseroff;
    tempHitRate_laseron = nnz(tempTrialMarker_laseron==1)/tempGoNum_laseron;
    tempFalseRate_laseron = nnz(tempTrialMarker_laseron==3)/tempNogoNum_laseron;
    if tempHitRate_laseroff == 1
        a = 1-1/(2*tempGoNum_laseroff);
    elseif tempHitRate_laseroff == 0
        a = 1/(2*tempGoNum_laseroff);
    else
        a = norminv(tempHitRate_laseroff);
    end
    if tempFalseRate_laseroff == 1
        b = 1-1/(2*tempNogoNum_laseroff);
    elseif tempFalseRate_laseroff == 0
        b = 1/(2*tempNogoNum_laseroff);
    else
        b = norminv(tempFalseRate_laseroff);
    end
    if tempHitRate_laseron == 1
        c = 1-1/(2*tempGoNum_laseron);
    elseif tempHitRate_laseron == 0
        c = 1/(2*tempGoNum_laseron);
    else
        c = norminv(tempHitRate_laseron);
    end
    if tempFalseRate_laseron == 1
        d = 1-1/(2*tempNogoNum_laseron);
    elseif tempFalseRate_laseron == 0
        d = 1/(2*tempNogoNum_laseron);
    else
        d = norminv(tempFalseRate_laseron);
    end
    C_dprime(j,1) = a-b; % laser off
    C_dprime(j,2) = c-d; % laser on
end
% experimental group
for j = 1:size(E_Dataofmice,2) % Mice
    if iscell(E_Dataofmice{1,j}.DRTTrialMarker)
        tempTrialMarker = E_Dataofmice{1,j}.DRTTrialMarker{1,WelltrainedDayID};
    else
        tempTrialMarker = E_Dataofmice{1,j}.DRTTrialMarker;
    end
    tempBlockNum = length(tempTrialMarker) / TrialNumPerBlock;
    tempTrialMarker_laseroff = []; tempTrialMarker_laseron = [];
    for k = 1:2:tempBlockNum
        tempTrialMarker_laseroff = [tempTrialMarker_laseroff tempTrialMarker(TrialNumPerBlock*(k-1)+1:TrialNumPerBlock*k)];
    end
    for k = 2:2:tempBlockNum
        tempTrialMarker_laseron = [tempTrialMarker_laseron tempTrialMarker(TrialNumPerBlock*(k-1)+1:TrialNumPerBlock*k)];
    end
    tempGoNum_laseroff = nnz(tempTrialMarker_laseroff==1 | tempTrialMarker_laseroff==2);
    tempNogoNum_laseroff = nnz(tempTrialMarker_laseroff==3 | tempTrialMarker_laseroff==4);
    tempGoNum_laseron = nnz(tempTrialMarker_laseron==1 | tempTrialMarker_laseron==2);
    tempNogoNum_laseron = nnz(tempTrialMarker_laseron==3 | tempTrialMarker_laseron==4);
    tempHitRate_laseroff = nnz(tempTrialMarker_laseroff==1)/tempGoNum_laseroff;
    tempFalseRate_laseroff = nnz(tempTrialMarker_laseroff==3)/tempNogoNum_laseroff;
    tempHitRate_laseron = nnz(tempTrialMarker_laseron==1)/tempGoNum_laseron;
    tempFalseRate_laseron = nnz(tempTrialMarker_laseron==3)/tempNogoNum_laseron;
    if tempHitRate_laseroff == 1
        a = 1-1/(2*tempGoNum_laseroff);
    elseif tempHitRate_laseroff == 0
        a = 1/(2*tempGoNum_laseroff);
    else
        a = norminv(tempHitRate_laseroff);
    end
    if tempFalseRate_laseroff == 1
        b = 1-1/(2*tempNogoNum_laseroff);
    elseif tempFalseRate_laseroff == 0
        b = 1/(2*tempNogoNum_laseroff);
    else
        b = norminv(tempFalseRate_laseroff);
    end
    if tempHitRate_laseron == 1
        c = 1-1/(2*tempGoNum_laseron);
    elseif tempHitRate_laseron == 0
        c = 1/(2*tempGoNum_laseron);
    else
        c = norminv(tempHitRate_laseron);
    end
    if tempFalseRate_laseron == 1
        d = 1-1/(2*tempNogoNum_laseron);
    elseif tempFalseRate_laseron == 0
        d = 1/(2*tempNogoNum_laseron);
    else
        d = norminv(tempFalseRate_laseron);
    end
    E_dprime(j,1) = a-b; % laser off
    E_dprime(j,2) = c-d; % laser on
end
PlotWellTrainedResults(C_dprime,E_dprime,0.1,0.3,OptoGroup,0,0,0);
ylabel('Performance ( d prime )','FontSize',18,'FontName','Arial');
box off;

%% Wilcoxon signed-rank test and Mann-Whitney U test
p1_dprime=signrank(C_dprime(:,1),C_dprime(:,2));
p2_dprime=signrank(E_dprime(:,1),E_dprime(:,2));
diffinCtrl = C_dprime(:,2)-C_dprime(:,1);
diffinExp = E_dprime(:,2)-E_dprime(:,1);
p3_dprime = ranksum(diffinCtrl,diffinExp);
title(sprintf('p1=%d p2=%d p3=%d',p1_dprime,p2_dprime,p3_dprime));
set(gcf, 'Renderer', 'Painter'); saveas(gcf,['Welltrained d prime_' Mani],'fig'); close;

%% Ratio of aborted trials
c_aborted = []; e_aborted = [];
for j = 1:size(C_Dataofmice,2) % control
    if length(C_Dataofmice{1,j}.AbortedTrialsRatio_laseron) > 1
        c_aborted(j,1) = C_Dataofmice{1,j}.AbortedTrialsRatio_laseroff(WelltrainedDayID); % laser off
        c_aborted(j,2) = C_Dataofmice{1,j}.AbortedTrialsRatio_laseron(WelltrainedDayID); % laser on
    else
        c_aborted(j,1) = C_Dataofmice{1,j}.AbortedTrialsRatio_laseroff; % laser off
        c_aborted(j,2) = C_Dataofmice{1,j}.AbortedTrialsRatio_laseron; % laser on
    end
end
for j = 1:size(E_Dataofmice,2) % experimental
    if length(E_Dataofmice{1,j}.AbortedTrialsRatio_laseron) > 1
        e_aborted(j,1) = E_Dataofmice{1,j}.AbortedTrialsRatio_laseroff(WelltrainedDayID); % laser off
        e_aborted(j,2) = E_Dataofmice{1,j}.AbortedTrialsRatio_laseron(WelltrainedDayID); % laser on
    else
        e_aborted(j,1) = E_Dataofmice{1,j}.AbortedTrialsRatio_laseroff; % laser off
        e_aborted(j,2) = E_Dataofmice{1,j}.AbortedTrialsRatio_laseron; % laser on
    end
end
PlotWellTrainedResults(c_aborted,e_aborted,0.1,0.3,OptoGroup,0,0,0);
p1_aborted=signrank(c_aborted(:,1),c_aborted(:,2));
p2_aborted=signrank(e_aborted(:,1),e_aborted(:,2));
diffinCtrl = c_aborted(:,2)-c_aborted(:,1);
diffinExp = e_aborted(:,2)-e_aborted(:,1);
p3_aborted = ranksum(diffinCtrl,diffinExp);
title(sprintf('p1=%d p2=%d p3=%d',p1_aborted,p2_aborted,p3_aborted));
ylabel('Ratio of aborted trials ( % )','FontSize',18,'FontName','Arial');
ylim([0 0.5]);
box off;
set(gcf, 'Renderer', 'Painter'); saveas(gcf,['Welltrained Ratio of aborted trials-' Mani],'fig'); close;
save('TwoGroupsMiceIDandPValue','C_Filename','E_Filename');
