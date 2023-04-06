%% Plot OGNG behavioral performance in the well-trained phase.

clear; clc; close all;

%% Enter path
CurrPath = uigetdir;
cd(CurrPath);

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
Task = 'OGNG';
Mani = 'Suppress mPFC';
OptoGroup = 'ChR2'; 
C1 = [0 0 0]; 
if strcmp(OptoGroup,'NpHR') 
    C2 = [0 125 0]/255;
else
    C2 = [67 106 178]/255;
end
WelltrainedDayID = 1; TrialNumPerBlock = 20;

%% Hit rate and CR rate
HitRate_Ctrl = []; 
HitRate_Exp = [];
CRRate_Ctrl = []; 
CRRate_Exp = [];
for j = 1:size(C_Dataofmice,2) % control group
    if strcmp(class(C_Dataofmice{1,1}.AllDayBlockHitRate),'cell')
        TotalBlockNum = size(C_Dataofmice{1,j}.AllDayBlockHitRate{1,WelltrainedDayID},2);
        HitRate_Ctrl(j,1) =  mean(C_Dataofmice{1,j}.AllDayBlockHitRate{1,WelltrainedDayID}(:,1:2:TotalBlockNum-1)); % laser off
        CRRate_Ctrl(j,1) =  mean(C_Dataofmice{1,j}.AllDayBlockCRRate{1,WelltrainedDayID}(:,1:2:TotalBlockNum-1)); 
        HitRate_Ctrl(j,2) =  mean(C_Dataofmice{1,j}.AllDayBlockHitRate{1,WelltrainedDayID}(:,2:2:TotalBlockNum)); % laser on
        CRRate_Ctrl(j,2) =  mean(C_Dataofmice{1,j}.AllDayBlockCRRate{1,WelltrainedDayID}(:,2:2:TotalBlockNum));
    else
        TotalBlockNum = size(C_Dataofmice{1,j}.AllDayBlockHitRate,2);
        HitRate_Ctrl(j,1) =  mean(C_Dataofmice{1,j}.AllDayBlockHitRate(:,1:2:TotalBlockNum-1)); % laser off
        CRRate_Ctrl(j,1) =  mean(C_Dataofmice{1,j}.AllDayBlockCRRate(:,1:2:TotalBlockNum-1));
        HitRate_Ctrl(j,2) =  mean(C_Dataofmice{1,j}.AllDayBlockHitRate(:,2:2:TotalBlockNum)); % laser on
        CRRate_Ctrl(j,2) =  mean(C_Dataofmice{1,j}.AllDayBlockCRRate(:,2:2:TotalBlockNum));
    end
end
for j = 1:size(E_Dataofmice,2) % experimental group
    if strcmp(class(E_Dataofmice{1,1}.AllDayBlockHitRate),'cell')
        TotalBlockNum = size(E_Dataofmice{1,j}.AllDayBlockHitRate{1,WelltrainedDayID},2);
        HitRate_Exp(j,1) =  mean(E_Dataofmice{1,j}.AllDayBlockHitRate{1,WelltrainedDayID}(:,1:2:TotalBlockNum-1)); % laser off
        CRRate_Exp(j,1) =  mean(E_Dataofmice{1,j}.AllDayBlockCRRate{1,WelltrainedDayID}(:,1:2:TotalBlockNum-1)); 
        HitRate_Exp(j,2) =  mean(E_Dataofmice{1,j}.AllDayBlockHitRate{1,WelltrainedDayID}(:,2:2:TotalBlockNum)); % laser on
        CRRate_Exp(j,2) =  mean(E_Dataofmice{1,j}.AllDayBlockCRRate{1,WelltrainedDayID}(:,2:2:TotalBlockNum));
    else
        TotalBlockNum = size(E_Dataofmice{1,j}.AllDayBlockHitRate,2);
        HitRate_Expl(j,1) =  mean(E_Dataofmice{1,j}.AllDayBlockHitRate(:,1:2:TotalBlockNum-1)); % laser off
        CRRateExp(j,1) =  mean(E_Dataofmice{1,j}.AllDayBlockCRRate(:,1:2:TotalBlockNum-1));
        HitRate_Exp(j,2) =  mean(E_Dataofmice{1,j}.AllDayBlockHitRate(:,2:2:TotalBlockNum)); % laser on
        CRRate_Exp(j,2) =  mean(E_Dataofmice{1,j}.AllDayBlockCRRate(:,2:2:TotalBlockNum));
    end
end
PlotWellTrainedResults(HitRate_Ctrl,HitRate_Exp,0.1,0.3,OptoGroup,0,0,1);
hold on
PlotWellTrainedResults(CRRate_Ctrl,CRRate_Exp,0.1,0.3,OptoGroup,0,0,0);
hold on

%% Wilcoxon signed-rank test and Mann-Whitney U test
Phit_wic = signrank(HitRate_Ctrl(:,1),HitRate_Ctrl(:,2)); % hit rate
Phit_wie = signrank(HitRate_Exp(:,1),HitRate_Exp(:,2));
Diff_HitRate_Ctrl = HitRate_Ctrl(:,2)-HitRate_Ctrl(:,1);
Diff_HitRate_Exp = HitRate_Exp(:,2)-HitRate_Exp(:,1);
Phit_ce = ranksum(Diff_HitRate_Ctrl,Diff_HitRate_Exp); 
Pcr_wic = signrank(CRRate_Ctrl(:,1),CRRate_Ctrl(:,2)); % corr. rej. rate
Pcr_wie = signrank(CRRate_Exp(:,1),CRRate_Exp(:,2));
Diff_CRRate_Ctrl = CRRate_Ctrl(:,2)-CRRate_Ctrl(:,1);
Diff_CRRate_Exp = CRRate_Exp(:,2)-CRRate_Exp(:,1);
Pcr_ce = ranksum(Diff_CRRate_Ctrl,Diff_CRRate_Exp); 
title(sprintf('Phit_wic=%d Phit_wie=%d Phit_ce=%d Pcr_wic=%d Pcr_wie=%d Pcr_ce=%d',Phit_wic,Phit_wie,Phit_ce,Pcr_wic,Pcr_wie,Pcr_ce));
ylabel('Correct Rate ( % )','FontSize',18,'FontName','Arial');
ylim([0 1]);
box off;
set(gcf, 'Renderer', 'Painter'); saveas(gcf,['Hit and CR rates-' Mani '-WellTrained-' Task],'fig'); close;
