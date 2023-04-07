%% Lick rate and lick raster, following optogenetic manipulation.

clear; clc; close all;

%% Load files
% control group
[C_Filename Pathway] = uigetfile({'*.mat','Matlab files(*.mat)';},...
    'Pick some file','MultiSelect','on');
for iMouse = 1:size(C_Filename,2)
    C_Dataofmice{iMouse} = load(C_Filename{1,iMouse});
end
% experimental group
[E_Filename Pathway] = uigetfile({'*.mat','Matlab files(*.mat)';},...
    'Pick some file','MultiSelect','on');
for iMouse = 1:size(E_Filename,2)
    E_Dataofmice{iMouse} = load(E_Filename{1,iMouse});
end

%% Assignment
Beforetrial = 2; SampLen = 1; DelayLen = 10; TestLen = 0.5; RespLen = 1; Aftertrial = 3;
Mani = 'Activate mPFC-aAIC';
C1 = [0 0 0];
OptoGroup = 'ChR2';
if strcmp(OptoGroup,'NpHR')
    C2 = [0 125 0]/255;
else
    C2 = [67 106 178]/255;
end
LearningDay = 6; 
TimeGain = 10;
RespLickNumThres_firsthalf = 2; % go trials
RespLickNumThres = 6; % go trials

%% Lick rate and lick raster
% control group
LickPSTH_Go_C = cell(1,LearningDay);
LickPSTH_NoGo_C = cell(1,LearningDay);
LickRaster_Go_C = cell(1,LearningDay);
LickRaster_NoGo_C = cell(1,LearningDay);
for iSess = 1:LearningDay  % session
    for j = 1:size(C_Dataofmice,2) % mice
        if iSess <= size(C_Dataofmice{1,j}.HitRate,2) % task days
            tempTrialMarker = C_Dataofmice{1,j}.DRTTrialMarker{iSess};
            tempLickPSTH = C_Dataofmice{1,j}.LickPSTH{iSess};
            tempLickRaster = C_Dataofmice{1,j}.LickinTrial{iSess};
            tempTestLickNum = cellfun(@(x) nnz(x>Beforetrial+SampLen+DelayLen & x<Beforetrial+SampLen+DelayLen+TestLen),tempLickRaster,'UniformOutput',1);
            tempRespLickNum_firsthalf = cellfun(@(x) nnz(x>Beforetrial+SampLen+DelayLen+TestLen & x<Beforetrial+SampLen+DelayLen+TestLen+RespLen/2),tempLickRaster,'UniformOutput',1);
            tempRespLickNum = cellfun(@(x) nnz(x>Beforetrial+SampLen+DelayLen+TestLen & x<Beforetrial+SampLen+DelayLen+TestLen+RespLen),tempLickRaster,'UniformOutput',1);
            % lick rate
            tempLickPSTH_Go = tempLickPSTH((tempTrialMarker(:)==1 | tempTrialMarker(:) == 2) & tempRespLickNum_firsthalf(:) >= RespLickNumThres_firsthalf & tempRespLickNum(:) <= RespLickNumThres,:);
            tempLickPSTH_NoGo = tempLickPSTH((tempTrialMarker(:)==3 | tempTrialMarker(:)==4) & tempRespLickNum(:) <= RespLickNumThres,:);
            % lick raster
            tempLickRaster_Go = tempLickRaster(:,(tempTrialMarker(:)==1 | tempTrialMarker(:) == 2) & tempRespLickNum_firsthalf(:) >= RespLickNumThres_firsthalf & tempRespLickNum(:) <= RespLickNumThres);
            tempLickRaster_NoGo = tempLickRaster(:,(tempTrialMarker(:)==3 | tempTrialMarker(:)==4) & tempRespLickNum(:) <= RespLickNumThres);
            if ~isempty(tempLickPSTH_Go)
                LickPSTH_Go_C{iSess} = [LickPSTH_Go_C{iSess}; TimeGain*mean(tempLickPSTH_Go,1)];
                LickRaster_Go_C{iSess} = [LickRaster_Go_C{iSess}; {tempLickRaster_Go}];
            end
            if ~isempty(tempLickPSTH_NoGo)
                LickPSTH_NoGo_C{iSess} = [LickPSTH_NoGo_C{iSess}; TimeGain*mean(tempLickPSTH_NoGo,1)];
                LickRaster_NoGo_C{iSess} = [LickRaster_NoGo_C{iSess}; {tempLickRaster_NoGo}];
            end
        end
    end
end
% experimental group
LickPSTH_Go_E = cell(1,LearningDay);
LickPSTH_NoGo_E = cell(1,LearningDay);
LickRaster_Go_E = cell(1,LearningDay);
LickRaster_NoGo_E = cell(1,LearningDay);
for iSess = 1:LearningDay  % session
    for j = 1:size(E_Dataofmice,2) % mice
        if iSess <= size(E_Dataofmice{1,j}.HitRate,2) % task days
            tempTrialMarker = E_Dataofmice{1,j}.DRTTrialMarker{iSess};
            tempLickPSTH = E_Dataofmice{1,j}.LickPSTH{iSess};
            tempLickRaster = E_Dataofmice{1,j}.LickinTrial{iSess};
            tempTestLickNum = cellfun(@(x) nnz(x>Beforetrial+SampLen+DelayLen & x<Beforetrial+SampLen+DelayLen+TestLen),tempLickRaster,'UniformOutput',1);
            tempRespLickNum_firsthalf = cellfun(@(x) nnz(x>Beforetrial+SampLen+DelayLen+TestLen & x<Beforetrial+SampLen+DelayLen+TestLen+RespLen/2),tempLickRaster,'UniformOutput',1);
            tempRespLickNum = cellfun(@(x) nnz(x>Beforetrial+SampLen+DelayLen+TestLen & x<Beforetrial+SampLen+DelayLen+TestLen+RespLen),tempLickRaster,'UniformOutput',1);
            % lick rate
            tempLickPSTH_Go = tempLickPSTH((tempTrialMarker(:)==1 | tempTrialMarker(:) == 2) & tempRespLickNum_firsthalf(:) >= RespLickNumThres_firsthalf & tempRespLickNum(:) <= RespLickNumThres,:);
            tempLickPSTH_NoGo = tempLickPSTH((tempTrialMarker(:)==3 | tempTrialMarker(:)==4) & tempRespLickNum(:) <= RespLickNumThres,:);
            % lick raster
            tempLickRaster_Go = tempLickRaster(:,(tempTrialMarker(:)==1 | tempTrialMarker(:) == 2) & tempRespLickNum_firsthalf(:) >= RespLickNumThres_firsthalf & tempRespLickNum(:) <= RespLickNumThres);
            tempLickRaster_NoGo = tempLickRaster(:,(tempTrialMarker(:)==3 | tempTrialMarker(:)==4) & tempRespLickNum(:) <= RespLickNumThres);
            if ~isempty(tempLickPSTH_Go)
                LickPSTH_Go_E{iSess} = [LickPSTH_Go_E{iSess}; TimeGain*mean(tempLickPSTH_Go,1)];
                LickRaster_Go_E{iSess} = [LickRaster_Go_E{iSess}; {tempLickRaster_Go}];
            end
            if ~isempty(tempLickPSTH_NoGo)
                LickPSTH_NoGo_E{iSess} = [LickPSTH_NoGo_E{iSess}; TimeGain*mean(tempLickPSTH_NoGo,1)];
                LickRaster_NoGo_E{iSess} = [LickRaster_NoGo_E{iSess}; {tempLickRaster_NoGo}];
            end
        end
    end
end

%% Licking rate in Go and NoGo trials
for iTrialtype = 1:2
    if iTrialtype == 1
        trialname = 'Go trials';
        C_LickPSTH = LickPSTH_Go_C;
        E_LickPSTH = LickPSTH_Go_E;
    elseif iTrialtype == 2
        trialname = 'No-Go trials';
        C_LickPSTH = LickPSTH_NoGo_C;
        E_LickPSTH = LickPSTH_NoGo_E;
    end
    WholeTrialLen = Beforetrial + SampLen + DelayLen + TestLen + RespLen + Aftertrial;
    for iTargetDay = 1:LearningDay
        interval = ceil(WholeTrialLen/6);
        XStartNum = 0;
        figure('OuterPosition',[219 303 750 600]);
        plotshadow(C_LickPSTH{iTargetDay},C1,2,3,0);
        plotshadow(E_LickPSTH{iTargetDay},C2,2,3,0);
        hold on
        time = 1:size(C_LickPSTH{iTargetDay},2);
        plot(time/TimeGain,smooth(mean(C_LickPSTH{iTargetDay},1),3)','color',C1,'linewidth',3);
        hold on
        plot(time/TimeGain,smooth(mean(E_LickPSTH{iTargetDay},1),3)','color',C2,'linewidth',3);
        maxlickrate = 1.5*max(horzcat(smooth(mean(C_LickPSTH{iTargetDay},1),3)',smooth(mean(E_LickPSTH{iTargetDay},1),3)'));
        % cluster-based permuatation test
        [SigTime_Real, SigTimebelowChance_Real] = ClusterBasedPermutation_ForBothReal(C_LickPSTH{iTargetDay},E_LickPSTH{iTargetDay});
        if ~isempty(SigTime_Real)
            for i=1:length(SigTime_Real)
                patch([SigTime_Real(i) SigTime_Real(i) SigTime_Real(i)+1 SigTime_Real(i)+1]./TimeGain,[maxlickrate-0.1 maxlickrate maxlickrate maxlickrate-0.1],C1,'edgecolor','none');
                hold on
            end
        end
        if ~isempty(SigTimebelowChance_Real)
            for i=1:length(SigTimebelowChance_Real)
                patch([SigTimebelowChance_Real(i) SigTimebelowChance_Real(i) SigTimebelowChance_Real(i)+1 SigTimebelowChance_Real(i)+1]./TimeGain,[maxlickrate-0.1 maxlickrate maxlickrate maxlickrate-0.1],C1,'edgecolor','none');
                hold on
            end
        end
        set(gca,'XTick',Beforetrial:interval:WholeTrialLen,'XTickLabel',{XStartNum,num2str(XStartNum+interval),...
            num2str(XStartNum+2*interval),num2str(XStartNum+3*interval),num2str(XStartNum+4*interval)},'FontName','Arial','FontSize',16,'xlim',[1 WholeTrialLen]);
        set(gca,'YTick',[0 2 4 6 8],'YTickLabel',{'0','2','4','6','8'},'FontName','Arial','FontSize',16,'ylim',[0 maxlickrate]);
        xlabel('Time (s)','FontSize',18,'FontName','Arial');
        ylabel('Lick Rate (Hz)','FontSize',18,'FontName','Arial');
        plot([Beforetrial Beforetrial],[0 8],'linestyle','--','color',C1);
        plot([Beforetrial+SampLen  Beforetrial+SampLen],[0 8],'linestyle','--','color',C1);
        plot([Beforetrial+SampLen+DelayLen Beforetrial+SampLen+DelayLen],[0 8],'linestyle','--','color',C1);
        plot([Beforetrial+SampLen+DelayLen+TestLen Beforetrial+SampLen+DelayLen+TestLen],[0 8],'linestyle','--','color',C1);
        box off;
        set(gcf, 'Renderer', 'Painter'); saveas(gcf,['Lick Rate in two groups-' trialname '-in Day' num2str(iTargetDay) '-' Mani],'fig');
        close;
    end
end

%% Licking raster in Go and NoGo trials
ShownTrialNum = 10; ExamDayID_Go = 1; ExamDayID_NoGo = 5; 
for iTrialtype = 1:2
    if iTrialtype == 1
        tempLickRaster_Go_C = LickRaster_Go_C{ExamDayID_Go};
        tempLickRaster_Go_E = LickRaster_Go_E{ExamDayID_Go};
        tempMiceNum_C = numel(tempLickRaster_Go_C);
        tempMiceNum_E = numel(tempLickRaster_Go_E);
        tempPermuMiceID_C = randperm(tempMiceNum_C);
        tempPermuMiceID_E = randperm(tempMiceNum_E);
        tempLickRaster_Go_C = tempLickRaster_Go_C{tempPermuMiceID_C(1)};
        tempLickRaster_Go_E = tempLickRaster_Go_E{tempPermuMiceID_E(1)};
        tempTrialNum_C = numel(tempLickRaster_Go_C);
        tempTrialNum_E = numel(tempLickRaster_Go_E);
        tempPermuTrialID_C = randperm(tempTrialNum_C);
        tempPermuTrialID_E = randperm(tempTrialNum_E);
        ExamLickRaster_Go_C = tempLickRaster_Go_C(:,tempPermuTrialID_C(1:ShownTrialNum));
        ExamLickRaster_Go_E = tempLickRaster_Go_E(:,tempPermuTrialID_E(1:ShownTrialNum));
        trial = 'Go trials';
        if strcmp(OptoGroup,'ChR2')
            Lickraster = horzcat(ExamLickRaster_Go_E,ExamLickRaster_Go_C);
        else
            Lickraster = horzcat(ExamLickRaster_Go_C,ExamLickRaster_Go_E);
        end
        TarDayID = ExamDayID_Go;
    else
        tempLickRaster_NoGo_C = LickRaster_NoGo_C{ExamDayID_NoGo};
        tempLickRaster_NoGo_E = LickRaster_NoGo_E{ExamDayID_NoGo};
        tempMiceNum_C = numel(tempLickRaster_NoGo_C);
        tempMiceNum_E = numel(tempLickRaster_NoGo_E);
        tempPermuMiceID_C = randperm(tempMiceNum_C);
        tempPermuMiceID_E = randperm(tempMiceNum_E);
        tempLickRaster_NoGo_C = tempLickRaster_NoGo_C{tempPermuMiceID_C(1)};
        tempLickRaster_NoGo_E = tempLickRaster_NoGo_E{tempPermuMiceID_E(1)};
        tempTrialNum_C = numel(tempLickRaster_NoGo_C);
        tempTrialNum_E = numel(tempLickRaster_NoGo_E);
        tempPermuTrialID_C = randperm(tempTrialNum_C);
        tempPermuTrialID_E = randperm(tempTrialNum_E);
        ExamLickRaster_NoGo_C = tempLickRaster_NoGo_C(:,tempPermuTrialID_C(1:ShownTrialNum));
        ExamLickRaster_NoGo_E = tempLickRaster_NoGo_E(:,tempPermuTrialID_E(1:ShownTrialNum));
        trial = 'No-Go trials';
        if strcmp(OptoGroup,'ChR2')
            Lickraster = horzcat(ExamLickRaster_NoGo_E,ExamLickRaster_NoGo_C);
        else
            Lickraster = horzcat(ExamLickRaster_NoGo_C,ExamLickRaster_NoGo_E);
        end
        TarDayID = ExamDayID_NoGo;
    end
    figure('OuterPosition',[219 303 420 534]);
    patch([Beforetrial Beforetrial Beforetrial+SampLen Beforetrial+SampLen],[0 length(Lickraster)+3 length(Lickraster)+3 0],'k','FaceAlpha',0.2,'edgecolor','none');
    hold on
    patch([Beforetrial+SampLen+DelayLen Beforetrial+SampLen+DelayLen Beforetrial+SampLen+DelayLen+TestLen Beforetrial+SampLen+DelayLen+TestLen],[0 length(Lickraster)+3 length(Lickraster)+3 0],'k','FaceAlpha',0.2,'edgecolor','none');
    hold on
    for iTrial = 1:numel(Lickraster) % trial ID
        if strcmp(OptoGroup,'ChR2')
            if iTrial <= ShownTrialNum
                for iLick = 1:size(Lickraster{iTrial},1) % experimental
                    plot([Lickraster{iTrial}(iLick,1) Lickraster{iTrial}(iLick,1)], [length(Lickraster)+2-iTrial-0.5 length(Lickraster)+2-iTrial+0.5],'color',C2);
                    hold on
                end
            else
                for iLick = 1:size(Lickraster{iTrial},1) % control
                    plot([Lickraster{iTrial}(iLick,1) Lickraster{iTrial}(iLick,1)], [length(Lickraster)+1-iTrial-0.5 length(Lickraster)+1-iTrial+0.5],'color',C1);
                    hold on
                end
            end
        else
            if iTrial <= ShownTrialNum
                for iLick = 1:size(Lickraster{iTrial},1) % control
                    plot([Lickraster{iTrial}(iLick,1) Lickraster{iTrial}(iLick,1)], [length(Lickraster)+2-iTrial-0.5 length(Lickraster)+2-iTrial+0.5],'color',C1);
                    hold on
                end
            else
                for iLick = 1:size(Lickraster{iTrial},1) % experimental
                    plot([Lickraster{iTrial}(iLick,1) Lickraster{iTrial}(iLick,1)], [length(Lickraster)+1-iTrial-0.5 length(Lickraster)+1-iTrial+0.5],'color',C2);
                    hold on
                end
            end
        end
    end
    set(gca,'XTick',0:1:Beforetrial+SampLen+DelayLen+TestLen+RespLen+3,'XTickLabel',{'0','5','10','15'},'FontName','Arial','FontSize',16,'xlim',[Beforetrial-1 Beforetrial+SampLen+DelayLen+TestLen+RespLen+3]);
    set(gca,'YTick',[1 ShownTrialNum length(Lickraster)+1],'YTickLabel',{'','',''},'FontName','Arial','FontSize',16,'ylim',[0 length(Lickraster)+2-1+0.5]);
    box off;
    set(gcf, 'Renderer', 'Painter'); saveas(gcf,['LickRasterBetweenTwoGroups-' trial '-LearningDay-' num2str(TarDayID) '-' Mani],'fig'); close all;
end

