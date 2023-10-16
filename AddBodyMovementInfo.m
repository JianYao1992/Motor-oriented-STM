%% Information about body movement, including pupil size, nose movement, and torso movement.

clear; clc; close all;

%% Assignment
ShownBaseLen = 0.5;
SampOdorLen = 1;
DelayLen = 4;
TestOdorLen = 0.5;
RespLen = 0;
ShownITIlen = 0;
TimeGain = 60;
AnalysisWindow = 0.2;
BinNumOfAnalysisWindow = AnalysisWindow*TimeGain;
BinNumThresForSigAnalysisWindow = (AnalysisWindow/4)*TimeGain;

%% Target directory
CurrPath = uigetdir;
AllPath = genpath(CurrPath);
SplitPath = strsplit(AllPath,';');
SubPath = SplitPath';
SubPath = SubPath(2:end-1);

%% Remove first-class pathway
FirstGradeFolder = dir(CurrPath);
ToRemoveFileID = [];
for i = 3:size(FirstGradeFolder,1)
    % pathway name to remove
    ToDeleteFolderName = strcat(CurrPath,'\',FirstGradeFolder(i).name);
    % pathway ID to remove
    for j = 1:length(SubPath)
        if strcmp(ToDeleteFolderName,SubPath(j)) == 1
            ToRemoveFileID = [ToRemoveFileID; j];
        end
    end
end
% remove pathway unrelated to data
SubPath(ToRemoveFileID)= [];

%% Analyze the significance of body movement
for iPath = 1:size(SubPath,1)
    Path = SubPath{iPath,1};
    cd(Path);
    % extracellular recording file
    ExtraCellularFiles = dir('*short*.mat');
    % torso movement file
    TrunkFiles = dir('*_Normed1505BMbT_trunk*.mat');
    % nose movement file
    NoseFiles = dir('*_Normed1505BMbT_nose*.mat');
    % pupil diameter file
    PupilFiles = dir('*_Normed1505BMbT_pupil*.mat');
    
    for iFile = 1:size(ExtraCellularFiles,1)
        % load extracellular data
        if size(ExtraCellularFiles,1) == 1
            filename = ExtraCellularFiles.name;
        elseif size(ExtraCellularFiles,1) > 1
            filename = ExtraCellularFiles(iFile,1).name;
        end
        load(filename);
        fprintf('Is processing file: %s\n',filename);
        % load trunk movement data
        if size(TrunkFiles,1) >= 1
            if size(TrunkFiles,1) == 1
                load(TrunkFiles.name);
            elseif size(TrunkFiles,1) > 1
                load(TrunkFiles(iFile,1).name);
            end
            BMinfo_Trunk = byTrialBM;
            timestamp_Trunk = tempxticks;
        end
        % load nose movement data
        if size(NoseFiles,1) >= 1
            if size(NoseFiles,1) == 1
                load(NoseFiles.name);
            else
                load(NoseFiles(iFile,1).name);
            end
            BMinfo_Nose = byTrialBM;
            timestamp_Nose = tempxticks;
        end
        % load pupil diameter data
        if size(PupilFiles,1) >= 1
            if size(PupilFiles,1) == 1
                load(PupilFiles.name);
            else
                load(PupilFiles(iFile,1).name);
            end
            BMinfo_Pupil = byTrialBM;
            timestamp_Pupil = tempxticks;
        end
        
        %% Judge analysis windows with statistically significant movement for each trial
        % trunk movement
        if exist('BMinfo_Trunk')
            [TrialBasedBM_Trunk,AverTrunkMove_Go,AverTrunkMove_NoGo,IsSigBM_Trunk,TrunkMoveProbability_Go,TrunkMoveProbability_NoGo,IsTrunkMoveWindow_delay,IsTrunkMoveWindow_delay_Go,IsTrunkMoveWindow_delay_NoGo,TimePeriod_Trunk] = JudgeSigMovementWindow('Trunk movement',BMinfo_Trunk,CompletedTrialID,TrialMark,timestamp_Trunk,TimeGain,ShownBaseLen,SampOdorLen,DelayLen,TestOdorLen,RespLen,ShownITIlen,[1 0 0],AnalysisWindow,BinNumOfAnalysisWindow,BinNumThresForSigAnalysisWindow,false);
        else
            TrialBasedBM_Trunk = []; AverTrunkMove_Go = []; AverTrunkMove_NoGo = []; IsSigBM_Trunk = []; TrunkMoveProbability_Go = []; TrunkMoveProbability_NoGo = []; IsTrunkMoveWindow_delay = []; IsTrunkMoveWindow_delay_Go = []; IsTrunkMoveWindow_delay_NoGo = []; TimePeriod_Trunk = [];
        end
        % nose movement
        if exist('BMinfo_Nose')
            [TrialBasedBM_Nose,AverNoseMove_Go,AverNoseMove_NoGo,IsSigBM_Nose,NoseMoveProbability_Go,NoseMoveProbability_NoGo,IsNoseMoveWindow_delay,IsNoseMoveWindow_delay_Go,IsNoseMoveWindow_delay_NoGo,TimePeriod_Nose] = JudgeSigMovementWindow('Nose movement',BMinfo_Nose,CompletedTrialID,TrialMark,timestamp_Nose,TimeGain,ShownBaseLen,SampOdorLen,DelayLen,TestOdorLen,RespLen,ShownITIlen,[1 0 0],AnalysisWindow,BinNumOfAnalysisWindow,BinNumThresForSigAnalysisWindow,false);
        else
            TrialBasedBM_Nose = []; AverNoseMove_Go = []; AverNoseMove_NoGo = []; IsSigBM_Nose = []; NoseMoveProbability_Go = []; NoseMoveProbability_NoGo = []; IsNoseMoveWindow_delay = []; IsNoseMoveWindow_delay_Go = []; IsNoseMoveWindow_delay_NoGo = []; TimePeriod_Nose = [];
        end
        % pupil diameter
        if exist('BMinfo_Pupil')
            [TrialBasedBM_Pupil,AverPupilDilation_Go,AverPupilDilation_NoGo,IsSigBM_Pupil,PupilDilationProbability_Go,PupilDilationProbability_NoGo,IsPupilDilationWindow_delay,IsPupilDilationWindow_delay_Go,IsPupilDilationWindow_delay_NoGo,TimePeriod_Pupil] = JudgeSigMovementWindow('Pupil dilation',BMinfo_Pupil,CompletedTrialID,TrialMark,timestamp_Pupil,TimeGain,ShownBaseLen,SampOdorLen,DelayLen,TestOdorLen,RespLen,ShownITIlen,[1 0 1],AnalysisWindow,BinNumOfAnalysisWindow,BinNumThresForSigAnalysisWindow,false);
        else
            TrialBasedBM_Pupil = []; AverPupilDilation_Go = []; AverPupilDilation_NoGo = []; IsSigBM_Pupil = []; PupilDilationProbability_Go = []; PupilDilationProbability_NoGo = []; IsPupilDilationWindow_delay = []; IsPupilDilationWindow_delay_Go = []; IsPupilDilationWindow_delay_NoGo = []; TimePeriod_Pupil = [];
        end
        %% Save results
        if ~exist('LaserResults')
            LaserResults = [];
            LaserRGResults = [];
        end
        save(strcat('Analysis window number-',num2str(DelayLen/AnalysisWindow),'-SumExtraCellularAndVideo',filename),'AllTrialsNumber','CompletedTrialID','LaserResults','LaserRGResults','LickRate','LickTime','Newdata','NewWaveForm','Results','RGResults','SingleUnitList','SpikeTime','TrialMark',...
            'TrialBasedBM_Trunk','AverTrunkMove_Go','AverTrunkMove_NoGo','IsSigBM_Trunk','TrunkMoveProbability_Go','TrunkMoveProbability_NoGo','IsTrunkMoveWindow_delay','IsTrunkMoveWindow_delay_Go','IsTrunkMoveWindow_delay_NoGo','TimePeriod_Trunk',...
            'TrialBasedBM_Nose','AverNoseMove_Go','AverNoseMove_NoGo','IsSigBM_Nose','NoseMoveProbability_Go','NoseMoveProbability_NoGo','IsNoseMoveWindow_delay','IsNoseMoveWindow_delay_Go','IsNoseMoveWindow_delay_NoGo','TimePeriod_Nose',...
            'TrialBasedBM_Pupil','AverPupilDilation_Go','AverPupilDilation_NoGo','IsSigBM_Pupil','PupilDilationProbability_Go','PupilDilationProbability_NoGo','IsPupilDilationWindow_delay','IsPupilDilationWindow_delay_Go','IsPupilDilationWindow_delay_NoGo','TimePeriod_Pupil','-v7.3');
        clear LaserResults LaserRGResults BMinfo_Trunk BMinfo_Nose BMinfo_Pupil
    end
end


function [TrialBasedBM,AverBM_Go,AverBM_NoGo,IsSigBM,MovementProbability_Go,MovementProbability_NoGo,IsBmWindow_delay,IsBmWindow_delay_Go,IsBmWindow_delay_NoGo,TimePeriod] = JudgeSigMovementWindow(MovementType,BMinfo,FinishTrialID,TrialMark,timestamp,TimeGain,ShownBaseLen,SampOdorLen,DelayLen,TestOdorLen,RespLen,ShownITIlen,ColorValue,AnalysisWindow,BinNumOfAnalysisWindow,BinNumThresForSigAnalysisWindow,ToPlot)

% video data in target time period
if size(BMinfo,1) > size(BMinfo,2)
    TrialBasedBM = BMinfo';
elseif size(BMinfo,1) < size(BMinfo,2)
    TrialBasedBM = BMinfo;
end
SampOnsetTs = find(timestamp==0);
TrialBasedBM = TrialBasedBM(:,SampOnsetTs-floor(TimeGain*ShownBaseLen):SampOnsetTs+floor(TimeGain*(SampOdorLen+DelayLen+TestOdorLen+RespLen+ShownITIlen)));
% ID of Go, NoGo, and aborted trials
TrialID_Go = FinishTrialID(TrialMark(:,3)==1);
TrialID_NoGo = FinishTrialID(TrialMark(:,3)==2);
TrialID_Aborted = setdiff(1:1:size(TrialBasedBM,1),FinishTrialID);
% smooth movement data
for iTrial = 1:size(TrialBasedBM,1)
    TrialBasedBM(iTrial,:) = smooth(TrialBasedBM(iTrial,:),2*floor(TimeGain/10))'; % 200-msec smooth
end
% averaged movement in Go and NoGo trials
AverBM_Go = nanmean(TrialBasedBM(TrialID_Go,:));
AverBM_NoGo = nanmean(TrialBasedBM(TrialID_NoGo,:));
% test significance of movement in each bin, and plot
IsSigBM = [];
for iTrial = 1:size(TrialBasedBM,1)
    tempBM = TrialBasedBM(iTrial,:);
    tempBaseLevel = prctile(tempBM,5);
    tempSD = std(tempBM);
    tempThreshold = tempBaseLevel + 1.5*tempSD;
    tempIsSigBM = tempBM > max([1 tempThreshold]);
    IsSigBM = [IsSigBM; tempIsSigBM];
    if ToPlot
        % plot
        figure('position',[200 200 500 300]);
        patch([floor(TimeGain*ShownBaseLen) floor(TimeGain*ShownBaseLen) floor(TimeGain*(ShownBaseLen+SampOdorLen)) floor(TimeGain*(ShownBaseLen+SampOdorLen))],[0 5 5 0],'k','FaceAlpha',0.2,'edgecolor','none'); hold on
        patch([floor(TimeGain*(ShownBaseLen+SampOdorLen+DelayLen)) floor(TimeGain*(ShownBaseLen+SampOdorLen+DelayLen)) floor(TimeGain*(ShownBaseLen+SampOdorLen+DelayLen+TestOdorLen)) floor(TimeGain*(ShownBaseLen+SampOdorLen+DelayLen+TestOdorLen))],[0 5 5 0],'k','FaceAlpha',0.2,'edgecolor','none'); hold on
        plot(0:1:numel(tempBM)-1,tempBM,'color',ColorValue); hold on
        plot(0:1:numel(tempBM)-1,tempIsSigBM-1,'k'); hold on
        plot([0 numel(tempBM)-1],[tempThreshold tempThreshold],'--k'); hold on
        % judge current trial type
        if ismember(iTrial,TrialID_Go)
            trialtype = 'Go trial';
        elseif ismember(iTrial,TrialID_NoGo)
            trialtype = 'NoGo trial';
        elseif ismember(iTrial,TrialID_Aborted)
            trialtype = 'Aborted trial';
        end
        set(gca,'XTick',0:TimeGain:numel(tempBM)-1,'XTickLabel',num2cell((0:TimeGain:numel(tempBM)-1)/TimeGain),'XLim',[0 numel(tempBM)-1]);
        set(gcf,'Render','Painter'); saveas(gcf,sprintf('%s movement-%s-TrialID %d',MovementType,trialtype,iTrial)); close all;
    end
end
% movement probability of Go trials
MovementProbability_Go = nanmean(IsSigBM(TrialID_Go,:));
MovementProbability_NoGo = nanmean(IsSigBM(TrialID_NoGo,:));
% test significance of movement in each analysis window duirng the delay period
IsSigBM_delay = IsSigBM(:,1+TimeGain*(ShownBaseLen+SampOdorLen):TimeGain*(ShownBaseLen+SampOdorLen+DelayLen));
IsBmWindow_delay = zeros(size(IsSigBM_delay,1),floor(DelayLen/AnalysisWindow));
for iTrial = 1:size(IsSigBM_delay,1)
    IsBM_CurrTrial = IsSigBM_delay(iTrial,:);
    for iWindow = 1:floor(DelayLen/AnalysisWindow)
        IsBM_CurrWindow = IsBM_CurrTrial(1+BinNumOfAnalysisWindow*(iWindow-1):BinNumOfAnalysisWindow*iWindow);
        if nnz(IsBM_CurrWindow) >= BinNumThresForSigAnalysisWindow
            IsBmWindow_delay(iTrial,iWindow) = 1;
        end
    end
end
IsBmWindow_delay_Go = IsBmWindow_delay(TrialID_Go,:);
IsBmWindow_delay_NoGo = IsBmWindow_delay(TrialID_NoGo,:);
TimePeriod = [-1*ShownBaseLen SampOdorLen+DelayLen+TestOdorLen+RespLen+ShownITIlen];
end