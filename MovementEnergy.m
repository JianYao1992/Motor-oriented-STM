%% Movement energy in Go and NoGo trials, including torso movement, nose movement, and pupil diameter.

clear; clc; close all;

%% Assignment
ShownBaseLen = 0.5;
SampOdorLen = 1;
DelayLen = 4;
TestOdorLen = 0.5;
SmoothValue = 3;
TimeGain = 60;
MovementResults_Trunk = cell(1,2);
MovementResults_Nose = cell(1,2);
MovementResults_Pupil = cell(1,2);
C1 = [1 0 0];
C2 = [0 0 1];

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

%% Analyze the movement probability in Go and NoGo trials
for iPath = 1:size(SubPath,1)
    Path = SubPath{iPath,1};
    cd(Path);
    % extracellular recording file
    DataFiles = dir('*SumExtraCellularAndVideo*.mat');
    for iFile = 1:size(DataFiles,1)
        % load extracellular data
        if size(DataFiles,1) == 1
            filename = DataFiles.name;
        else
            filename = DataFiles(iFile,1).name;
        end
        load(filename);
        % trunk movement
        if ~isempty(AverTrunkMove_Go)
            SingleSessMove_Trunk_Go = AverTrunkMove_Go;
            SingleSessMove_Trunk_Go = SingleSessMove_Trunk_Go./mean(SingleSessMove_Trunk_Go(:,1:TimeGain*ShownBaseLen),2);
            SingleSessMove_Trunk_NoGo = AverTrunkMove_NoGo;
            SingleSessMove_Trunk_NoGo = SingleSessMove_Trunk_NoGo./mean(SingleSessMove_Trunk_NoGo(:,1:TimeGain*ShownBaseLen),2);
        end
        % nose movement
        if ~isempty(AverNoseMove_Go)
            SingleSessMove_Nose_Go = AverNoseMove_Go;
            SingleSessMove_Nose_Go = SingleSessMove_Nose_Go./mean(SingleSessMove_Nose_Go(:,1:TimeGain*ShownBaseLen),2);
            SingleSessMove_Nose_NoGo = AverNoseMove_NoGo;
            SingleSessMove_Nose_NoGo = SingleSessMove_Nose_NoGo./mean(SingleSessMove_Nose_NoGo(:,1:TimeGain*ShownBaseLen),2);
        end
        % pupil diameter
        if ~isempty(AverPupilDilation_Go)
            SingleSessMove_Pupil_Go = AverPupilDilation_Go;
            SingleSessMove_Pupil_Go = SingleSessMove_Pupil_Go./mean(SingleSessMove_Pupil_Go(:,1:TimeGain*ShownBaseLen),2);
            SingleSessMove_Pupil_NoGo = AverPupilDilation_NoGo;
            SingleSessMove_Pupil_NoGo = SingleSessMove_Pupil_NoGo./mean(SingleSessMove_Pupil_NoGo(:,1:TimeGain*ShownBaseLen),2);
        end
        
        %% Intergrated into all-session results
        % torso
        if exist('SingleSessMove_Trunk_Go') && isempty(regexp(filename,'Day01'))
            MovementResults_Trunk{1} = [MovementResults_Trunk{1}; SingleSessMove_Trunk_Go];
            MovementResults_Trunk{2} = [MovementResults_Trunk{2}; SingleSessMove_Trunk_NoGo];
        end
        % nose
        if exist('SingleSessMove_Nose_Go') && isempty(regexp(filename,'Day01'))
            MovementResults_Nose{1} = [MovementResults_Nose{1}; SingleSessMove_Nose_Go];
            MovementResults_Nose{2} = [MovementResults_Nose{2}; SingleSessMove_Nose_NoGo];
        end
        % pupil
        if exist('SingleSessMove_Pupil_Go') && isempty(regexp(filename,'Day01'))
            MovementResults_Pupil{1} = [MovementResults_Pupil{1}; SingleSessMove_Pupil_Go];
            MovementResults_Pupil{2} = [MovementResults_Pupil{2}; SingleSessMove_Pupil_NoGo];
        end
        clear AverTrunkMove_Go AverTrunkMove_NoGo AverNoseMove_Go AverNoseMove_NoGo AverPupilDilation_Go AverPupilDilation_NoGo...
            SingleSessMove_Trunk_Go SingleSessMove_Trunk_NoGo SingleSessMove_Nose_Go SingleSessMove_Nose_NoGo SingleSessMove_Pupil_Go SingleSessMove_Pupil_NoGo
    end
end
%% Compare movement in Go and NoGo trials
cd(CurrPath);
% trunk
PlotGoNoGoTrialsMovement('Trunk movement',MovementResults_Trunk{1},MovementResults_Trunk{2},ShownBaseLen,SampOdorLen,DelayLen,TestOdorLen,TimeGain,SmoothValue,C1,C2);
% nose
PlotGoNoGoTrialsMovement('Nose movement',MovementResults_Nose{1},MovementResults_Nose{2},ShownBaseLen,SampOdorLen,DelayLen,TestOdorLen,TimeGain,SmoothValue,C1,C2);
% pupil
PlotGoNoGoTrialsMovement('Pupil dilation',MovementResults_Pupil{1},MovementResults_Pupil{2},ShownBaseLen,SampOdorLen,DelayLen,TestOdorLen,TimeGain,SmoothValue,C1,C2);



function PlotGoNoGoTrialsMovement(ParameterName,Result_S1,Result_S2,ShownBaseLen,SampOdorLen,DelayLen,TestOdorLen,TimeGain,SmoothValue,Color1,Color2)

figure('OuterPosition',[200 200 500 300]);
plotshadow(Result_S1,Color1,2,SmoothValue,0,TimeGain);
plotshadow(Result_S2,Color2,2,SmoothValue,0,TimeGain);
MaxValue = 1.5 * max(horzcat(smooth(mean(Result_S1,1),SmoothValue)',smooth(mean(Result_S2,1),SmoothValue)'));
% cluster-based permutation test
[SigTime_Real, SigTimebelowChance_Real] = ClusterBasedPermutation_ForBothReal(Result_S1,Result_S2);
if ~isempty(SigTime_Real)
    for i=1:length(SigTime_Real)
        patch([SigTime_Real(i) SigTime_Real(i) SigTime_Real(i)+1 SigTime_Real(i)+1]./TimeGain,[MaxValue-0.1 MaxValue MaxValue MaxValue-0.1],[0 0 0],'edgecolor','none');
        hold on
    end
end
if ~isempty(SigTimebelowChance_Real)
    for i=1:length(SigTimebelowChance_Real)
        patch([SigTimebelowChance_Real(i) SigTimebelowChance_Real(i) SigTimebelowChance_Real(i)+1 SigTimebelowChance_Real(i)+1]./TimeGain,[MaxValue-0.1 MaxValue MaxValue MaxValue-0.1],[0 0 0],'edgecolor','none');
        hold on
    end
end
plot([ShownBaseLen ShownBaseLen],[0 MaxValue],'--k');
plot([ShownBaseLen+SampOdorLen ShownBaseLen+SampOdorLen],[0 MaxValue],'--k');
plot([ShownBaseLen+SampOdorLen+DelayLen ShownBaseLen+SampOdorLen+DelayLen],[0 MaxValue],'--k');
plot([ShownBaseLen+SampOdorLen+DelayLen+TestOdorLen ShownBaseLen+SampOdorLen+DelayLen+TestOdorLen],[0 MaxValue],'--k');
set(gca,'XTick',0:1:ShownBaseLen+SampOdorLen+DelayLen+TestOdorLen,'XTickLabel',num2cell((0:1:ShownBaseLen+SampOdorLen+DelayLen+TestOdorLen)-ShownBaseLen),'FontName','Arial','FontSize',16,'xlim',[0 ShownBaseLen+SampOdorLen+DelayLen+TestOdorLen]);
set(gca,'YLim',[0 MaxValue]);
xlabel('Time from sample onset (s)','FontSize',18,'FontName','Arial');
ylabel([ParameterName  ' (%)'],'FontSize',18,'FontName','Arial');
title(sprintf('n = %d sessions',size(Result_S1,1)));
box off;
set(gcf, 'Renderer', 'Painter'); saveas(gcf,sprintf('%s in Go and NoGo trials',ParameterName),'fig'); close all;
end
