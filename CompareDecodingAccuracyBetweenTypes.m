%% Compare decoding accuracy of FCSP events between congruent active and incongruent inactive FC neuronal pairs.

clear; clc; close all;

%% Assignment
CurrPath = uigetdir;
File_CongAct = dir('*CongAct-FCSPsDecoding*.mat');
File_IncongInact = dir('*InconInact-FCSPsDecoding*.mat');

%% Load decoding result
% congruent active FC neuronal pairs
DecodingAccuracy_CongAct = zeros(100,size(File_CongAct,1));
for iFile = 1:size(File_CongAct,1)
    tempFileName = File_CongAct(iFile).name;
    StartPos = regexp(tempFileName,'=');
    EndPos = regexp(tempFileName,'-CtrlGroup');
    tempPos = StartPos+1:EndPos-1;
    tempPairNum = str2num(tempFileName(tempPos))-1;
    load(tempFileName);
    DecodingAccuracy_CongAct(:,tempPairNum) = data;
end
% incongruent inactive FC neuronal pairs
DecodingAccuracy_IncongInact = zeros(100,size(File_IncongInact,1));
for iFile = 1:size(File_IncongInact,1)
    tempFileName = File_IncongInact(iFile).name;
    StartPos = regexp(tempFileName,'=');
    EndPos = regexp(tempFileName,'-CtrlGroup');
    tempPos = StartPos+1:EndPos-1;
    tempPairNum = str2num(tempFileName(tempPos))-1;
    load(tempFileName);
    DecodingAccuracy_IncongInact(:,tempPairNum) = data;
end

%% Figure
figure('position',[500 200 500 300]);
plotshadow(DecodingAccuracy_CongAct,[1 0 0],3,5,0.1);
plotshadow(DecodingAccuracy_IncongInact,[0 0 1],3,5,0.1);
% bootstrap test
UnitsNum = 2:60;
for iNum = 1:size(DecodingAccuracy_CongAct,2)
    tempdata_CongAct = DecodingAccuracy_CongAct(:,iNum);
    tempdata_IncongInact = DecodingAccuracy_IncongInact(:,iNum);
    tempDiff = tempdata_CongAct - tempdata_IncongInact;
    if mean(tempDiff) >= 0
        tempP = nnz(tempDiff<=0)/numel(tempDiff);
    else
        tempP = nnz(tempDiff>=0)/numel(tempDiff);
    end
    if tempP < 0.05
        patch([UnitsNum(iNum)/10-0.05 UnitsNum(iNum)/10-0.05 UnitsNum(iNum)/10+0.05 UnitsNum(iNum)/10+0.05],[0.53 0.55 0.55 0.53],[0 0 0],'edgecolor','none'); hold on
    end
end
set(gca,'XTick',0:1:UnitsNum(end)/10,'xlim',[0 UnitsNum(end)/10],'YTick',0.5:0.1:1,'ylim',[0.5 1]);
set(gcf,'Render','Painter'); saveas(gcf,'Compare decoding accuracy of CongAct and InconInact FC','fig'); close all;









