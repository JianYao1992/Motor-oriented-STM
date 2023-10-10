%% This code plots spatial spread of blue laser effect on Vgat-ChR2 + mice.

clear; clc; close all;

%% Assignment
LaserPowerNum = 3;
Dist = [0 1 2 3];
TimeGain = 10;
C = [{[0 0 0]},{[0 0 1]},{[1 0 0]}];

%% Enter pathway
CurrPath = uigetdir;

%% Assignment
BeforeLaserLen = 9;
ShownBaselineLen = 4;
LaserLen = 4;

%% Information of all neurons
MatFile = dir('*.mat');
AllUnitsFR = cell(0);
AllUnitsDistID = [];
AllUnitsTypeID = [];
for iFile = 1:size(MatFile,1)
    load(MatFile(iFile,1).name);
    for iUnit = 1:size(SingleUnitList,1)
        tempUnitFR = cellfun(@(x) x(iUnit,:), LaserResults,'UniformOutput',0);
        AllUnitsFR{end+1,1} = vertcat(tempUnitFR{:});
    end
    AllUnitsDistID = [AllUnitsDistID; DistID];
    AllUnitsTypeID = [AllUnitsTypeID; UnitTypeID];
end

%% Averaged FR of neurons in each combination of distance and laser power
LaserTrialNum = numel(LaserResults);
PerLaserPowerTrialNum = LaserTrialNum/LaserPowerNum;
% putative pyramidal neurons
AverFR_Pyr = cell(LaserPowerNum,numel(Dist));
for itrDist = 1:numel(Dist)
    TarUnitsFR = AllUnitsFR(AllUnitsDistID==Dist(itrDist) & AllUnitsTypeID==1,:);
    for itrLaser = 1:LaserPowerNum
        temp = cellfun(@(x) mean(x(1+PerLaserPowerTrialNum*(itrLaser-1):PerLaserPowerTrialNum*itrLaser,:),1)*TimeGain,TarUnitsFR,'UniformOutput',0);
        AverFR_Pyr{itrLaser,itrDist} = vertcat(temp{:});
    end
end
AverFR_Pyr = cellfun(@(x) mean(x(:,1+TimeGain*BeforeLaserLen:TimeGain*(BeforeLaserLen+LaserLen)),2)./mean(x(:,1+TimeGain*(BeforeLaserLen-ShownBaselineLen):TimeGain*BeforeLaserLen),2),AverFR_Pyr,'UniformOutput',0);
% putative FS neurons
AverFR_Fs = cell(LaserPowerNum,numel(Dist));
for itrDist = 1:numel(Dist)
    TarUnitsFR = AllUnitsFR(AllUnitsDistID==Dist(itrDist) & AllUnitsTypeID==2,:);
    for itrLaser = 1:LaserPowerNum
        temp = cellfun(@(x) mean(x(1+PerLaserPowerTrialNum*(itrLaser-1):PerLaserPowerTrialNum*itrLaser,:),1)*TimeGain,TarUnitsFR,'UniformOutput',0);
        AverFR_Fs{itrLaser,itrDist} = vertcat(temp{:});
    end
end
AverFR_Fs = cellfun(@(x) mean(x(:,1+TimeGain*BeforeLaserLen:TimeGain*(BeforeLaserLen+LaserLen)),2)./mean(x(:,1+TimeGain*(BeforeLaserLen-ShownBaselineLen):TimeGain*BeforeLaserLen),2),AverFR_Fs,'UniformOutput',0);

%% Plot error bar
% putative pyramidal neurons
figure('position',[400 200 550 400]);
for iLaser = 1:size(AverFR_Pyr,1)
    tempFR = AverFR_Pyr(iLaser,:);
    tempAverFR = cellfun(@mean,tempFR);
    tempSem = cellfun(@(x) std(x)/sqrt(size(x,1)),tempFR,'UniformOutput',1);
    errorbar(1:1:size(tempFR,2),tempAverFR,tempSem,'color',C{iLaser},'marker','o','markerfacecolor',C{iLaser},'markeredgecolor','none'); hold on;
end;
yline(1,'--k'); hold on
box off
set(gca,'XTick',1:4,'XTickLabel',num2cell(0:3),'XLim',[0.7 size(tempFR,2)+0.3]);
set(gcf,'Render','Painter'); saveas(gcf,'Spatial spread of optogenetic effect on putative pyramidal neurons','fig'); close all;
% putative FS neurons
figure('position',[400 200 550 400]);
for iLaser = 1:size(AverFR_Fs,1)
    tempFR = AverFR_Fs(iLaser,:);
    tempAverFR = cellfun(@mean,tempFR);
    tempSem = cellfun(@(x) std(x)/sqrt(size(x,1)),tempFR,'UniformOutput',1);
    errorbar(1:1:size(tempFR,2),tempAverFR,tempSem,'color',C{iLaser},'marker','o','markerfacecolor',C{iLaser},'markeredgecolor','none'); hold on;
end;
yline(1,'--k'); hold on
box off
set(gca,'XTick',1:4,'XTickLabel',num2cell(0:3),'XLim',[0.7 size(tempFR,2)+0.3]);
set(gcf,'Render','Painter'); saveas(gcf,'Spatial spread of optogenetic effect on putative FS neurons','fig'); close all;

%% Tw-ANOVA-md
% putative pyramidal neurons
Frame_Pyr = ConstructTwANOVAmdTestMatrix(AverFR_Pyr);
[SSQs, DFs, MSQs, Fs, Ps] = mixed_between_within_anova(Frame_Pyr,0,'Statistics for spatial spread of optogenetic effect on putative pyramidal neurons_Tw-ANOVA-md');
% putative pyramidal neurons
Frame_Fs = ConstructTwANOVAmdTestMatrix(AverFR_Fs);
[SSQs, DFs, MSQs, Fs, Ps] = mixed_between_within_anova(Frame_Fs,0,'Statistics for spatial spread of optogenetic effect on putative FS neurons_Tw-ANOVA-md');
















