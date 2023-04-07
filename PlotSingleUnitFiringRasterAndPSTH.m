%% Plot single case unit firing raster and PSTH 

clear; clc; close all;

%% Assignment
TargetBrain = 'aAIC';
UnitID = 93;
Group = 'CtrlGroup';
TrialNum = 10;
C1 = [1 0 0];
C2 = [0 0 1]; 
TimeGain = 10;

%% Load results of spike rasters, FR, and selectivity of target neurons
FileName_FR = strcat(TargetBrain,'FRinS1S2_',Group,'.mat');
load(FileName_FR);
FileName_Selec = strcat(TargetBrain,'SelectivityData_',Group,'.mat');
load(FileName_Selec);

%% Spike rasters in example trials for example neuron
RandomTrialID = randperm(size(TargetBrainRGinS1{1,1},1));
RandomTrialID = RandomTrialID(1,1:TrialNum);
SingleUnitRGinS1 = TargetBrainRGinS1{1,UnitID}(RandomTrialID,1);
SingleUnitRGinS2 = TargetBrainRGinS2{1,UnitID}(RandomTrialID,1);
SingleUnitRG =  [SingleUnitRGinS1; SingleUnitRGinS2];

%% Raster plot
figure('OuterPosition',[219 303 420 534]);
patch([1 1 2 2],[0 size(SingleUnitRG,1)+3 size(SingleUnitRG,1)+3 0],'k','FaceAlpha',0.2,'edgecolor','none'); hold on
patch([6 6 6.5 6.5],[0 size(SingleUnitRG,1)+3 size(SingleUnitRG,1)+3 0],'k','FaceAlpha',0.2,'edgecolor','none'); hold on
for itr = 1:size(SingleUnitRG,1) % trial
    if itr <= TrialNum
        for itr1 = 1:size(SingleUnitRG{itr,1},1) % spike rasters in S1 trials
            plot([SingleUnitRG{itr,1}(itr1,1)-3 SingleUnitRG{itr,1}(itr1,1)-3], [size(SingleUnitRG,1)+2-itr-0.5 size(SingleUnitRG,1)+2-itr+0.5],'r');
            hold on
        end
    else
        for itr1 = 1:size(SingleUnitRG{itr,1},1) % spike rasters in S2 trials
            plot([SingleUnitRG{itr,1}(itr1,1)-3 SingleUnitRG{itr,1}(itr1,1)-3], [size(SingleUnitRG,1)+1-itr-0.5 size(SingleUnitRG,1)+1-itr+0.5],'b');
            hold on
        end
    end
end
set(gca,'XTick',1:1:9,'XTickLabel',num2cell(0:1:8),'FontName','Arial','FontSize',16,'xlim',[0.1 7]);
set(gca,'YTick',[0 TrialNum size(SingleUnitRG,1)],'YTickLabel',{'','',''},'FontName','Arial','FontSize',16,'ylim',[0 size(SingleUnitRG,1)+2-1+0.5]);
box off
set(gcf,'Renderer','Painter'); saveas(gcf,['Raster_' FileName_FR(1:4) '_Unit ID' num2str(UnitID)],'fig'); close all;

%% PSTH plot
figure('OuterPosition',[219 303 420 534]);
TargetBrainUnitsFRinS1 = TargetBrainUnitsFRinS1{UnitID}(:,11:end);
TargetBrainUnitsFRinS2 = TargetBrainUnitsFRinS2{UnitID}(:,11:end);
plot((1:size(TargetBrainUnitsFRinS1,2))/TimeGain, smooth(mean(TargetBrainUnitsFRinS1,1),3)', 'color',C1,'linewidth',2);
hold on
plot((1:size(TargetBrainUnitsFRinS2,2))/TimeGain, smooth(mean(TargetBrainUnitsFRinS2,1),3)', 'color',C2,'linewidth',2);
hold on
patch([1 1 2 2],[0 5 5 0],'k','FaceAlpha',0.2,'edgecolor','none');
hold on
patch([6 6 6.5 6.5],[0 5 5 0],'k','FaceAlpha',0.2,'edgecolor','none');
hold on
plotshadow(TargetBrainUnitsFRinS1,C1,2,3,0);
hold on
plotshadow(TargetBrainUnitsFRinS2,C2,2,3,0);
if ~isempty(TargetBrainUnitsSigBinID{UnitID})
    hold on
    for j=1:length(TargetBrainUnitsSigBinID{UnitID})
        patch([TargetBrainUnitsSigBinID{UnitID}(j)-1 TargetBrainUnitsSigBinID{UnitID}(j)-1 TargetBrainUnitsSigBinID{UnitID}(j) TargetBrainUnitsSigBinID{UnitID}(j)],[0.5 1 1 0.5],[0 0 0],'edgecolor','none');
        hold on
    end
end
set(gca,'XTick',1:1:9,'XTickLabel',num2cell(0:1:8),'FontName','Arial','FontSize',16,'xlim',[0.1 7]);
box off;
set(gcf,'Renderer','Painter'); saveas(gcf,['PSTH_' FileName_FR(1:4) '_Unit ID' num2str(UnitID)],'fig'); close all;