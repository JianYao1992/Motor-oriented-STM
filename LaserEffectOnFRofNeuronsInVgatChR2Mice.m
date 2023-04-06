%% This code is used for analysing laser effect on FR of mPFC neurons in Vgat-ChR2 mice.

clear; clc; close all;

%% Enter path
pwd = uigetdir;
CurrentPath=pwd; AllPath=genpath(CurrentPath); SplitPath=strsplit(AllPath,';'); SubPath=SplitPath'; SubPath=SubPath(2:end-1);

%% Assignment
TarMiceID = {'M20'};
beforelaser = 9;
laserduration = 4; 
afterlaser = 9; 
plus = 2;
LaserBin = 10:13;
TimeGain=10;
LeftmPFCUnitsRG = [];
LeftmPFCUnitsFR = [];
RightmPFCUnitsRG = [];
RightmPFCUnitsFR = [];
LeftmPFCUnitsWaveform = [];
RightmPFCUnitsWaveform = [];

%% FR of target neurons
for iPath = 1:size(SubPath,1) 
    IsTargetPath = [];
    for k = 1:length(TarMiceID)
        IsTargetPath = [IsTargetPath regexp(SubPath{iPath,1},TarMiceID{1,k})];
    end
    if ~isempty(IsTargetPath)
        Path=SubPath{iPath,1};
        cd(Path);
        JAVAFiles = dir('*.mat');
        for j = 1:size(JAVAFiles,1)
            load(JAVAFiles(j,1).name);
            if ~isempty(SingleUnitList)
                for iUnit = 1:size(SingleUnitList,1) % Unit
                    if SingleUnitList(iUnit,1) <= 8
                        LeftmPFCUnitsRG = [LeftmPFCUnitsRG {(LaserRGResults(iUnit,1:size(LaserRGResults,2)/2))'}];
                        tempSingleUnitFR = [];
                        for iTrial = 1:size(LaserResults,2)/2 % Trial
                            tempSingleUnitFR = [tempSingleUnitFR; LaserResults{1,iTrial}(iUnit,:)];
                        end
                        LeftmPFCUnitsFR = [LeftmPFCUnitsFR {tempSingleUnitFR}];
                        LeftmPFCUnitsWaveform = [LeftmPFCUnitsWaveform; NewWaveForm(iUnit,:)];
                    else
                        RightmPFCUnitsRG = [RightmPFCUnitsRG {(LaserRGResults(iUnit,1+size(LaserRGResults,2)/2:end))'}];
                        tempSingleUnitFR = [];
                        for iTrial = 1+size(LaserResults,2)/2:length(LaserResults) % Trial
                            tempSingleUnitFR = [tempSingleUnitFR; LaserResults{1,iTrial}(iUnit,:)];
                        end
                        RightmPFCUnitsFR = [RightmPFCUnitsFR {tempSingleUnitFR}];
                        RightmPFCUnitsWaveform = [RightmPFCUnitsWaveform; NewWaveForm(iUnit,:)];
                    end
                end
            end
        end
    end
end
cd(pwd);
mPFCUnitsFR = [LeftmPFCUnitsFR RightmPFCUnitsFR];
mPFCUnitsRG = [LeftmPFCUnitsRG RightmPFCUnitsRG];
mPFCUnitsWaveform = [LeftmPFCUnitsWaveform; RightmPFCUnitsWaveform];

%% Define ID of neurons showing activated,inhibited, and unchanged modulation by laser
SuppressedUnitID = [];
ActivatedUnitID = [];
for iUnit = 1:size(mPFCUnitsFR,2) 
    mPFCUnitsFR{1,iUnit} = TimeGain*mPFCUnitsFR{1,iUnit};
    for iBin = LaserBin
        p(iUnit,find(LaserBin==iBin)) = ranksum(mean(mPFCUnitsFR{1,iUnit}(:,1+10*(iBin-1):10*iBin),2),mean(mPFCUnitsFR{1,iUnit}(:,TimeGain*(beforelaser-1)+1:beforelaser*TimeGain),2));
    end
    corrected_p(iUnit,:) = p(iUnit,:)*length(LaserBin);
end
for iUnit = 1:size(mPFCUnitsFR,2)
    for iBin=LaserBin(1) % the first delay second
        if corrected_p(iUnit,find(LaserBin==iBin)) <= 0.05
            if mean(mean(mPFCUnitsFR{1,iUnit}(:,1+10*(iBin-1):10*iBin),2)) > mean(mean(mPFCUnitsFR{1,iUnit}(:,1+(beforelaser-1)*TimeGain:beforelaser*TimeGain),2)) & isempty(find(ActivatedUnitID==iUnit))
                ActivatedUnitID = [ActivatedUnitID iUnit];
            elseif mean(mean(mPFCUnitsFR{1,iUnit}(:,1+10*(iBin-1):10*iBin),2)) < mean(mean(mPFCUnitsFR{1,iUnit}(:,1+(beforelaser-1)*TimeGain:beforelaser*TimeGain),2)) & isempty(find(SuppressedUnitID==iUnit))
                SuppressedUnitID = [SuppressedUnitID iUnit];
            end
        end
    end
end
ModulatedUnitID = [ActivatedUnitID SuppressedUnitID];
UnchangedUnitID = setdiff(1:length(mPFCUnitsFR),ModulatedUnitID);

%% Proportion of suppressed and activated neurons
colormap([1 0 0; 0 0 1; 0 0 0]);
pie(horzcat(SuppressedUnitID,ActivatedUnitID,UnchangedUnitID));
set(gcf,'Render','Painter'); saveas(gcf,'Distribution of responsive aAIC neurons','fig');

%% Plot z-scored FR of suppressed and activated neurons
% activation
Z_on = []; % laser on period
ActivatedUnitsFR = mPFCUnitsFR(1,ActivatedUnitID);
for itr=1:size(ActivatedUnitsFR,2) % unit
    Z_on(itr,:) = (mean(ActivatedUnitsFR{itr})-mean(mean(ActivatedUnitsFR{itr}(:,1+(beforelaser-2)*TimeGain:beforelaser*TimeGain),2)))/(0.00000000001+std(mean(ActivatedUnitsFR{itr}(:,1+(beforelaser-2)*TimeGain:beforelaser*TimeGain),1)));
end
higher_on = (smooth(mean(Z_on,1),3))' + std(Z_on,0,1)/sqrt(size(Z_on,1)); higher_on = higher_on((beforelaser-2)*TimeGain+1:(beforelaser+laserduration+plus)*TimeGain);
lower_on = (smooth(mean(Z_on,1),3))' - std(Z_on,0,1)/sqrt(size(Z_on,1));  lower_on = lower_on((beforelaser-2)*TimeGain+1:(beforelaser+laserduration+plus)*TimeGain);
Z_off = []; % laser off period
ActivatedUnitsFR = mPFCUnitsFR(1,ActivatedUnitID);
for itr=1:size(ActivatedUnitsFR,2) % unit
    Z_off(itr,:) = (mean(ActivatedUnitsFR{itr})-mean(mean(ActivatedUnitsFR{itr}(:,1:2*TimeGain),2)))/(0.00000000001+std(mean(ActivatedUnitsFR{itr}(:,1:2*TimeGain),1)));
end
higher_off = (smooth(mean(Z_off,1),3))' + std(Z_off,0,1)/sqrt(size(Z_off,1)); higher_off = higher_off(1:(2+laserduration+plus)*TimeGain);
lower_off = (smooth(mean(Z_off,1),3))' - std(Z_off,0,1)/sqrt(size(Z_off,1)); lower_off = lower_off(1:(2+laserduration+plus)*TimeGain);
figure;
x = 1:length(higher_on); % time bin
X = [x/10,fliplr(x/10)];
Y_off = [higher_off, fliplr(lower_off)]; Y_on = [higher_on,fliplr(lower_on)];
d = fill(X,Y_on,[255 100 18]/255,'edgecolor','none');
alpha(d,0.5);
hold on
plot(x/10,smooth(mean(Z_on(:,(beforelaser-2)*TimeGain+1:(beforelaser+laserduration+plus)*TimeGain),1),3),'color',[1 0 0],'linewidth',4);
hold on
e = fill(X,Y_off,[0 0 0],'edgecolor','none');
alpha(e,0.5);
hold on
plot(x/10,smooth(mean(Z_off(:,1:(2+laserduration+plus)*TimeGain),1),3),'color',[0 0 0],'linewidth',4);
hold on
MaxNormFR = 1.1*max(higher_on); MinNormFR = 0.9*min(lower_on);
rectangle('Position',[2,MaxNormFR-1,laserduration,0.75],'LineStyle','none','faceColor',[0 0 1]);
box('off');
set(gca,'XTick',[2 2+length(LaserBin)],'XTickLabel',{'0','4'},'FontName','Arial','xlim',[1 2+length(LaserBin)+plus-1]);
set(gca,'YTick',0:5:15,'YTickLabel',{'0','5','10','15'},'FontName','Arial','FontSize',16,'ylim',[MinNormFR MaxNormFR]);
hold on
plot([2 2],[MinNormFR MaxNormFR],'k--','linewidth',2);
hold on
plot([2+length(LaserBin) 2+length(LaserBin)],[MinNormFR MaxNormFR],'k--','linewidth',2);
ylabel('Normalized Firing Rate','FontSize',16,'FontName','Arial');
set(gcf,'Renderer','Painter'); saveas(gcf,'ZscoreFiringRate_ActivatedUnit','fig'); close;
% suppression
Z_on = []; % laser on period
SuppressedUnitsFR = mPFCUnitsFR(1,SuppressedUnitID);
for itr=1:size(SuppressedUnitsFR,2) % unit
    Z_on(itr,:)=(mean(SuppressedUnitsFR{itr})-mean(mean(SuppressedUnitsFR{itr}(:,1+(beforelaser-2)*TimeGain:beforelaser*TimeGain),2)))/(0.00000000001+std(mean(SuppressedUnitsFR{itr}(:,1+(beforelaser-2)*TimeGain:beforelaser*TimeGain),1)));
end
higher_on = (smooth(mean(Z_on,1),3))' + std(Z_on,0,1)/sqrt(size(Z_on,1)); higher_on = higher_on((beforelaser-2)*TimeGain+1:(beforelaser+laserduration+plus)*TimeGain);
lower_on = (smooth(mean(Z_on,1),3))' - std(Z_on,0,1)/sqrt(size(Z_on,1));  lower_on = lower_on((beforelaser-2)*TimeGain+1:(beforelaser+laserduration+plus)*TimeGain);
Z_off = []; % laser off period
SuppressedUnitsFR = mPFCUnitsFR(1,SuppressedUnitID);
for itr = 1:size(SuppressedUnitsFR,2) % unit
    Z_off(itr,:) = (mean(SuppressedUnitsFR{itr})-mean(mean(SuppressedUnitsFR{itr}(:,1:2*TimeGain),2)))/(0.00000000001+std(mean(SuppressedUnitsFR{itr}(:,1:2*TimeGain),1)));
end
higher_off = (smooth(mean(Z_off,1),3))' + std(Z_off,0,1)/sqrt(size(Z_off,1)); higher_off = higher_off(1:(2+laserduration+plus)*TimeGain);
lower_off = (smooth(mean(Z_off,1),3))' - std(Z_off,0,1)/sqrt(size(Z_off,1)); lower_off = lower_off(1:(2+laserduration+plus)*TimeGain);
figure;
x = 1:length(higher_on); % time bin
X = [x/10,fliplr(x/10)];
Y_off = [higher_off,fliplr(lower_off)]; Y_on = [higher_on,fliplr(lower_on)];
d=fill(X,Y_on,[67 106 178]/255,'edgecolor','none');
alpha(d,0.5);
hold on
plot(x/10,smooth(mean(Z_on(:,(beforelaser-2)*TimeGain+1:(beforelaser+laserduration+plus)*TimeGain),1),3),'color',[67 106 178]/255,'linewidth',4);
hold on
e=fill(X,Y_off,[0 0 0],'edgecolor','none');
alpha(e,0.5);
hold on
plot(x/10,smooth(mean(Z_off(:,1:(2+laserduration+plus)*TimeGain),1),3),'color',[0 0 0],'linewidth',4);
hold on
MaxNormFR = 1.1*max(higher_on); MinNormFR = 1.1*min(lower_on);
rectangle('Position',[2,MaxNormFR-1,laserduration,0.75],'LineStyle','none','faceColor',[0 0 1]);
box('off');
set(gca,'XTick',[2 2+length(LaserBin)],'XTickLabel',{'0','4'},'FontName','Arial','xlim',[1 2+length(LaserBin)+plus-1]);
set(gca,'YTick',-4:2:2,'YTickLabel',{'-4','-2','0','2'},'FontName','Arial','FontSize',16,'ylim',[MinNormFR MaxNormFR]);
hold on
plot([2 2],[MinNormFR MaxNormFR],'k--','linewidth',2);
hold on
plot([2+length(LaserBin) 2+length(LaserBin)],[MinNormFR MaxNormFR],'k--','linewidth',2);
ylabel('Normalized Firing Rate','FontSize',16,'FontName','Arial'); 
set(gcf,'Renderer','Painter'); saveas(gcf,'ZscoreFiringRate_SuppressedUnit','fig'); close;

%% Plot spike rasters of target neurons
trialnum = 10;
for iUnit = 1:size(mPFCUnitsRG,2)
    tempUnitsRG = cell(2*trialnum,1);
    SelectedTrial = randperm(size(mPFCUnitsRG{1,iUnit},1)); SelectedTrial(trialnum+1:end-trialnum) = [];
    figure;
    for itr = 1:2*trialnum % laser off and on trials
        if itr <= trialnum
            timeoffset = 2;
        else
            timeoffset = 8;
        end
        for itr1 = 1:size(mPFCUnitsRG{1,iUnit}{SelectedTrial(itr),1},1) % each RG
            plot([mPFCUnitsRG{1,iUnit}{SelectedTrial(itr),1}(itr1,1)-timeoffset mPFCUnitsRG{1,iUnit}{SelectedTrial(itr),1}(itr1,1)-timeoffset], [2*trialnum+1-itr-0.5 2*trialnum+1-itr+0.5],'k');
            hold on
        end
    end
    xlim([0 7]);
    ylim([0.5 2*trialnum+0.5]);
    box off;
    saveas(gcf,['mPFCUnit#' num2str(iUnit) 'Raster'],'fig');
    close;
end

%% PSTH plots of target neurons
Baseline = 1:90; 
LaserOn = 71:160;
for iUnit = 1:size(mPFCUnitsFR,2)
    tempUnitsFR = mPFCUnitsFR{iUnit};
    figure;
    x = 1:90;
    % laser off
    highervalue_off = smooth(mean(tempUnitsFR(:,Baseline),1),3)' + std(tempUnitsFR(:,Baseline),0,1)/sqrt(size(tempUnitsFR,1));
    lowervalue_off = smooth(mean(tempUnitsFR(:,Baseline),1),3)' - std(tempUnitsFR(:,Baseline),0,1)/sqrt(size(tempUnitsFR,1));
    % laser on
    highervalue_on = smooth(mean(tempUnitsFR(:,LaserOn),1),3)' + std(tempUnitsFR(:,LaserOn),0,1)/sqrt(size(tempUnitsFR,1));
    lowervalue_on = smooth(mean(tempUnitsFR(:,LaserOn),1),3)' - std(tempUnitsFR(:,LaserOn),0,1)/sqrt(size(tempUnitsFR,1));
    X = [x/10,fliplr(x/10)];
    Y_off = [highervalue_off, fliplr(lowervalue_off)]; Y_on = [highervalue_on, fliplr(lowervalue_on)];
    d=fill(X,Y_off,[0 0 0],'edgecolor','none');
    alpha(d,0.5);
    hold on
    plot(x/10,smooth(mean(tempUnitsFR(:,Baseline),1),3)','k','linewidth',5);
    hold on
    e=fill(X,Y_on,[67 106 178]/255,'edgecolor','none');
    alpha(e,0.5);
    hold on
    plot(x/10,smooth(mean(tempUnitsFR(:,LaserOn),1),3)','color',[67 106 178]/255,'linewidth',5);
    set(gca,'XTick',1:2+laserduration+plus,'XTickLabel',{'0','1','2','3','4','5','6','7'},'FontName','Arial','xlim',[1 2+laserduration+plus]);
    box off; saveas(gcf,['mPFCUnit#' num2str(iUnit) 'FiringRate'],'fig'); close all;
end

%% calculate averaged waveform of target neurons
for itr1 = 1:size(mPFCUnitsWaveform,1) % neuron
    figure;
    time = 1:size(mPFCUnitsWaveform,2); % time
    for itr2 = 1:4
        step = size(mPFCUnitsWaveform,2)/4;
        plot(time((1+(itr2-1)*step):(itr2*step)),mPFCUnitsWaveform(itr1,(1+(itr2-1)*step):(itr2*step)),'color',[0 0 0],'linewidth',3);
        hold on
    end
    box('off');
    hold on
    plot([210 210+58*1000/1450],[-100 -100],'k');
    hold on
    plot([210 210],[-100 -50],'k');
    axis off;
    saveas(gcf,['mPFCWaveFormUnit#' num2str(itr1)],'fig');
    close;
end






