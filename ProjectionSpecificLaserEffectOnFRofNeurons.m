%% Laser effect on  FR, following projection-specific manipulation.

clear; clc; close all;

%% Assignment
TargetBrain = 'aAIC'; 
Group  = 'NpHRGroup';
beforeneuronnum = 0; FileID = 1; BilateralmPFCMiceID = [{'M19'} {'M20'}];
AllunitRG=[]; AllunitFR=[]; mPFC_ID = []; aAIC_ID=[];

%% Target path
pwd = uigetdir; 
CurrPath=pwd; AllPath=genpath(CurrPath); SplitPath=strsplit(AllPath,';'); SubPath=SplitPath'; SubPath=SubPath(2:end-1);

%% FR and ID of all neurons
for iPath = 1:size(SubPath,1)
    Path=SubPath{iPath,1};
    cd(Path);
    JAVAFiles = dir('*.mat');
    for j = 1:size(JAVAFiles,1)
        Filename{1,FileID} = JAVAFiles(j,1).name;
        load(Filename{1,FileID});
        if ~isempty(SingleUnitList)
            % FR of all neurons
            for iUnit = 1:size(SingleUnitList,1) % neuron
                AllunitRG = [AllunitRG {(LaserRGResults(iUnit,:))'}];
                tempSingleUnitFR = [];
                for iTrial = 1:size(LaserResults,2) % trial
                    tempSingleUnitFR = [tempSingleUnitFR; LaserResults{1,iTrial}(iUnit,:)];
                end
                AllunitFR = [AllunitFR {tempSingleUnitFR}];
            end
            tempmPFCID = []; tempAIID = [];
            % ID of mPFC neurons
            Position = [];
            for k = 1:length(BilateralmPFCMiceID)
                Position = [Position regexp(Filename{1,FileID},BilateralmPFCMiceID{1,k})]
            end
            if ~isempty(Position)
                tempmPFCID = find(SingleUnitList(:,1)>=1 & SingleUnitList(:,1)<=16);  % mice with bilateral mPFC recording
            else
                tempmPFCID = find((SingleUnitList(:,1)>=9 & SingleUnitList(:,1)<=16) | (SingleUnitList(:,1)>=25 & SingleUnitList(:,1)<=32));
            end
            if ~isempty(tempmPFCID)
                mPFC_ID = [mPFC_ID tempmPFCID'+beforeneuronnum];
            end
            % ID of aAIC neurons
            if ~isempty(Position)
                tempAIID = [];
            else
                tempAIID = find((SingleUnitList(:,1)>=1 & SingleUnitList(:,1)<=8) | (SingleUnitList(:,1)>=17 & SingleUnitList(:,1)<=24));
            end
            if ~isempty(tempAIID)
                aAIC_ID = [aAIC_ID tempAIID'+beforeneuronnum];
            end
            beforeneuronnum = beforeneuronnum + size(SingleUnitList,1);
        end
        FileID = FileID + 1;
    end
end
cd(pwd);
if strcmp(TargetBrain,'mPFC')
    TarUnitsID = mPFC_ID;
else
    TarUnitsID = aAIC_ID;
end
TarUnitsFR = AllunitFR(:,TarUnitsID);
TarUnitsRG = AllunitRG(:,TarUnitsID);
LaserFR = cellfun(@(x) x(:,31:90),TarUnitsFR,'UniformOutput', 0);
beforelaser = 1; laserduration = 4; afterlaser = 1; plus = 0; TimeGain=10;

%% ID of neurons showing suppressed, activated, and unchanged modulation by laser
SuppressedUnitID = [];
ActivatedUnitID = [];
for iUnit = 1:size(LaserFR,2)
    LaserFR{1,iUnit}=TimeGain*LaserFR{1,iUnit};
    for iBin = 1:floor(size(LaserFR{1,iUnit},2)/TimeGain)
        p(iUnit,iBin) = ranksum(mean(LaserFR{1,iUnit}(:,1+10*(iBin-1):10*iBin),2),mean(LaserFR{1,iUnit}(:,1:beforelaser*TimeGain),2));
    end
    corrected_p(iUnit,:) = p(iUnit,:)*4;
end
for iUnit = 1:size(LaserFR,2)
    for iBin = 2 % the first delay second
        if corrected_p(iUnit,iBin) <= 0.05
            if mean(mean(LaserFR{1,iUnit}(:,1+10*(iBin-1):10*iBin),2)) > mean(mean(LaserFR{1,iUnit}(:,1:10),2)) & isempty(find(ActivatedUnitID==iUnit))
                ActivatedUnitID = [ActivatedUnitID iUnit];
            elseif mean(mean(LaserFR{1,iUnit}(:,1+10*(iBin-1):10*iBin),2)) < mean(mean(LaserFR{1,iUnit}(:,1:10),2)) & isempty(find(SuppressedUnitID==iUnit))
                SuppressedUnitID = [SuppressedUnitID iUnit];
            end
        end
    end
end
ModulatedUnitID = [ActivatedUnitID SuppressedUnitID];
UnchangedUnitID = setdiff(1:length(TargetBrainUnitsID),ModulatedUnitID);

%% Proportion of suppressed, activated, and unchanged neurons
figure;
cm = [0 0 1; 0 1 0; [100 100 100]/255];
colormap(cm);
pie([length(ActivatedUnitID),length(SuppressedUnitID),length(UnchangedUnitID)]);
set(gcf,'Render','Painter'); saveas(gcf,['Proportion of responsive ' TargetBrain ' neurons'],'fig'); close all;

%% Plot Z-scored firing rate
% activation
Z = [];
unitdisplay = LaserFR(1,ActivatedUnitID);
for itr = 1:size(unitdisplay,2) % neuron
    Z(itr,:) = (mean(unitdisplay{itr})-mean(mean(unitdisplay{itr}(:,1:10),2)))/(0.00000000001+std(mean(unitdisplay{itr}(:,1:10),1)));
end
figure;
x = 1:size(Z,2); % time bin
higher_act = (smooth(mean(Z,1),3))' + std(Z,0,1)/sqrt(size(Z,1));
lower_act = (smooth(mean(Z,1),3))' - std(Z,0,1)/sqrt(size(Z,1));
X_score = [x/10,fliplr(x/10)];
Y_score = [higher_act, fliplr(lower_act)];
d = fill(X_score,Y_score,[255 100 18]/255,'edgecolor','none');
alpha(d,0.5); hold on
plot(x/10,smooth(mean(Z,1),3),'color',[1 0 0],'linewidth',4); hold on
% suppression
Z = [];
unitdisplay = LaserFR(1,SuppressedUnitID);
for itr = 1:size(unitdisplay,2) % neuron
    Z(itr,:) = (mean(unitdisplay{itr})-mean(mean(unitdisplay{itr}(:,1:10),2)))/(0.00000000001+std(mean(unitdisplay{itr}(:,1:10),1)));
end
x = 1:size(Z,2); % time bin
higher_sup = (smooth(mean(Z,1),3))' + std(Z,0,1)/sqrt(size(Z,1));
lower_sup = (smooth(mean(Z,1),3))' - std(Z,0,1)/sqrt(size(Z,1));
X_score = [x/10,fliplr(x/10)];
Y_score = [higher_sup, fliplr(lower_sup)];
d = fill(X_score,Y_score,[0 0 0],'edgecolor','none');
MinNormFR = 1.1*min([min(lower_act) min(lower_sup)]);
MaxNormFR = 1.1*max([max(lower_act) max(lower_sup)]);
alpha(d,0.2);
hold on
plot(x/10,smooth(mean(Z,1),3),'color',[0 0 0],'linewidth',4);
hold on
rectangle('Position',[beforelaser,MaxNormFR-0.5,laserduration,0.5],'LineStyle','none','faceColor',[67 106 178]/255);
box('off');
set(gca,'XTick',[1 5],'XTickLabel',{'0','4'},'FontName','Arial','xlim',[0 beforelaser+laserduration+afterlaser+plus]);
xlabel('Time (sec)','FontSize',12,'FontName','Arial');
set(gca,'YTick',MinNormFR:5:MaxNormFR,'YTickLabel',{'-5','0','5','10','15'},'FontName','Arial','FontSize',16,'ylim',[MinNormFR MaxNormFR]);
hold on
plot([1 1],[MinNormFR MaxNormFR-0.5],'k--','linewidth',2);
hold on
plot([5 5],[MinNormFR MaxNormFR-0.5],'k--','linewidth',2);
ylabel('Normalized Firing Rate','FontSize',16,'FontName','Arial'); % can be changed at any time.
set(gcf,'Renderer','Painter'); saveas(gcf,['ZscoreFiringRate_' TargetBrain '_' Group],'fig'); close;

%% Spike raster and PSTH polts of target neurons
trialnum = 20;
for iUnit = 1:size(TarUnitsRG,2)
    NewTargetBrainUnitsRG{1,iUnit} = cellfun(@(x) x(x>=3 & x<=9)-3,TarUnitsRG{1,iUnit},'UniformOutput', 0);
    figure1 = figure;
    axe_upper = axes('parent',figure1,'position', [0.1 0.5 0.8 0.45],'XTickLabel',[],'YTickLabel',[]);
    hold(axe_upper,'all');
    SelectedTrialID = randperm(length(NewTargetBrainUnitsRG{1,iUnit}));
    for itr = 1:trialnum % trial
        for itr1 = 1:size(NewTargetBrainUnitsRG{1,iUnit}{SelectedTrialID(itr),1},1) % each RG
            plot([NewTargetBrainUnitsRG{1,iUnit}{SelectedTrialID(itr),1}(itr1,1) NewTargetBrainUnitsRG{1,iUnit}{SelectedTrialID(itr),1}(itr1,1)], [trialnum+1-itr-0.5 trialnum+1-itr+0.5],'k');
            hold on
        end
    end
    set(axe_upper,'XTick',zeros(1,0));
    xlim([0 6]);
    ylim([0.5 trialnum+0.5]);
    box off
    SingleUnitFR = LaserFR{1,iUnit};
    axe_down = axes('parent',figure1,'position', [0.1 0.05 0.8 0.4]);
    % hold(axe_down,'all');
    time = 1:size(SingleUnitFR,2);
    average = mean(SingleUnitFR,1);
    hold on
    plot(time/10, smooth(average,3), 'color','k','linewidth',5);
    hold on
    xlim([0 6]);
    box off
    saveas(gcf,[TargetBrain 'Unit' num2str(iUnit) 'FiringinLaser'],'fig');
    close;
end

