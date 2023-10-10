%% Proportion of memory neurons over the course of learning days.

clear; clc; close all;

%% Assignment
Group = 'CtrlGroup';
TarReg = 'aAIC';

%% Load result
load(['NeuronOriginandSpikeRasterInformationfor' Group '.mat']); % neuron origin information
load(['Go-NoGo-PreferredUnitsID-' Group '.mat']); % Go- and NoGo-preferred neurons
load([TarReg 'SelectivityData_' Group '.mat']); % sustained and transient neurons
if strcmp(TarReg,'mPFC')
    TarRegUnitsOrigInfo = mPFCOriginInfo;
else
    TarRegUnitsOrigInfo = aAICOriginInfo;
end
TarRegUnitsOrigInfo = fliplr(TarRegUnitsOrigInfo);

%% ID of Go-preferred, NoGo-preferred, sustained, transient, and memory neurons
% Go-preferred and NoGo-preferred neurons
if strcmp(TarReg,'mPFC')
    GoPrefUnitsID = PreferCodingNeuronID.mPFC.GoPrefUnitsID; 
    NoGoPrefUnitsID = PreferCodingNeuronID.mPFC.NoGoPrefUnitsID; 
elseif strcmp(TarReg,'aAIC')
    GoPrefUnitsID = PreferCodingNeuronID.aAIC.GoPrefUnitsID;
    NoGoPrefUnitsID = PreferCodingNeuronID.aAIC.NoGoPrefUnitsID;
end
GoPrefUnitsID = GoPrefUnitsID(:);
NoGoPrefUnitsID = NoGoPrefUnitsID(:);
% sustained and transient neurons
TransientCodingNeuronID = [TransientCodingNeuronID; SwitchedSustainedUnitID];
SustainedUnitID = SustainedUnitID(:);
TransientCodingNeuronID = TransientCodingNeuronID(:);
% memory neurons
MemoryUnitID = vertcat(SustainedUnitID,TransientCodingNeuronID);
clearvars -except Group TarReg GoPrefUnitsID NoGoPrefUnitsID SustainedUnitID TransientCodingNeuronID MemoryUnitID TarRegUnitsOrigInfo

%% Origin information (learning day and mice) of Go-preferred, NoGo-preferred, sustained, transient, and memory neurons
GoPrefUnitsInfo = GetOriginInfo(GoPrefUnitsID,TarRegUnitsOrigInfo);
NoGoPrefUnitsInfo = GetOriginInfo(NoGoPrefUnitsID,TarRegUnitsOrigInfo);
SustUnitsInfo = GetOriginInfo(SustainedUnitID,TarRegUnitsOrigInfo);
TranUnitsInfo = GetOriginInfo(TransientCodingNeuronID,TarRegUnitsOrigInfo);
MemoryUnitsInfo = GetOriginInfo(MemoryUnitID,TarRegUnitsOrigInfo);

%% Proportion of memory neurons over the course of learning
UniqLearningDayID = unique(TarRegUnitsOrigInfo(:,1));
UniqMiceID = unique(TarRegUnitsOrigInfo(:,2));
Prop_GoPrefUnits = cell(1,numel(UniqLearningDayID));
Prop_NoGoPrefUnits = cell(1,numel(UniqLearningDayID));
Prop_SustUnits = cell(1,numel(UniqLearningDayID));
Prop_TranUnits = cell(1,numel(UniqLearningDayID));
Prop_MemoryUnits = cell(1,numel(UniqLearningDayID));
for iDay = 1:numel(UniqLearningDayID)
    for iMouse = 1:numel(UniqMiceID) 
        tempUnitsNum = nnz(TarRegUnitsOrigInfo(:,1)==UniqLearningDayID(iDay) & TarRegUnitsOrigInfo(:,2)==UniqMiceID(iMouse));
        % Go-preferred neurons
        tempGoPrefUnitsNum = nnz(GoPrefUnitsInfo(:,1)==UniqLearningDayID(iDay) & GoPrefUnitsInfo(:,2)==UniqMiceID(iMouse));
        Prop_GoPrefUnits{iDay} = [Prop_GoPrefUnits{iDay}; tempGoPrefUnitsNum/tempUnitsNum];
        % NoGo-preferred neurons
        tempNoGoPrefUnitsNum = nnz(NoGoPrefUnitsInfo(:,1)==UniqLearningDayID(iDay) & NoGoPrefUnitsInfo(:,2)==UniqMiceID(iMouse));
        Prop_NoGoPrefUnits{iDay} = [Prop_NoGoPrefUnits{iDay}; tempNoGoPrefUnitsNum/tempUnitsNum];
        % sustained neurons
        tempSustUnitsNum = nnz(SustUnitsInfo(:,1)==UniqLearningDayID(iDay) & SustUnitsInfo(:,2)==UniqMiceID(iMouse));
        Prop_SustUnits{iDay} = [Prop_SustUnits{iDay}; tempSustUnitsNum/tempUnitsNum];
        % transient neurons
        tempTranUnitsNum = nnz(TranUnitsInfo(:,1)==UniqLearningDayID(iDay) & TranUnitsInfo(:,2)==UniqMiceID(iMouse));
        Prop_TranUnits{iDay} = [Prop_TranUnits{iDay}; tempTranUnitsNum/tempUnitsNum];
        % memory neurons
        tempMemoryUnitsNum = nnz(MemoryUnitsInfo(:,1)==UniqLearningDayID(iDay) & MemoryUnitsInfo(:,2)==UniqMiceID(iMouse));
        Prop_MemoryUnits{iDay} = [Prop_MemoryUnits{iDay}; tempMemoryUnitsNum/tempUnitsNum];
    end
end

%% Figure (error bar plot)
% proportion of Go-preferred and NoGo-preferred neurons over learning
figure('position',[200 200 300 500]);
AverProp_GoPrefUnits = cellfun(@mean,Prop_GoPrefUnits);
SemProp_GoPrefUnits = cellfun(@(x) std(x)/sqrt(numel(x)),Prop_GoPrefUnits);
errorbar(1:numel(AverProp_GoPrefUnits),AverProp_GoPrefUnits,SemProp_GoPrefUnits,'r','marker','o','markerfacecolor','r','markeredgecolor','none'); hold on
AverProp_NoGoPrefUnits = cellfun(@mean,Prop_NoGoPrefUnits);
SemProp_NoGoPrefUnits = cellfun(@(x) std(x)/sqrt(numel(x)),Prop_NoGoPrefUnits);
errorbar(1:numel(AverProp_NoGoPrefUnits),AverProp_NoGoPrefUnits,SemProp_NoGoPrefUnits,'b','marker','o','markerfacecolor','b','markeredgecolor','none'); hold on
set(gca,'XTick',1:1:numel(UniqLearningDayID),'xlim',[0.6 numel(UniqLearningDayID)+0.4],'FontName','Arial','FontSize',16,'ylim',[0 1.25*max(max(AverProp_GoPrefUnits),max(AverProp_NoGoPrefUnits))]);
ylabel('Proportion of Go- and NoGo-preferred neurons (%)','FontSize',18,'FontName','Arial'); box off;
set(gcf,'Renderer','Painter'); saveas(gcf,sprintf('Proportion of %s Go- and NoGo-preferred neurons over learning_%s',TarReg,Group),'fig'); close all;
stat = TwRepeatedMeasureANOVA(Prop_GoPrefUnits,Prop_NoGoPrefUnits,{'day','type'});
save(sprintf('Proportion of %s Go- and NoGo-preferred neurons over learning_%s',TarReg,Group),'stat','-v7.3');
% proportion of sustained and transient neurons over learning
figure('position',[200 200 300 500]);
AverProp_SustUnits = cellfun(@mean,Prop_SustUnits);
SemProp_SustUnits = cellfun(@(x) std(x)/sqrt(numel(x)),Prop_SustUnits);
errorbar(1:numel(AverProp_SustUnits),AverProp_SustUnits,SemProp_SustUnits,'k','marker','o','markerfacecolor','k','markeredgecolor','none'); hold on
AverProp_TranUnits = cellfun(@mean,Prop_TranUnits);
SemProp_TranUnits = cellfun(@(x) std(x)/sqrt(numel(x)),Prop_TranUnits);
errorbar(1:numel(AverProp_TranUnits),AverProp_TranUnits,SemProp_TranUnits,'k--','marker','o','markerfacecolor','k','markeredgecolor','none'); hold on
set(gca,'XTick',1:1:numel(UniqLearningDayID),'xlim',[0.6 numel(UniqLearningDayID)+0.4],'FontName','Arial','FontSize',16,'ylim',[0 1.25*max(max(AverProp_SustUnits),max(AverProp_TranUnits))]);
ylabel('Proportion of sustained and transient neurons (%)','FontSize',18,'FontName','Arial'); box off;
set(gcf,'Renderer','Painter'); saveas(gcf,sprintf('Proportion of %s sustained and transient neurons over learning_%s',TarReg,Group),'fig'); close all;
stat = TwRepeatedMeasureANOVA(Prop_SustUnits,Prop_TranUnits,{'day','type'});
save(sprintf('Proportion of %s sustained and transient neurons over learning_%s',TarReg,Group),'stat','-v7.3');
% proportion of memory neurons over learning
figure('position',[200 200 300 500]);
AverProp_MemoryUnits = cellfun(@mean,Prop_MemoryUnits);
SemProp_MemoryUnits = cellfun(@(x) std(x)/sqrt(numel(x)),Prop_MemoryUnits);
errorbar(1:numel(AverProp_MemoryUnits),AverProp_MemoryUnits,SemProp_MemoryUnits,'k','marker','o','markerfacecolor','k','markeredgecolor','none'); hold on
set(gca,'XTick',1:1:numel(UniqLearningDayID),'xlim',[0.6 numel(UniqLearningDayID)+0.4],'FontName','Arial','FontSize',16,'ylim',[0 1.25*max(AverProp_MemoryUnits)]);
ylabel('Proportion of memory neurons (%)','FontSize',18,'FontName','Arial'); box off;
set(gcf,'Renderer','Painter'); saveas(gcf,sprintf('Proportion of %s memory neurons over learning_%s',TarReg,Group),'fig'); close all;




function stat = TwRepeatedMeasureANOVA(data_g1,data_g2,VarName)

%% Two-factor, within-subject repeated measures ANOVA.
% For designs with two within-subject factors.
DayNum = size(data_g1,2);
BetweenGroup = [];
MiceID = [];
BetweenDay = [];
for iMouse = 1:size(data_g1{1},1)
    MiceID = [MiceID; iMouse*ones(2*DayNum,1)];
    BetweenGroup = [BetweenGroup; vertcat(ones(DayNum,1),2*ones(DayNum,1))];
    BetweenDay = [BetweenDay; (horzcat(1:DayNum,1:DayNum))'];
end
data_g1 = horzcat(data_g1{:});
data_g2 = horzcat(data_g2{:});
Perf = (horzcat(data_g1,data_g2))';
Perf = Perf(:);
stat = rm_anova2(Perf,MiceID,BetweenDay,BetweenGroup,VarName);
end

function UnitsOriginInfo = GetOriginInfo(UnitsID,OriginInfo)

matrix = zeros(numel(UnitsID),2);
for iUnit = 1:size(UnitsID,1)
    tempUnitID = UnitsID(iUnit);
    matrix(iUnit,:) = OriginInfo(tempUnitID,1:2);
end
UnitsOriginInfo = horzcat(matrix,UnitsID);
end
