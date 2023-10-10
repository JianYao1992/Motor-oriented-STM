%% Selectivity for S1 and S2 trials at population level

clear; clc; close all;

%% Assignment
TargetBrain = 'aAIC'; 
Group = 'CtrlGroup';
BaseLen = 9; 
ShownBaseLen = 2; 
SampOdorLen = 1; 
DelayLen = 4; 
TestOdorLen = 0.5; 
RespLen = 1;
AllUnitBinnedFR = [];
AllUnitBinnedRG = [];
AllUnitTrialID_S1 = [];
AllUnitTrialID_S2 = [];
ID_mPFC = []; 
ID_aAIC = [];
beforeneuronnum = 0;
TimeGain = 10;
Waveform = [];
C1 = [1 0 0]; 
C2 = [0 0 1]; 
C3 = [0 0 0]; 
C4 = [0 1 0]; 
BilateralaAICMiceID = [{'M0021'} {'M0022'} {'M0023'} {'M0024'}];
FileID = 1;

%% Target directory
CurrPath = uigetdir; 
AllPath = genpath(CurrPath);
SplitPath = strsplit(AllPath,';');
SubPath = SplitPath';
SubPath = SubPath(2:end-1);

%% ID and FR of trials for all neurons
for iPath = 1:size(SubPath,1)
    Path = SubPath{iPath,1};
    cd(Path);
    MatFiles = dir('*short*.mat');
    for j = 1:size(MatFiles,1)
        Filename{1,FileID} = MatFiles(j,1).name;
        load(Filename{1,FileID});
        if ~isempty(SingleUnitList)
            for itru = 1:size(SingleUnitList,1) % neuron
                AllUnitBinnedRG = [AllUnitBinnedRG {(RGResults(itru,:))'}];
                tempsingleunitBinnedFR = [];
                for itrt = 1:size(TrialMark,1) % trial
                    tempsingleunitBinnedFR = [tempsingleunitBinnedFR; Results{1,itrt}(itru,:)];
                end
                AllUnitBinnedFR = [AllUnitBinnedFR {tempsingleunitBinnedFR}];
            end
            TrialID_S1 = repmat({find(TrialMark(:,3)==1)},1,size(SingleUnitList,1));
            TrialID_S2 = repmat({find(TrialMark(:,3)==2)},1,size(SingleUnitList,1));
            AllUnitTrialID_S1 = [AllUnitTrialID_S1 TrialID_S1];
            AllUnitTrialID_S2 = [AllUnitTrialID_S2 TrialID_S2];
            % ID of mPFC neurons
            tempID_mPFC = [];
            Position = [];
            for k = 1:length(BilateralaAICMiceID)
                Position = [Position regexp(Filename{1,FileID},BilateralaAICMiceID{1,k})];
            end
            if ~isempty(Position)
                tempID_mPFC = []; 
            else
                tempID_mPFC = find(SingleUnitList(:,1)>=1 & SingleUnitList(:,1)<=16);
            end
            if ~isempty(tempID_mPFC)
                ID_mPFC = [ID_mPFC tempID_mPFC'+beforeneuronnum];
            end
            % ID of aAIC neurons
            tempID_aAIC = [];
            if ~isempty(Position)
                tempID_aAIC = find(SingleUnitList(:,1)>=1 & SingleUnitList(:,1)<=16);
            else
                tempID_aAIC = find(SingleUnitList(:,1)>=17 & SingleUnitList(:,1)<=32);
            end
            if ~isempty(tempID_aAIC)
                ID_aAIC = [ID_aAIC tempID_aAIC'+beforeneuronnum];
            end
            beforeneuronnum = beforeneuronnum + size(SingleUnitList,1);
            % waveform
            Waveform = [Waveform; NewWaveForm];
        end
        FileID = FileID + 1;
    end
end
cd(CurrPath);
AllUnitTrialID = [AllUnitTrialID_S1; AllUnitTrialID_S2];

%% Binned FR of target neurons
windowsize = 0.1; sliding = 0.1; % 100ms window, sliding 100ms
if strcmp(TargetBrain,'mPFC')
    TarUnitsID = ID_mPFC;
else
    TarUnitsID = ID_aAIC;
end
TarUnitBinnedFR = AllUnitBinnedFR(1,TarUnitsID);
TarRGResults = AllUnitBinnedRG(1,TarUnitsID);
TarUnitsTrialID = AllUnitTrialID(:,TarUnitsID);
NewTarUnitBinnedFR = cell(1,size(TarUnitBinnedFR,2));
[TargetBrainRGinS1,~,~,TargetBrainRGinS2] = ExtractSpecTrialNumForEachNeuron(TarRGResults,TarUnitsTrialID,[],0,0);
for i = 1:size(TarUnitBinnedFR,2) % neuron
    for j = 1:size(TarUnitBinnedFR{1,i},1) % trial
        for k = 1:sliding*TimeGain:size(TarUnitBinnedFR{1,i},2)-(windowsize*TimeGain-1) % time bin
            NewTarUnitBinnedFR{1,i}(j,k) = sum(TarUnitBinnedFR{1,i}(j,k:k+windowsize*TimeGain-1))/windowsize;
        end
    end
end

% %% Plot FR of target neurons
% TarPopuFR = cellfun(@mean,NewTarUnitBinnedFR,'UniformOutput', 0);
% TarPopuFR = vertcat(TarPopuFR{:});
% TarPopuFR = TarPopuFR(:,31:100);
% figure;
% plotshadow(TarPopuFR,[0 0 0],2,3,0);
% hold on
% time = 1:size(TarPopuFR,2);
% plot(time/10, smooth(mean(TarPopuFR),3), 'color',[0 0 0],'linewidth',2);
% plot([1 1],[0 6.3],'color','k','linestyle','--','linewidth',1.5); %  sample odor
% plot([2 2],[0 6.3],'color','k','linestyle','--','linewidth',1.5);
% plot([6 6],[0 6.3],'color','k','linestyle','--','linewidth',1.5); %  response odor
% plot([6.5 6.5],[0 6.3],'color','k','linestyle','--','linewidth',1.5);
% set(gca,'XTick',1:1:7,'XTickLabel',{'0','1','2','3','4','5','6'},'xlim',[0.1 7.4]);
% set(gca,'YTick',5:0.5:6.5,'YTickLabel',{'5','5.5','6','6.5'},'FontName','Arial','FontSize',16,'ylim',[4 7]);
% ylabel('Population Averaged Firing Rate ( Hz )','FontSize',18,'FontName','Arial');
% box off;
% set(gcf,'Renderer','Painter'); saveas(gcf,[TargetBrain '-population averaged FR-' Group],'fig'); close;

%% Plot real and shuffled selectivity
[RealSelectivity,ShuffleSelectivity] = GetRealShuffleData(NewTarUnitBinnedFR,TarUnitsTrialID,[]);
figure;
% real selectivity
time = 1:size(RealSelectivity,2);
Realaverage = mean(RealSelectivity,1);
Highervalue = (smooth(Realaverage,3))' + std(RealSelectivity,0,1)/(0.0001+sqrt(size(RealSelectivity,1)));
Lowervalue = (smooth(Realaverage,3))' - std(RealSelectivity,0,1)/(0.0001+sqrt(size(RealSelectivity,1)));
Time = [time/10, fliplr(time/10)];
value = [Highervalue, fliplr(Lowervalue)];
a = fill(Time, value, [0 0 0], 'edgecolor','none');
alpha(a,0.2);
hold on
plot(time/10, smooth(Realaverage,3), 'color',[0 0 0],'linewidth',2);
hold on
% shuffle selectivity
time = 1:size(ShuffleSelectivity,2);
Shuffleaverage = mean(ShuffleSelectivity,1);
temp95Percentile= prctile(ShuffleSelectivity,[2.5 97.5],1);  % 2.5th and 97.5th
Highervalue = temp95Percentile(2,:);
Lowervalue = temp95Percentile(1,:);
Time = [time/10, fliplr(time/10)];
value = [Highervalue, fliplr(Lowervalue)];
a = fill(Time, value, [0 0 0], 'edgecolor','none');
alpha(a,0.2);
hold on
plot(time/10, smooth(Shuffleaverage,3), 'color',[0 0 0],'linestyle','--','linewidth',2);
hold on
% cluster-based permutation test
[SigTime,~] = ClusterPermutationTest_ForRealAndShuffle(RealSelectivity, ShuffleSelectivity);
plot([BaseLen+1/(2*TimeGain) BaseLen+1/(2*TimeGain)],[0 0.22],'color','k','linestyle','--','linewidth',1.5); %  sample odor
plot([BaseLen+SampOdorLen+1/(2*TimeGain) BaseLen+SampOdorLen+1/(2*TimeGain)],[0 0.22],'color','k','linestyle','--','linewidth',1.5);
plot([BaseLen+SampOdorLen+DelayLen+1/(2*TimeGain) BaseLen+SampOdorLen+DelayLen+1/(2*TimeGain)],[0 0.22],'color','k','linestyle','--','linewidth',1.5); %  response odor
plot([BaseLen+SampOdorLen+DelayLen+TestOdorLen+1/(2*TimeGain) BaseLen+SampOdorLen+DelayLen+TestOdorLen+1/(2*TimeGain)],[0 0.22],'color','k','linestyle','--','linewidth',1.5);
hold on
for i=1:length(SigTime)
    patch([SigTime(i)-0.5 SigTime(i)-0.5 SigTime(i)+0.5 SigTime(i)+0.5]./TimeGain,[0.235 0.24 0.24 0.235],C3,'edgecolor','none');
    hold on
end
set(gca,'XTick',BaseLen+1/(2*TimeGain):1:BaseLen+SampOdorLen+DelayLen+TestOdorLen+1/(2*TimeGain),'XTickLabel',num2cell(floor(BaseLen+1/(2*TimeGain):1:BaseLen+SampOdorLen+DelayLen+TestOdorLen+1/(2*TimeGain))),'xlim',[BaseLen-ShownBaseLen+1/(2*TimeGain) BaseLen+SampOdorLen+DelayLen+TestOdorLen+RespLen+1/(2*TimeGain)]);
set(gca,'YTick',0.1:0.1:0.4,'YTickLabel',{'0.1','0.2','0.3','0.4'},'FontName','Arial','FontSize',16,'ylim',[0.09 0.4]);
box off;
set(gcf,'Renderer','Painter'); saveas(gcf,[TargetBrain '-SelectivityPSTH-' Group],'fig'); close;

%% ID of bins with significant selectivity for each neuron
[TargetBrainUnitsFRinS1,~,~,TargetBrainUnitsFRinS2] = ExtractSpecTrialNumForEachNeuron(NewTarUnitBinnedFR,TarUnitsTrialID,[],0,0);
[~,TargetBrainUnitsSigBinID] = GetPutativeSigUnitID(TargetBrainUnitsFRinS1,TargetBrainUnitsFRinS2,TimeGain,BaseLen,SampOdorLen,DelayLen,TestOdorLen,ShownBaseLen,RespLen/2,0);

%% Heatmap of selectivity for target neurons ///not plot for significant bins of all neurons, sorted by FR(S1)-FR(S2) during delay///
SortingPeriod = 1+TimeGain*(BaseLen+SampOdorLen):TimeGain*(BaseLen+SampOdorLen+DelayLen);
IsVisible = 'off';
SortedID = SortingByDiffBetweenTwoSampleTrials(TargetBrainUnitsFRinS1,TargetBrainUnitsFRinS2,SortingPeriod,1,BaseLen,SampOdorLen,DelayLen,TestOdorLen,TimeGain);
set(gcf,'Renderer','Painter'); saveas(gcf,[TargetBrain 'SelectivityHeatMap_' Group],'fig'); close;

%% Heatmap of normalzied FR in S1 and S2 trials ///same sorting order as heatmap of selectivity///
TargetBrainNormFRinS1 = [];
TargetBrainNormFRinS2 = [];
for itr = 1:size(TargetBrainUnitsFRinS1,2) % neuron
    TargetBrainNormFRinS1 = [TargetBrainNormFRinS1; (mean(TargetBrainUnitsFRinS1{1,itr})-mean(mean(TargetBrainUnitsFRinS1{1,itr}(:,1+TimeGain*(BaseLen-ShownBaseLen):TimeGain*BaseLen))))/std(mean(TargetBrainUnitsFRinS1{1,itr}(:,1+TimeGain*(BaseLen-ShownBaseLen):TimeGain*BaseLen),1))]; % compress cross-trial into one, and std for baseline time but not for cross-trial one!
    TargetBrainNormFRinS2 = [TargetBrainNormFRinS2; (mean(TargetBrainUnitsFRinS2{1,itr})-mean(mean(TargetBrainUnitsFRinS2{1,itr}(:,1+TimeGain*(BaseLen-ShownBaseLen):TimeGain*BaseLen))))/std(mean(TargetBrainUnitsFRinS2{1,itr}(:,1+TimeGain*(BaseLen-ShownBaseLen):TimeGain*BaseLen),1))];
end
for i = 1:length(SortedID)
    SortedTargetBrainNormFRinS1(i,:) = TargetBrainNormFRinS1(SortedID(i),:);
    SortedTargetBrainNormFRinS2(i,:) = TargetBrainNormFRinS2(SortedID(i),:);
end
SmoothSortedTargetBrainNormFRinS1 = [];
SmoothSortedTargetBrainNormFRinS2 = [];
for i = 1:size(SortedTargetBrainNormFRinS1,1)
    SmoothSortedTargetBrainNormFRinS1 = [SmoothSortedTargetBrainNormFRinS1; smooth(SortedTargetBrainNormFRinS1(i,:),3)'];
    SmoothSortedTargetBrainNormFRinS2 = [SmoothSortedTargetBrainNormFRinS2; smooth(SortedTargetBrainNormFRinS2(i,:),3)'];
end
% S1 trials
figure;
colormap('jet');
imagesc(SmoothSortedTargetBrainNormFRinS1,[-5 5]);
box('off');
hold on
unitnum = size(SmoothSortedTargetBrainNormFRinS1,1);
plot([TimeGain*(BaseLen+1/(2*TimeGain)) TimeGain*(BaseLen+1/(2*TimeGain))],[0.5 unitnum + 0.5],'color','k','linestyle','--','linewidth',3); % DPA sample odor
plot([TimeGain*(BaseLen+SampOdorLen+1/(2*TimeGain)) TimeGain*(BaseLen+SampOdorLen+1/(2*TimeGain))],[0.5 unitnum + 0.5],'color','k','linestyle','--','linewidth',3);
plot([TimeGain*(BaseLen+SampOdorLen+DelayLen+1/(2*TimeGain)) TimeGain*(BaseLen+SampOdorLen+DelayLen+1/(2*TimeGain))],[0.5 unitnum + 0.5],'color','k','linestyle','--','linewidth',3); % DPA test odor
plot([TimeGain*(BaseLen+SampOdorLen+DelayLen+TestOdorLen+1/(2*TimeGain)) TimeGain*(BaseLen+SampOdorLen+DelayLen+TestOdorLen+1/(2*TimeGain))],[0.5 unitnum + 0.5],'color','k','linestyle','--','linewidth',3);
set(gca,'XTick',[TimeGain*(BaseLen+1/(2*TimeGain)) TimeGain*(BaseLen+SampOdorLen+1/(2*TimeGain)) TimeGain*(BaseLen+SampOdorLen+DelayLen+1/(2*TimeGain)) TimeGain*(BaseLen+SampOdorLen+DelayLen+TestOdorLen+1/(2*TimeGain))],'XTickLabel',{'0','1','5','5.5'},'FontName','Arial','FontSize',16);
set(gca,'XTickLabel',num2cell([0 SampOdorLen SampOdorLen+DelayLen SampOdorLen+DelayLen+TestOdorLen]),'FontName','Arial','FontSize',16);
xlabel('Time (sec)','FontSize',18,'FontName','Arial');
box off;
cd(CurrPath);
set(gcf,'Renderer','Painter'); saveas(gcf,[TargetBrain 'NormFRinS1-' Group],'fig'); close all;
% S2 trials
figure;
colormap('jet');
imagesc(SmoothSortedTargetBrainNormFRinS2,[-5 5]);
box('off');
hold on
unitnum = size(SmoothSortedTargetBrainNormFRinS2,1);
plot([TimeGain*(BaseLen+1/(2*TimeGain)) TimeGain*(BaseLen+1/(2*TimeGain))],[0.5 unitnum + 0.5],'color','k','linestyle','--','linewidth',3); % DPA sample odor
plot([TimeGain*(BaseLen+SampOdorLen+1/(2*TimeGain)) TimeGain*(BaseLen+SampOdorLen+1/(2*TimeGain))],[0.5 unitnum + 0.5],'color','k','linestyle','--','linewidth',3);
plot([TimeGain*(BaseLen+SampOdorLen+DelayLen+1/(2*TimeGain)) TimeGain*(BaseLen+SampOdorLen+DelayLen+1/(2*TimeGain))],[0.5 unitnum + 0.5],'color','k','linestyle','--','linewidth',3); % DPA test odor
plot([TimeGain*(BaseLen+SampOdorLen+DelayLen+TestOdorLen+1/(2*TimeGain)) TimeGain*(BaseLen+SampOdorLen+DelayLen+TestOdorLen+1/(2*TimeGain))],[0.5 unitnum + 0.5],'color','k','linestyle','--','linewidth',3);
set(gca,'XTick',[TimeGain*(BaseLen+1/(2*TimeGain)) TimeGain*(BaseLen+SampOdorLen+1/(2*TimeGain)) TimeGain*(BaseLen+SampOdorLen+DelayLen+1/(2*TimeGain)) TimeGain*(BaseLen+SampOdorLen+DelayLen+TestOdorLen+1/(2*TimeGain))],'XTickLabel',{'0','1','5','5.5'},'FontName','Arial','FontSize',16);
set(gca,'XTickLabel',num2cell([0 SampOdorLen SampOdorLen+DelayLen SampOdorLen+DelayLen+TestOdorLen]),'FontName','Arial','FontSize',16);
set(gca,'XTickLabel',{'0','1','5','5.5'},'FontName','Arial','FontSize',16);
xlabel('Time (sec)','FontSize',18,'FontName','Arial');
box off;
cd(CurrPath);
set(gcf,'Renderer','Painter'); saveas(gcf,[TargetBrain 'NormFRinS2-' Group],'fig'); close all;

%% Matrix of 'significant' selectivity for target neurons
TarSigSelecBinNum = zeros(length(TarUnitsID),1);
p = []; 
TargetBrainSigSelectivity = [];
FDR = [];
for i = 1:size(TargetBrainUnitsFRinS1,2) % neuron
    for j = 1:floor(size(TargetBrainUnitsFRinS1{1,1},2)/TimeGain) % window
        p(i,j) = ranksum(mean(TargetBrainUnitsFRinS1{1,SortedID(i)}(:,2+TimeGain*(j-1):TimeGain*j+1),2),mean(TargetBrainUnitsFRinS2{1,SortedID(i)}(:,2+TimeGain*(j-1):TimeGain*j+1),2));
    end
    FDR(i,:) = (ShownBaseLen+SampOdorLen+DelayLen+TestOdorLen+RespLen/2)*p(i,:);
    for j = 1:floor(size(TargetBrainUnitsFRinS1{1,1},2)/TimeGain) % window
        if FDR(i,j) <= 0.05
            TargetBrainSigSelectivity(i,1+TimeGain*(j-1):TimeGain*j) = (mean(mean(TargetBrainUnitsFRinS1{1,SortedID(i)}(:,2+TimeGain*(j-1):TimeGain*j+1),2))-mean(mean(TargetBrainUnitsFRinS2{1,SortedID(i)}(:,2+TimeGain*(j-1):TimeGain*j+1),2)))/(mean(mean(TargetBrainUnitsFRinS1{1,SortedID(i)}(:,2+TimeGain*(j-1):TimeGain*j+1),2))+mean(mean(TargetBrainUnitsFRinS2{1,SortedID(i)}(:,2+TimeGain*(j-1):TimeGain*j+1),2)));
            if j >= BaseLen+SampOdorLen+1 && j <= BaseLen+SampOdorLen+DelayLen
                TarSigSelecBinNum(i,1) = TarSigSelecBinNum(i,1) + 1;
            end
        else
            TargetBrainSigSelectivity(i,1+TimeGain*(j-1):TimeGain*j) = 0;
        end
    end
end

%% ID of sustained neurons
PutativeSustainedUnitID = SortedID(find(TarSigSelecBinNum==4),1);
% exclude neurons with switched selectivity
SustainedUnitID = []; 
SwitchedSustainedUnitID = [];
for i=1:length(PutativeSustainedUnitID)
    SelectivityValue = [];
    for j=1:(length(SortingPeriod)/TimeGain) % delay window
        if mean(TargetBrainSigSelectivity(find(SortedID==PutativeSustainedUnitID(i)),(1+TimeGain*(BaseLen+SampOdorLen+j-1)):(TimeGain*(BaseLen+SampOdorLen+j))),2) > 0
            SelectivityValue(1,j) = 1;
        else
            SelectivityValue(1,j) = -1;
        end
    end
    if isempty(setdiff(SelectivityValue,ones(1,(length(SortingPeriod)/TimeGain)))) || isempty(setdiff(SelectivityValue,-1*ones(1,(length(SortingPeriod)/TimeGain))))
        SustainedUnitID = [SustainedUnitID; PutativeSustainedUnitID(i)];
    else
        SwitchedSustainedUnitID = [SwitchedSustainedUnitID; PutativeSustainedUnitID(i)];
    end
end

%% ID of transient neurons
PutativeTransientUnitID = SortedID(find(TarSigSelecBinNum >= 1 & TarSigSelecBinNum < DelayLen),1);
NoTransientCodingNeuronID =[];
% Z-stat value and corresponding p value for target neurons
ZStatistics = cell(1,length(TarUnitsID));
P = cell(1,length(TarUnitsID));
for i=1:length(PutativeTransientUnitID)
    [tempNoTransientCodingNeuronID,tempSigBinID,tempZStatistics,tempP] = TransientCodingTest(PutativeTransientUnitID(i),TargetBrainUnitsFRinS1{PutativeTransientUnitID(i)},TargetBrainUnitsFRinS2{PutativeTransientUnitID(i)},TargetBrainUnitsSigBinID{PutativeTransientUnitID(i)},BaseLen,SampOdorLen,DelayLen,TimeGain);
    NoTransientCodingNeuronID = [NoTransientCodingNeuronID; tempNoTransientCodingNeuronID];
    TargetBrainUnitsSigBinID{PutativeTransientUnitID(i)} = tempSigBinID;
    ZStatistics{PutativeTransientUnitID(i)} = tempZStatistics;
    P{PutativeTransientUnitID(i)} = tempP;
    disp(strcat('///Finish testing ',num2str(i),'th transient coding unit///'));
end
% change selectivity value of NoTransientCodingNeuronID to 0
for i=1:length(NoTransientCodingNeuronID)
    TargetBrainSigSelectivity(find(SortedID==NoTransientCodingNeuronID(i)),1+TimeGain*(BaseLen+SampOdorLen):TimeGain*(BaseLen+SampOdorLen+DelayLen)) = zeros(1,TimeGain*DelayLen);
end
% change selectivity value of new non-significant bins to 0 for TransientCodingNeuronID
TransientCodingNeuronID = setdiff(PutativeTransientUnitID,NoTransientCodingNeuronID);
for i=1:length(TransientCodingNeuronID)
    tempRowID = find(SortedID==TransientCodingNeuronID(i));
    tempNoSigBinID = setdiff(1:floor(size(TargetBrainUnitsFRinS1{1,1},2)/TimeGain),TargetBrainUnitsSigBinID{TransientCodingNeuronID(i)});
    for j = 1:length(tempNoSigBinID)
        TargetBrainSigSelectivity(tempRowID,1+TimeGain*(tempNoSigBinID(j)-1):TimeGain*tempNoSigBinID(j)) = zeros(1,TimeGain);
    end
end

%% PSTH plots of target neurons
for i=1:length(TarUnitsID) % neurons
    figure;
    plot((1:size(TargetBrainUnitsFRinS1{1,i},2))/TimeGain, smooth(mean(TargetBrainUnitsFRinS1{1,i},1),3)', 'color',C1,'linewidth',2);
    hold on
    plot((1:size(TargetBrainUnitsFRinS2{1,i},2))/TimeGain, smooth(mean(TargetBrainUnitsFRinS2{1,i},1),3)', 'color',C2,'linewidth',2);
    hold on
    patch([BaseLen BaseLen BaseLen+SampOdorLen BaseLen+SampOdorLen],[0 100 100 0],'k','FaceAlpha',0.2,'edgecolor','none');
    hold on
    patch([BaseLen+SampOdorLen+DelayLen BaseLen+SampOdorLen+DelayLen BaseLen+SampOdorLen+DelayLen+TestOdorLen BaseLen+SampOdorLen+DelayLen+TestOdorLen],[0 100 100 0],'k','FaceAlpha',0.2,'edgecolor','none');
    hold on
    plotshadow(TargetBrainUnitsFRinS1{1,i},C1,2,3,0);
    hold on
    plotshadow(TargetBrainUnitsFRinS2{1,i},C2,2,3,0);
    tempmaxFRinS1Trials = max(mean(TargetBrainUnitsFRinS1{1,i},1));
    tempmaxFRinS2Trials = max(mean(TargetBrainUnitsFRinS2{1,i},1));
    tempmaxFR = max([tempmaxFRinS1Trials tempmaxFRinS2Trials])+5;
    if ~isempty(TargetBrainUnitsSigBinID{i})
        hold on
        for j=1:length(TargetBrainUnitsSigBinID{i})
            patch([TargetBrainUnitsSigBinID{i}(j)-1 TargetBrainUnitsSigBinID{i}(j)-1 TargetBrainUnitsSigBinID{i}(j) TargetBrainUnitsSigBinID{i}(j)],[tempmaxFR-2.5 tempmaxFR-2 tempmaxFR-2 tempmaxFR-2.5],C3,'edgecolor','none');
            hold on
        end
    end
    set(gca,'XTick',1:1:floor(size(TargetBrainUnitsFRinS1{1,1},2)/TimeGain),'FontName','Arial','FontSize',16,'Xlim',[BaseLen-0.5 floor(size(TargetBrainUnitsFRinS1{1,1},2)/TimeGain)]);
    ylim([0 tempmaxFR]);
    xlabel('Time (sec)','FontSize',18,'FontName','Arial');
    ylabel('FR (Hz)','FontSize',18,'FontName','Arial');
    box off;
    set(gcf,'Renderer','Painter'); saveas(gcf,[TargetBrain '#' num2str(i) ],'fig'); close all;
end

%% Neurons with 1-, 2-, and 3-second STM-encoding ability during delay.
UnitIDWith1sCoding = [];
UnitIDWith2sCoding = []; 
UnitIDWith3sCoding = [];
for iUnit=1:length(TransientCodingNeuronID)
    SigNum = length(TargetBrainUnitsSigBinID{TransientCodingNeuronID(iUnit)}(TargetBrainUnitsSigBinID{TransientCodingNeuronID(iUnit)}>=BaseLen+SampOdorLen+1 & TargetBrainUnitsSigBinID{TransientCodingNeuronID(iUnit)}<=BaseLen+SampOdorLen+DelayLen));
    if SigNum == 1
        UnitIDWith1sCoding = [UnitIDWith1sCoding TransientCodingNeuronID(iUnit)];
    elseif SigNum == 2
        UnitIDWith2sCoding = [UnitIDWith2sCoding TransientCodingNeuronID(iUnit)];
    elseif SigNum == 3
        UnitIDWith3sCoding = [UnitIDWith3sCoding TransientCodingNeuronID(iUnit)];
    end
end
RatioUnitWith1sCoding = length(UnitIDWith1sCoding)/length(TarUnitsID);
RatioUnitWith2sCoding = length(UnitIDWith2sCoding)/length(TarUnitsID);
RatioUnitWith3sCoding = length(UnitIDWith3sCoding)/length(TarUnitsID);
RatioUnitWith4sCoding = (length(SustainedUnitID)+length(SwitchedSustainedUnitID))/length(TarUnitsID);

%% Proportion of neurons with STM-encoding ability over the time course of delay.
DelayBinCodingUnitID = cell(1,length(SortingPeriod)/TimeGain);
for iUnit=1:length(TarUnitsID)
    if ~isempty(find(TargetBrainUnitsSigBinID{iUnit}==BaseLen+SampOdorLen+1)) % 1st bin
        DelayBinCodingUnitID{1} = [DelayBinCodingUnitID{1} iUnit];
    end
    if ~isempty(find(TargetBrainUnitsSigBinID{iUnit}==BaseLen+SampOdorLen+2)) % 2nd bin
        DelayBinCodingUnitID{2} = [DelayBinCodingUnitID{2} iUnit];
    end
    if ~isempty(find(TargetBrainUnitsSigBinID{iUnit}==BaseLen+SampOdorLen+3)) % 3rd bin
        DelayBinCodingUnitID{3} = [DelayBinCodingUnitID{3} iUnit];
    end
    if ~isempty(find(TargetBrainUnitsSigBinID{iUnit}==BaseLen+SampOdorLen+4)) % 4th bin
        DelayBinCodingUnitID{4} = [DelayBinCodingUnitID{4} iUnit];
    end
end

%% ID of switched neurons with  2- and 3-second STM-encoding ability
% 2-second STM-encoding ability
Switched2sCodingUnitID = [];
for i=1:length(UnitIDWith2sCoding)
    SelectivityValue = [];
    for j=1:(length(SortingPeriod)/TimeGain) % window of delay
        if mean(TargetBrainSigSelectivity(find(SortedID==UnitIDWith2sCoding(i)),(1+TimeGain*(BaseLen+SampOdorLen+j-1)):(TimeGain*(BaseLen+SampOdorLen+j))),2) > 0
            SelectivityValue(1,j) = 1;
        elseif mean(TargetBrainSigSelectivity(find(SortedID==UnitIDWith2sCoding(i)),(1+TimeGain*(BaseLen+SampOdorLen+j-1)):(TimeGain*(BaseLen+SampOdorLen+j))),2) < 0
            SelectivityValue(1,j) = -1;
        else
            SelectivityValue(1,j) = 0;
        end
    end
    if sum(SelectivityValue) == 0
        Switched2sCodingUnitID = [Switched2sCodingUnitID; UnitIDWith2sCoding(i)];
    end
end
UnitIDWith2sConsistentCoding = setdiff(UnitIDWith2sCoding,Switched2sCodingUnitID);
% 3-second STM-encoding ability
UnitIDWith3sConsistentCoding = [];
for i=1:length(UnitIDWith3sCoding)
    SelectivityValue = [];
    for j=1:(length(SortingPeriod)/TimeGain) % window of delay
        if mean(TargetBrainSigSelectivity(find(SortedID==UnitIDWith3sCoding(i)),(1+TimeGain*(BaseLen+SampOdorLen+j-1)):(TimeGain*(BaseLen+SampOdorLen+j))),2) > 0
            SelectivityValue(1,j) = 1;
        elseif mean(TargetBrainSigSelectivity(find(SortedID==UnitIDWith3sCoding(i)),(1+TimeGain*(BaseLen+SampOdorLen+j-1)):(TimeGain*(BaseLen+SampOdorLen+j))),2) < 0
            SelectivityValue(1,j) = -1;
        else
            SelectivityValue(1,j) = 0;
        end
    end
    if sum(SelectivityValue) == 3 || sum(SelectivityValue) == -3
        UnitIDWith3sConsistentCoding = [UnitIDWith3sConsistentCoding; UnitIDWith3sCoding(i)];
    end
end
Switched3sCodingUnitID = setdiff(UnitIDWith3sCoding,UnitIDWith3sConsistentCoding);
save([TargetBrain 'SelectivityData_' Group],'SortedID','TargetBrainUnitsSigBinID','TargetBrainSigSelectivity','DelayBinCodingUnitID','SustainedUnitID',...
    'SwitchedSustainedUnitID','TransientCodingNeuronID','NoTransientCodingNeuronID','P','ZStatistics','UnitIDWith1sCoding','UnitIDWith2sCoding','UnitIDWith2sConsistentCoding',...
    'Switched2sCodingUnitID','UnitIDWith3sCoding','UnitIDWith3sConsistentCoding','Switched3sCodingUnitID','TargetBrainRGinS1','TargetBrainRGinS2','TargetBrainUnitsFRinS1','TargetBrainUnitsFRinS2','TarUnitsID');

% %% Heatmap of significant selectivity for memory neurons
% PlotHeatMapOfSigSelecForMemoryNeurons([TargetBrain 'SelectivityData_' Group '.mat']);