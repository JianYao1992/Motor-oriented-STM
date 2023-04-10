%% Trajectory and distance in PCA space

clear; clc; close all;

%% Assignment
TarReg = 'aAIC';
IsRawFRorNormFR = 1; % 1: raw FR; 2: norm FR
UnitsNumInAnalysis = 300;
TrialNuminEachCondition = 60;
ResamplingNum_Real = 100;
ResamplingNum_Shuffle = 1000;
C = [{[0 0 0]} {[0 1 0]} {[0 0 1]}];
% event duration
BaseLen = 2;
SampOdorLen = 1;
DelayLen = 4; 
TestOdorLen = 0.5; 
TimeGain =10;

%% Target directory
CurrentPath = uigetdir;
cd(CurrentPath);

%% Load result of FR
load(sprintf('%sFRinS1S2_CtrlGroup.mat',TarReg));

%% Trajectory
% input matrix
if IsRawFRorNormFR == 1
    FR_S1 = cellfun(@(x) x(:,1+BaseLen*TimeGain:(BaseLen+SampOdorLen+DelayLen)*TimeGain),TargetBrainUnitsFRinS1,'UniformOutput', 0);
    FR_S2 = cellfun(@(x) x(:,1+BaseLen*TimeGain:(BaseLen+SampOdorLen+DelayLen)*TimeGain),TargetBrainUnitsFRinS2,'UniformOutput', 0);
    MeanFR_S1 = cellfun(@mean,FR_S1,'UniformOutput',0);
    MeanFR_S1 = cellfun(@(x) smooth(x,5)',MeanFR_S1,'UniformOutput',0);
    MeanFR_S1 = (vertcat(MeanFR_S1{:}))';
    MeanFR_S2 = cellfun(@mean,FR_S2,'UniformOutput', 0);
    MeanFR_S2 = cellfun(@(x) smooth(x,5)',MeanFR_S2,'UniformOutput',0);
    MeanFR_S2 = (vertcat(MeanFR_S2{:}))';
elseif IsRawFRorNormFR == 2
    MeanFR_S1 = cellfun(@mean,TargetBrainUnitsFRinS1,'UniformOutput',0);
    MeanFR_S1 = cellfun(@(x) (x-mean(x(1:BaseLen*TimeGain)))/std(x(1:BaseLen*TimeGain)),MeanFR_S1,'UniformOutput',0); 
    MeanFR_S1 = cellfun(@(x) x(:,1+(BaseLen+SampOdorLen)*TimeGain:(BaseLen+SampOdorLen+DelayLen)*TimeGain),MeanFR_S1,'UniformOutput',0);
    MeanFR_S1 = (vertcat(MeanFR_S1{:}))';
    MeanFR_S2 = cellfun(@mean,TargetBrainUnitsFRinS2,'UniformOutput',0);
    MeanFR_S2 = cellfun(@(x) (x-mean(x(1:BaseLen*TimeGain)))/std(x(1:BaseLen*TimeGain)),MeanFR_S2,'UniformOutput',0);
    MeanFR_S2 = cellfun(@(x) x(:,1+(BaseLen+SampOdorLen)*TimeGain:(BaseLen+SampOdorLen+DelayLen)*TimeGain),MeanFR_S2,'UniformOutput',0);
    MeanFR_S2 = (vertcat(MeanFR_S2{:}))';
end
MeanFR = [MeanFR_S1; MeanFR_S2];
% eigenvector
[EigenVector,score,EigenValue,~,explained] = pca(MeanFR);
TotalBinNum = size(MeanFR,1);
Sample1Score = score(1:TotalBinNum/2,:);
Sample2Score = score(TotalBinNum/2+1:TotalBinNum,:);
% plot 3-D trajectory
figure;
% sample 1 trajectory
plot3(Sample1Score(:,1), Sample1Score(:,2),Sample1Score(:,3),'r','marker','o','markerfacecolor','r'); % delay period
hold on
% sample 2 trajectory
plot3(Sample2Score(:,1), Sample2Score(:,2),Sample2Score(:,3),'b','marker','o','markerfacecolor','b'); % delay period
hold on
xlabel('PC1','fontsize',14);
ylabel('PC2','fontsize',14);
zlabel('PC3','fontsize',14);
set(gca,'FontSize',14);
set(gcf,'Render','Painter'); saveas(gcf,[TarReg 'AllNeurons_3DPCATrajectory'],'fig'); close;

%% Distance between trajectories
tic
WorkerNumber = 10; NullDistributionNum = 10; num_resample_runs = 100;  % 100 cross validation
poolobj = gcp('nocreate'); % If no pool, do not create new one.
if isempty(poolobj)
    myCluster = parcluster('local'); myCluster.NumWorkers = WorkerNumber; parpool(myCluster,WorkerNumber);
end
% real data
disp('-----Step 1: analysis of trajectory distance for real data-----');
ResampDistance_Real = [];
for iNullDis=1:NullDistributionNum
    f(iNullDis) = parfeval(@TrajectoryDistanceCalculation,1,TargetBrainUnitsFRinS1,TargetBrainUnitsFRinS2,TrialNuminEachCondition,UnitsNumInAnalysis,ResamplingNum_Real/NullDistributionNum,0,BaseLen,SampOdorLen,DelayLen);
    disp(['---' num2str(iNullDis) '-Computing End---'])
end
for iNullDis = 1:NullDistributionNum
    [~,Distance] = fetchNext(f);
    ResampDistance_Real = [ResampDistance_Real; Distance];
end
% shuffled data
disp('-----Step 2: analysis of trajectory distance for shuffled data-----')
ResampDistance_Shuffle = [];
for iNullDis = 1:NullDistributionNum
    f_shuffle(iNullDis) = parfeval(@TrajectoryDistanceCalculation,1,TargetBrainUnitsFRinS1,TargetBrainUnitsFRinS2,TrialNuminEachCondition,UnitsNumInAnalysis,ResamplingNum_Shuffle/NullDistributionNum,1,BaseLen,SampOdorLen,DelayLen);
    disp(['---' num2str(iNullDis) '-Computing End---'])
end
for iNullDis=1:NullDistributionNum
    [~,ShuffleDistance] = fetchNext(f_shuffle);  % Collect the results as they become available.
    ResampDistance_Shuffle = [ResampDistance_Shuffle; ShuffleDistance];
end

%% Plot distance between trajectories
disp('-----Step 3: plot distance between trajectories for real and shuffled data-----');
% distance in real data
figure;
time = 1:size(ResampDistance_Real,2);
average = mean(ResampDistance_Real,1);
plot(time/TimeGain,smooth(mean(ResampDistance_Real,1),3),'r','linewidth',3);
hold on
temp95Percentile = prctile(ResampDistance_Real,[2.5 97.5],1);
Time = [time/TimeGain, fliplr(time/TimeGain)];
value = [smooth(temp95Percentile(2,:),3)', fliplr(smooth(temp95Percentile(1,:),3)')];
a = fill(Time,value,[1 0 0],'edgecolor','none');
alpha(a,0.2);
% distance in shuffled data
time = 1:size(ResampDistance_Shuffle,2);
average = mean(ResampDistance_Shuffle,1);
plot(time/TimeGain,smooth(mean(ResampDistance_Shuffle,1),3),'color',[0 0 0],'linestyle','--','linewidth',3);
hold on
temp95Percentile= prctile(ResampDistance_Shuffle,[2.5 97.5],1);
Time = [time/TimeGain, fliplr(time/TimeGain)];
value = [temp95Percentile(2,:), fliplr(temp95Percentile(1,:))];
a = fill(Time,value,[0 0 0],'edgecolor','none');
alpha(a,0.2);
ylim([5 max(mean(ResampDistance_Real,1))+10]);

%% Label ID of bin with significant distance between trajectories
[SigBinID,~] = ClusterBasedPermutationTest_YJ(ResampDistance_Real,ResampDistance_Shuffle,1,size(ResampDistance_Shuffle,1));
SigBinID = find(SigBinID==1);
for i = 1:length(SigBinID)
    patch([(SigBinID(i)-1)/TimeGain (SigBinID(i)-1)/TimeGain SigBinID(i)/TimeGain SigBinID(i)/TimeGain],[max(mean(ResampDistance_Real,1))+6 max(mean(ResampDistance_Real,1))+7 max(mean(ResampDistance_Real,1))+7 max(mean(ResampDistance_Real,1))+6],[0 0 0],'edgecolor','none');
    hold on
end
box off
plot([BaseLen BaseLen],[5 max(mean(ResampDistance_Real,1))+5],'color','k','linestyle','--','linewidth',1.5); %  sample odor
plot([BaseLen+SampOdorLen BaseLen+SampOdorLen],[5 max(mean(ResampDistance_Real,1))+5],'color','k','linestyle','--','linewidth',1.5);
plot([BaseLen+SampOdorLen+DelayLen BaseLen+SampOdorLen+DelayLen],[5 max(mean(ResampDistance_Real,1))+5],'color','k','linestyle','--','linewidth',1.5); %  response odor
plot([BaseLen+SampOdorLen+DelayLen+TestOdorLen BaseLen+SampOdorLen+DelayLen+TestOdorLen],[5 max(mean(ResampDistance_Real,1))+5],'color','k','linestyle','--','linewidth',1.5);
xlim([0.2 BaseLen+SampOdorLen+DelayLen+TestOdorLen+0.6])
set(gcf,'Render','Painter'); saveas(gcf,[TarReg 'AllNeurons_3DPCATrajectoryDistance_n=' num2str(UnitsNumInAnalysis)],'fig'); close;
toc
