%% Plot accuracy of cross-temporal decoding

clear; clc; close all;

%% Assignment
Group = 'CtrlGroup';
Target = 'aAIC'; 
NeuronCodingType = 0; % 0: no remove; 1: remove sustained neurons; 2: remove transient neurons 3: only sustained 4: only transient 5: only non-memory neuron
if NeuronCodingType == 0
    RegionCoding = 'All';
elseif NeuronCodingType == 1
    RegionCoding = 'Remove sustained';
elseif NeuronCodingType == 2
    RegionCoding = 'Remove transient';
elseif NeuronCodingType == 3
    RegionCoding = 'Only sustained';
elseif NeuronCodingType == 4
    RegionCoding = 'Only transient';
elseif NeuronCodingType == 5
    RegionCoding = 'Only nonmemory';
end
X = 6:70;
ShuffleTimesForTest = 1000;

%% Load decoding result
switch(NeuronCodingType)
    case 0
        File_Recording = dir('*All_CTDecoding_Recording*.mat');
        File_Shuffled = dir('*All_CTDecoding_Shuffled*.mat');
    case 1
        File_Recording = dir('*Remove sustained_CTDecoding_Recording*.mat');
        File_Shuffled = dir('*Remove sustained_CTDecoding_Shuffled*.mat');
    case 2
        File_Recording = dir('*Remove transient_CTDecoding_Recording*.mat');
        File_Shuffled = dir('*Remove transient_CTDecoding_Shuffled*.mat');
    case 3
        File_Recording = dir('*Only sustained_CTDecoding_Recording*.mat');
        File_Shuffled = dir('*Only sustained_CTDecoding_Shuffled*.mat');
    case 4
        File_Recording = dir('*Only transient_CTDecoding_Recording*.mat');
        File_Shuffled = dir('*Only transient_CTDecoding_Shuffled*.mat');
    case 5
        File_Recording = dir('*Only nonmemory_CTDecoding_Recording*.mat');
        File_Shuffled = dir('*Only nonmemory_CTDecoding_Shuffled*.mat');
end
load(File_Recording.name);
RealDecodingResults = data;
load(File_Shuffled.name);
ShuffleDecodingResults = data;

%% Averaged result of cross-temporal decoding 
RealDecodingResults = mean(RealDecodingResults);
RealDecodingResults = reshape(RealDecodingResults,size(RealDecodingResults,2),size(RealDecodingResults,3));

%% Bin with significant decoding accuracy, based on cluster-based permutation test
TestBinNum = size(ShuffleDecodingResults,2);
[IsSig,ClusterBasedPermuTestIsSig,ShuffleIsSig,ShuffleSigClusterSize,RealDecodingResultsSigClusterSize,RealDecodingResultsSigClusterID,IsSigLessThanChance,ClusterBasedPermuTestIsSigLessChance,RealDecodingResultsSigLessClusterSize]...
    = ClusterBasedPermutationTest_CTD(RealDecodingResults,ShuffleDecodingResults,TestBinNum,ShuffleTimesForTest,2);

%% Plot CTD result and outline bins with significant decoding accuracy
figure1 = figure;
axes1 = axes('Parent',figure1,'position',[0.1 0.1 0.8 0.8],'FontSize',14,'FontName','Times New Roman');
box(axes1,'off');
hold(axes1,'all');
interval = 0;
colormap('jet');
RealDecodingResults = RealDecodingResults(X,X);
ClusterBasedPermuTestIsSig = ClusterBasedPermuTestIsSig(X,X);
ClusterBasedPermuTestIsSigLessChance = ClusterBasedPermuTestIsSigLessChance(X,X);
imagesc(RealDecodingResults,[0.35 0.85]);
hold on 
if exist('ClusterBasedPermuTestIsSignificant','var')
    contour(ClusterBasedPermuTestIsSig,[1 1],'-','color',[1 0 0],'linewidth',3)
    contour(ClusterBasedPermuTestIsSigLessChance,[1 1],'-','color',[1 0 0],'linewidth',3)
end
hold on
plot([5.5 5.5],[1 length(X)],'k','LineStyle','--','linewidth',3);
plot([15.5 15.5],[1 length(X)],'k','LineStyle','--','linewidth',3);
plot([55.5 55.5],[1 length(X)],'k','LineStyle','--','linewidth',3);
plot([60.5 60.5],[1 length(X)],'k','LineStyle','--','linewidth',3);
plot([1 length(X)],[5.5 5.5],'k','LineStyle','--','linewidth',3);
plot([1 length(X)],[15.5 15.5],'k','LineStyle','--','linewidth',3);
plot([1 length(X)],[55.5 55.5],'k','LineStyle','--','linewidth',3);
plot([1 length(X)],[60.5 60.5],'k','LineStyle','--','linewidth',3);
set(gca,'YDir','normal','XTick',[5.5 25.5 45.5 65.5],'XTickLabel',{'0','2','4','6'},'YTick',[5.5 25.5 45.5 65.5],'YTickLabel',{'0','2','4','6'},'FontName','Arial','FontSize',16);
xlabel('Test time ( s )','FontSize',18,'FontName','Arial');
ylabel('Train time ( s )','FontSize',18,'FontName','Arial'); 
location1 = strfind(File_Recording.name,'_n=');
location2 = strfind(File_Recording.name,'.mat');
set(gcf, 'Renderer', 'Painter'); saveas(gcf,[Target ' Neurons_' RegionCoding '_CTDecodingConfusionMatrix_' Group File_Recording.name(location1:location2-1)],'fig'); close;
save([Target ' Neurons_' RegionCoding '_ClusterBasedPermuTestIsSigBetterAndLowerthanShuffle_' Group File_Recording.name(location1:location2-1)],'ClusterBasedPermuTestIsSignificant','ClusterBasedPermuTestIsSigLessThanChance');
