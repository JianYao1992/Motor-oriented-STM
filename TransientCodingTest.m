%% Judge whether individual neuron is transient neuron during delay.
% input data is FR of putative transient neurons

function [NoTransientCodingNeuronID,tempSigBinID,tempZStatistics,tempP] = TransientCodingTest(NeuronID,FRinS1,FRinS2,PutativeSigBinID)
BinnedFiringRateinS1 = []; BinnedFiringRateinS2 = []; BinID = [3 4 5 6]; beforetrial = 1; sample = 1;TimeGain = 10;
DiffBinFRinS1 =[]; DiffBinFRinS2 =[]; tempSigBinID=[]; WorkerNum =10;

%% Binned firing rate,choose only delay bin
for iDelayBin = 1:length(BinID)
    BinnedFiringRateinS1(:,iDelayBin) = mean(FRinS1(:,TimeGain*(beforetrial+sample+iDelayBin-1)+1:TimeGain*(beforetrial+sample+iDelayBin)),2);
    BinnedFiringRateinS2(:,iDelayBin) = mean(FRinS2(:,TimeGain*(beforetrial+sample+iDelayBin-1)+1:TimeGain*(beforetrial+sample+iDelayBin)),2);
end

%% FR difference between Real bins in two sample trials
for iTrial = 1:size(FRinS1,1)  % S1 trials
    for iReferenceBin = 1:length(BinID)
        for iTargetBin = 1:length(BinID)
            DiffBinFRinS1(iReferenceBin,iTargetBin,iTrial) = BinnedFiringRateinS1(iTrial,iReferenceBin)-BinnedFiringRateinS1(iTrial,iTargetBin);
        end
    end
end
for iTrial = 1:size(FRinS2,1)  % S2 trials
    for iReferenceBin = 1:length(BinID)
        for iTargetBin = 1:length(BinID)
            DiffBinFRinS2(iReferenceBin,iTargetBin,iTrial) = BinnedFiringRateinS2(iTrial,iReferenceBin)-BinnedFiringRateinS2(iTrial,iTargetBin);
        end
    end
end
tempZStatistics = zeros(length(BinID),length(BinID));
for iReferenceBin = 1:length(BinID)
    for iTargetBin = 1:length(BinID)
        if iReferenceBin~=iTargetBin
            [~,~,stats] = ranksum(reshape(DiffBinFRinS1(iReferenceBin,iTargetBin,:),1,size(FRinS1,1)),reshape(DiffBinFRinS2(iReferenceBin,iTargetBin,:),1,size(FRinS2,1)));
            tempZStatistics(iReferenceBin,iTargetBin) = stats.zval;
        end
    end
end

%% FR difference between permutaed bins in two sample trials
RepeatedPermutationTimes = 2000;
poolobj = gcp('nocreate'); % if no pool, do not create new one.
if isempty(poolobj)
    myCluster = parcluster('local'); myCluster.NumWorkers = WorkerNum; parpool(myCluster,WorkerNum);
end
ShuffledZStatistics=zeros(length(BinID),length(BinID),RepeatedPermutationTimes);
for iShuffleTimes=1:RepeatedPermutationTimes
    f(iShuffleTimes) = parfeval(@ComputeFRDiffScoreForPermutation,1,BinnedFiringRateinS1,BinnedFiringRateinS2,length(BinID));
end
for iShuffleTimes=1:RepeatedPermutationTimes
    [~,tempShuffledZStatistics] = fetchNext(f);  % collect the results as they become available.
    ShuffledZStatistics(:,:,iShuffleTimes)= tempShuffledZStatistics;
end

%% Permutation test
IsSignificant=zeros(length(BinID),length(BinID));
IsSigLessThanChance=zeros(length(BinID),length(BinID));
tempP=zeros(length(BinID),length(BinID));
for iReferenceBin=1:length(BinID)
    for iTargetBin = 1:length(BinID)
        if tempZStatistics(iReferenceBin,iTargetBin)>=mean(ShuffledZStatistics(iReferenceBin,iTargetBin,:))
            tempP(iReferenceBin,iTargetBin)=sum(ShuffledZStatistics(iReferenceBin,iTargetBin,:)>tempZStatistics(iReferenceBin,iTargetBin))/RepeatedPermutationTimes;
            if tempP(iReferenceBin,iTargetBin)<=0.05
                IsSignificant(iReferenceBin,iTargetBin)=1;
            end
        else
            tempP(iReferenceBin,iTargetBin)=sum(ShuffledZStatistics(iReferenceBin,iTargetBin,:)<tempZStatistics(iReferenceBin,iTargetBin))/RepeatedPermutationTimes;
            if tempP(iReferenceBin,iTargetBin)<=0.05
                IsSigLessThanChance(iReferenceBin,iTargetBin)=1;
            end
        end
        
    end
end
for i=1:length(BinID)
    for j=1:length(BinID)
        if i==j
            tempP(i,j) = 1;
        end
    end
end
PutativeDelaySigBinID = PutativeSigBinID(PutativeSigBinID>=BinID(1) & PutativeSigBinID<=BinID(end));
NoSigBinID = setdiff(BinID,PutativeDelaySigBinID); % get no-selectivity bin ID
SubPutativeDelaySigBinID = [];

%% Investigate whether there are bins with statistically significant selectivity that actually endowed no selectivity
for i=1:length(PutativeDelaySigBinID) % putative sig bin ID
    for j=1:length(NoSigBinID) % no selectivity bin ID
        if tempP(NoSigBinID(j)-2,PutativeDelaySigBinID(i)-2) <=0.05 | tempP(PutativeDelaySigBinID(i)-2,NoSigBinID(j)-2) <=0.05
            SubPutativeDelaySigBinID = [SubPutativeDelaySigBinID PutativeDelaySigBinID(i)];
        end
    end
end
if isempty(SubPutativeDelaySigBinID)
    NoTransientCodingNeuronID = NeuronID;
    PutativeSigBinID(PutativeSigBinID>=BinID(1) & PutativeSigBinID<=BinID(end))=[];
    tempSigBinID = PutativeSigBinID;
else
    NoTransientCodingNeuronID = []; 
    NewNoSigBinID = setdiff(PutativeDelaySigBinID,SubPutativeDelaySigBinID);
    if ~isempty(NewNoSigBinID)
        NoSigBinID = [NoSigBinID NewNoSigBinID];
    end
    
    %% According to no-selectivity bin ID, investigate whether SubPutativeDelaySigBin show significant selectivity
    for i = 1:length(NoSigBinID)  % no selectivity bin ID
        for j = 1:length(SubPutativeDelaySigBinID) % selectivity bin ID
            if (tempP(NoSigBinID(i)-2,SubPutativeDelaySigBinID(j)-2)<=0.05 | tempP(SubPutativeDelaySigBinID(j)-2,NoSigBinID(i)-2)<=0.05)
                tempSigBinID = [tempSigBinID SubPutativeDelaySigBinID(j)];
            end
        end
    end
    tempSigBinID = intersect(PutativeSigBinID,tempSigBinID);
    PutativeSigBinID(PutativeSigBinID>=BinID(1) & PutativeSigBinID<=BinID(end))=[];
    tempSigBinID = sortrows([PutativeSigBinID tempSigBinID]');
    tempSigBinID = tempSigBinID';
end



