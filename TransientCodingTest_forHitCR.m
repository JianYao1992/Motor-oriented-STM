%% Judge whether single neuron is transient neuron that has significant selectivity for hit vs. CR
% input data is FR of individual 'transient' neuron

function [NoTransientCodingNeuronID,tempSigBinID,tempZStatistics,tempP] = TransientCodingTest_forHitCR(NeuronID,FR_Hit,FR_CR,PutativeSigBinID,EleBinsNum,BaseLen,ShownBaseLen,SampOdorLen,DelayLen,TestOdorLen,RespLen)

BinFR_Hit = []; 
BinFR_CR = []; 
BinFRDiff_Hit =[]; 
BinFRDiff_CR =[];
BinsRange = (BaseLen-ShownBaseLen)*10/EleBinsNum+1:1:(BaseLen+SampOdorLen+DelayLen+TestOdorLen+RespLen)*10/EleBinsNum;
tempSigBinID=[];
WorkerNum =10;

%% Binned FR, choosing bins in range
for iBin = 1:length(BinsRange)
    BinFR_Hit(:,iBin) = mean(FR_Hit(:,(BinsRange(iBin)-1)*EleBinsNum+2:BinsRange(iBin)*EleBinsNum+1),2);
    BinFR_CR(:,iBin) = mean(FR_CR(:,(BinsRange(iBin)-1)*EleBinsNum+2:BinsRange(iBin)*EleBinsNum+1),2);
end

%% Difference in FR across bins in hit and CR trials
for iTrial = 1:size(FR_Hit,1) % hit trials
    for iReferenceBin = 1:length(BinsRange)
        for iTargetBin = 1:length(BinsRange)
            BinFRDiff_Hit(iReferenceBin,iTargetBin,iTrial) = BinFR_Hit(iTrial,iReferenceBin) - BinFR_Hit(iTrial,iTargetBin);
        end
    end
end
for iTrial = 1:size(FR_CR,1)  % CR trials
    for iReferenceBin = 1:length(BinsRange)
        for iTargetBin = 1:length(BinsRange)
            BinFRDiff_CR(iReferenceBin,iTargetBin,iTrial) = BinFR_CR(iTrial,iReferenceBin)-BinFR_CR(iTrial,iTargetBin);
        end
    end
end
tempZStatistics = zeros(length(BinsRange),length(BinsRange));
for iReferenceBin = 1:length(BinsRange)
    for iTargetBin = 1:length(BinsRange)
        if iReferenceBin~=iTargetBin
            [~,~,stats] = ranksum(reshape(BinFRDiff_Hit(iReferenceBin,iTargetBin,:),1,size(FR_Hit,1)),reshape(BinFRDiff_CR(iReferenceBin,iTargetBin,:),1,size(FR_CR,1)));
            tempZStatistics(iReferenceBin,iTargetBin) = stats.zval;
        end
    end
end

%% Difference in FR across permutaed bins in Hit and CR trials
RepeatedPermutationTimes = 2000;
poolobj = gcp('nocreate'); % if no pool, do not create new one.
if isempty(poolobj)
    myCluster = parcluster('local'); myCluster.NumWorkers = WorkerNum; parpool(myCluster,WorkerNum);
end
ShuffledZStatistics = zeros(length(BinsRange),length(BinsRange),RepeatedPermutationTimes);
for iShuffleTimes=1:RepeatedPermutationTimes
    f(iShuffleTimes) = parfeval(@ComputeFRDiffScoreForPermutation,1,BinFR_Hit,BinFR_CR,length(BinsRange));
end
for iShuffleTimes=1:RepeatedPermutationTimes
    [~,tempShuffledZStatistics] = fetchNext(f);  % collect the results as they become available.
    ShuffledZStatistics(:,:,iShuffleTimes)= tempShuffledZStatistics;
end

%% Permutation test
IsSignificant = zeros(length(BinsRange),length(BinsRange));
IsSigLessThanChance = zeros(length(BinsRange),length(BinsRange));
tempP = zeros(length(BinsRange),length(BinsRange));
for iReferenceBin = 1:length(BinsRange)
    for iTargetBin = 1:length(BinsRange)
        if tempZStatistics(iReferenceBin,iTargetBin) >= mean(ShuffledZStatistics(iReferenceBin,iTargetBin,:))
            tempP(iReferenceBin,iTargetBin) = sum(ShuffledZStatistics(iReferenceBin,iTargetBin,:)>tempZStatistics(iReferenceBin,iTargetBin))/RepeatedPermutationTimes;
            if tempP(iReferenceBin,iTargetBin) <= 0.05
                IsSignificant(iReferenceBin,iTargetBin) = 1;
            end
        else
            tempP(iReferenceBin,iTargetBin) = sum(ShuffledZStatistics(iReferenceBin,iTargetBin,:)<tempZStatistics(iReferenceBin,iTargetBin))/RepeatedPermutationTimes;
            if tempP(iReferenceBin,iTargetBin) <= 0.05
                IsSigLessThanChance(iReferenceBin,iTargetBin) = 1;
            end
        end
        
    end
end
for i=1:length(BinsRange)
    for j=1:length(BinsRange)
        if i==j
            tempP(i,j) = 1;
        end
    end
end
PutativeSigBinID_inRange = PutativeSigBinID(PutativeSigBinID>=BinsRange(1) & PutativeSigBinID<=BinsRange(end));
NoSigBinID = setdiff(BinsRange,PutativeSigBinID_inRange); % ID of non-selective bin
SubPutativeSigBinID_inRange = [];

%% Investigate whether there is bin with significant selectivity that don't show selectivity
for i=1:length(PutativeSigBinID_inRange) % ID of 'significant' bin
    for j=1:length(NoSigBinID) % ID of bin showing no selectivity
        if tempP(find(BinsRange==NoSigBinID(j)),find(BinsRange==PutativeSigBinID_inRange(i))) <= 0.05 | tempP(find(BinsRange==PutativeSigBinID_inRange(i)),find(BinsRange==NoSigBinID(j))) <= 0.05
            SubPutativeSigBinID_inRange = [SubPutativeSigBinID_inRange PutativeSigBinID_inRange(i)];
        end
    end
end
if isempty(SubPutativeSigBinID_inRange)
    NoTransientCodingNeuronID = NeuronID;
    PutativeSigBinID(PutativeSigBinID>=BinsRange(1) & PutativeSigBinID<=BinsRange(end))=[];
    tempSigBinID = PutativeSigBinID;
else
    NoTransientCodingNeuronID = [];  
    NewNoSigBinID = setdiff(PutativeSigBinID_inRange,SubPutativeSigBinID_inRange);
    if ~isempty(NewNoSigBinID)
        NoSigBinID = [NoSigBinID NewNoSigBinID];
    end
    
    %% According to ID of bin showing no selectivity, investigate whether SubPutativeDelaySigBin show significant selectivity
    for i = 1:length(NoSigBinID)  % ID of bin showing no selectivity
        for j = 1:length(SubPutativeSigBinID_inRange) % ID of bin with significant selectivity
            if (tempP(find(BinsRange==NoSigBinID(i)),find(BinsRange==SubPutativeSigBinID_inRange(j)))<=0.05 | tempP(find(BinsRange==SubPutativeSigBinID_inRange(j)),find(BinsRange==NoSigBinID(i)))<=0.05)
                tempSigBinID = [tempSigBinID SubPutativeSigBinID_inRange(j)];
            end
        end
    end
    tempSigBinID = intersect(PutativeSigBinID,tempSigBinID);
    PutativeSigBinID(PutativeSigBinID>=BinsRange(1) & PutativeSigBinID<=BinsRange(end))=[];
    tempSigBinID = sortrows([PutativeSigBinID tempSigBinID]');
    tempSigBinID = tempSigBinID';
end



