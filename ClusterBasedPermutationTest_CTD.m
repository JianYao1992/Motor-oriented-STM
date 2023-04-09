function [IsSig,ClusBasedPermuTestIsSig,ShuffleIsSig,ShuffleSigClusSize,RealDecodingResultsSigClusSize,...
    RealDecodingResultsSigClusID,IsSigLessThanChance,ClusBasedPermuTestIsSigLessChance,RealDecodingResultsSigLessClusSize]...
    = ClusterBasedPermutationTest_CTD(RealDecodingResults,ShuffleTCTDecodingNullDistribution,TestBinNum,ShuffleDecodingTimes,SingleOrTwoSidesTest)

if nargin == 4
    IsSingleOrTwoSidesTest = 1;
else
    IsSingleOrTwoSidesTest = SingleOrTwoSidesTest;
end
IsSig = zeros(TestBinNum,TestBinNum);
ClusBasedPermuTestIsSig = zeros(TestBinNum,TestBinNum);
IsSigLessThanChance = zeros(TestBinNum,TestBinNum);
ClusBasedPermuTestIsSigLessChance = zeros(TestBinNum,TestBinNum);
ShuffleIsSig = zeros(ShuffleDecodingTimes,TestBinNum,TestBinNum);
ShuffleSigClusSize = zeros(TestBinNum,ShuffleDecodingTimes);
RealDecodingResultsSigClusSize = cell(TestBinNum,1);
RealDecodingResultsSigClusID = cell(TestBinNum,1);
RealDecodingResultsSigLessClusSize = cell(TestBinNum,1);
for iTrainBin=1:TestBinNum
    for iTestBin = 1:TestBinNum
        %% Permutation test
        if IsSingleOrTwoSidesTest == 1 % single-side test
            temp95Percentile = prctile(ShuffleTCTDecodingNullDistribution(:,iTrainBin,iTestBin),95);
            ShuffleIsSig(ShuffleTCTDecodingNullDistribution(:,iTrainBin,iTestBin)>temp95Percentile,iTrainBin,iTestBin) = 1;
            if RealDecodingResults(iTrainBin,iTestBin) > temp95Percentile
                IsSig(iTrainBin,iTestBin) = 1;
            end
        else % two-sides test
            % shuffled data
            temp95Percentile = prctile(ShuffleTCTDecodingNullDistribution(:,iTrainBin,iTestBin),[97.5 2.5]);
            ShuffleIsSig(ShuffleTCTDecodingNullDistribution(:,iTrainBin,iTestBin)>temp95Percentile(1),iTrainBin,iTestBin) = 1;
            ShuffleIsSig(ShuffleTCTDecodingNullDistribution(:,iTrainBin,iTestBin)<temp95Percentile(2),iTrainBin,iTestBin) = -1;
            % real data
            if RealDecodingResults(iTrainBin,iTestBin) > temp95Percentile(1)
                IsSig(iTrainBin,iTestBin) = 1;
            elseif RealDecodingResults(iTrainBin,iTestBin) < temp95Percentile(2)
                IsSigLessThanChance(iTrainBin,iTestBin) = -1;
            end
        end
    end
    %% Size of cluster for decoding result of real data, at specific training point
    % higher than chance
    [tempRealDecodingResultsSigClusSize,tempRealDecodingResultsSigClusD] = SigClusterSize(IsSig(iTrainBin,:));
    RealDecodingResultsSigClusSize{iTrainBin} = tempRealDecodingResultsSigClusSize;
    RealDecodingResultsSigClusID{iTrainBin} = tempRealDecodingResultsSigClusD;
    % lower than chance
    [tempRealDecodingResultsSigLessClusSize,tempRealDecodingResultsSigLessClusID] = SigClusterSize(IsSigLessThanChance(iTrainBin,:));
    RealDecodingResultsSigLessClusSize{iTrainBin} = tempRealDecodingResultsSigLessClusSize;
    
    %% Maximal size of cluster for decoding result of shuffled data, at specific training point
    tempTrainAllShuffleTestIsSig = ShuffleIsSig(:,iTrainBin,:);
    Size = size(tempTrainAllShuffleTestIsSig);
    tempTrainAllShuffleTestIsSig = reshape(tempTrainAllShuffleTestIsSig,Size(1),Size(3));
    for iShuffle = 1:ShuffleDecodingTimes
        ShuffletempTrainIsSignificant = tempTrainAllShuffleTestIsSig(iShuffle,:);
        [tempAllSigClusterSize,~] = SigClusterSize(ShuffletempTrainIsSignificant);
        if ~isempty(tempAllSigClusterSize)
            ShuffleSigClusSize(iTrainBin,iShuffle) = max(tempAllSigClusterSize);
        end
    end
    
    %% Cluster-based permutation test
    % higher than chance
    ShuffleSigClusSizeNullDis = ShuffleSigClusSize(iTrainBin,:);
    tempShuffle95PercentileClusSize = prctile(ShuffleSigClusSizeNullDis,95);
    if ~isempty(tempRealDecodingResultsSigClusSize)
        for iCluster = 1:length(tempRealDecodingResultsSigClusSize)
            if tempRealDecodingResultsSigClusSize(iCluster) > tempShuffle95PercentileClusSize
                ClusBasedPermuTestIsSig(iTrainBin,tempRealDecodingResultsSigClusD{iCluster}) = 1;
            end
        end
    end
    % lower than chance
    if IsSingleOrTwoSidesTest ~= 1
        if ~isempty(tempRealDecodingResultsSigLessClusSize)
            for iCluster = 1:length(tempRealDecodingResultsSigLessClusSize)
                if tempRealDecodingResultsSigLessClusSize(iCluster) > tempShuffle95PercentileClusSize
                    ClusBasedPermuTestIsSigLessChance(iTrainBin,tempRealDecodingResultsSigLessClusID{iCluster}) = 1;
                end
            end
        end
    end
end