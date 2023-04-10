function [ClusBasedPermuTestIsSig,ClusBasedPermuTestIsSig_Below] = PermutationTest_ClusterBasedOrNot(RealResults,NullDistributions,IsClusBasedPermuTest,ShuffleTimes)

RealResults = mean(RealResults);
TestBinNum = size(RealResults,2);
if IsClusBasedPermuTest == 1 
    %% Significance of bin based on 95% CI for real and shuffle results
    IsSig = zeros(1,TestBinNum);
    IsSig_Below = zeros(1,TestBinNum);
    ShuffleIsSig = zeros(ShuffleTimes,TestBinNum);
    ShuffleIsSig_Below = zeros(ShuffleTimes,TestBinNum);
    for iTestBin = 1:TestBinNum
        temp95Percentile = prctile(NullDistributions(:,iTestBin),[2.5 97.5]);
        ShuffleIsSig(NullDistributions(:,iTestBin) > temp95Percentile(2),iTestBin) = 1; % > 97.5% CI
        ShuffleIsSig_Below(NullDistributions(:,iTestBin) < temp95Percentile(1),iTestBin) = 1; % < 2.5% CI
        if RealResults(iTestBin) > temp95Percentile(2)
            IsSig(1,iTestBin) = 1;
        elseif RealResults(iTestBin) < temp95Percentile(1)
            IsSig_Below(1,iTestBin) = 1;
        end
    end
    
    %% Cluster size for real data
    [RealSigClusSize,RealSigClusBinID] = SigClusterSize(IsSig);
    [RealSigClusSize_Below,RealSigClusBinID_Below] = SigClusterSize(IsSig_Below);
    
    %% Cluster size for shuffled data
    ShuffleSigClusSizeNullDis = zeros(ShuffleTimes,1);
    ShuffleSigClusSizeNullDis_Below = zeros(ShuffleTimes,1);
    for iShuffle = 1:size(NullDistributions,1)
        [tempShuffleSigClusSize,~] = SigClusterSize(ShuffleIsSig(iShuffle,:)); % > 97.5% CI
        if ~isempty(tempShuffleSigClusSize)
            ShuffleSigClusSizeNullDis(iShuffle) = max(tempShuffleSigClusSize);
        end
        [tempShuffleSigClusSize_Below,~] = SigClusterSize(ShuffleIsSig_Below(iShuffle,:)); % < 2.5% CI
        if ~isempty(tempShuffleSigClusSize_Below)
            ShuffleSigClusSizeNullDis_Below(iShuffle) = max(tempShuffleSigClusSize_Below);
        end
    end
    
    %% Significance of real-data bin, after cluster-based permutation test
    % > 97.5% CI
    Shuf95PrctSigClusSize = prctile(ShuffleSigClusSizeNullDis,95); 
    ClusBasedPermuTestIsSig = zeros(1,TestBinNum); 
    if ~isempty(RealSigClusSize)
        for iCluster = 1:length(RealSigClusSize) % check each significant cluster in the real result
            if RealSigClusSize(iCluster) > Shuf95PrctSigClusSize % exceed the 95th percentile of shuffle data
                ClusBasedPermuTestIsSig(1,RealSigClusBinID{iCluster}) = 1;
            end
        end
    end
    % < 2.5% CI
    Shuf95PrctSigClusSize_Below = prctile(ShuffleSigClusSizeNullDis_Below,95);
    ClusBasedPermuTestIsSig_Below = zeros(1,TestBinNum); 
    if ~isempty(RealSigClusSize_Below)
        for iCluster = 1:length(RealSigClusSize_Below) 
            if RealSigClusSize_Below(iCluster) > Shuf95PrctSigClusSize_Below 
                ClusBasedPermuTestIsSig_Below(1,RealSigClusBinID_Below{iCluster}) = 1;
            end
        end
    end
else
    % permutation test
    ClusBasedPermuTestIsSig = zeros(1,TestBinNum);
    ClusBasedPermuTestIsSig_Below = [];
    for iBin = 1:TestBinNum
        p1 = length(find(NullDistributions(:,iBin) <  RealResults(iBin)))./size(NullDistributions, 1); 
        p2 = length(find(NullDistributions(:,iBin) >= RealResults(iBin)))./size(NullDistributions, 1);
        if min([p1 p2]) < 0.05
            ClusBasedPermuTestIsSig(iBin) = 1;
        end
    end
end