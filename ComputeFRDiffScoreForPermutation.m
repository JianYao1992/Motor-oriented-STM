function  zvalue_perm = ComputeFRDiffScoreForPermutation(BinnedFiringRateinS1,BinnedFiringRateinS2,BinNum)

    for iTrial = 1:size(BinnedFiringRateinS1,1) % trial
        tempTimeBin = randperm(BinNum);
        for iReferenceBin = 1:BinNum
            for iTargetBin = 1:BinNum
                DiffBinFRinS1(iReferenceBin,iTargetBin,iTrial) = BinnedFiringRateinS1(iTrial,tempTimeBin(iReferenceBin))-BinnedFiringRateinS1(iTrial,tempTimeBin(iTargetBin));
            end
        end
    end
    for iTrial = 1:size(BinnedFiringRateinS2,1)
        tempTimeBin = randperm(BinNum);
        for iReferenceBin = 1:BinNum
            for iTargetBin = 1:BinNum
                DiffBinFRinS2(iReferenceBin,iTargetBin,iTrial) = BinnedFiringRateinS2(iTrial,tempTimeBin(iReferenceBin))-BinnedFiringRateinS2(iTrial,tempTimeBin(iTargetBin));
            end
        end
    end
    zvalue_perm = zeros(BinNum,BinNum);
    for iReferenceBin = 1:BinNum
        for iTargetBin = 1:BinNum
            if iReferenceBin~=iTargetBin
                [~,~,stats] = ranksum(reshape(DiffBinFRinS1(iReferenceBin,iTargetBin,:),1,size(BinnedFiringRateinS1,1)),reshape(DiffBinFRinS2(iReferenceBin,iTargetBin,:),1,size(BinnedFiringRateinS2,1)));
                zvalue_perm(iReferenceBin,iTargetBin) = stats.zval;
            end
        end
    end