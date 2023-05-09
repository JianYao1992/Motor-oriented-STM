function SignificanceLevel = SigLevelAfterPermutation(AUC,ShuffledAUC)

LowConfidenceInterval = prctile(ShuffledAUC,2.5);
HighConfidenceInterval = prctile(ShuffledAUC,97.5);
if AUC < LowConfidenceInterval || AUC > HighConfidenceInterval
    SignificanceLevel = 1;
else
    SignificanceLevel = 0;
end

