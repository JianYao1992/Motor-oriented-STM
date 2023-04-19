function FcDensity = CalculateSessionBasedFcDensity(TotalPairsNumPerSecond,FcPairsNumPerSec,PairsNumthres,BinNumCriteria)

FcDensity = [];
for iSess = 1:size(TotalPairsNumPerSecond,1) % session
    TarBinID = find(TotalPairsNumPerSecond(iSess,:)>=PairsNumthres);
    if length(TarBinID) >= BinNumCriteria
        temp_Hit = FcPairsNumPerSec(iSess,TarBinID)./TotalPairsNumPerSecond(iSess,TarBinID);
        FcDensity(end+1,:) = mean(temp_Hit);
    end
end