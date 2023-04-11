function NormalizedFR = NormFrCalculation(RawFR,BaselinePer)

BaselineMean = mean(RawFR(:,BaselinePer),2);
BaselineSTD = std(RawFR(:,BaselinePer),[],2);
NormalizedFR = (RawFR-repmat(BaselineMean,1,size(RawFR,2)))./repmat(BaselineSTD,1,size(RawFR,2));
if ~isempty(isnan(NormalizedFR))
    NanID = +isnan(NormalizedFR(:,1));
    NanID = NanID==1;
    NormalizedFR(NanID,:) = [];
end