function [NormShuffleCorrTriFR,NormShuffleErrorTriFR] = NormFrCalculationForShuffledData(ShuffleCorrTriFR,ShuffleErrorTriFR,SPlen,BaselinePer)

ShuffleTimes = size(ShuffleCorrTriFR,2);
ShuffleCorrTriFR = reshape(ShuffleCorrTriFR,ShuffleTimes,SPlen-1);
ShuffleErrorTriFR = reshape(ShuffleErrorTriFR,ShuffleTimes,SPlen-1);

MeanFR = (ShuffleCorrTriFR+ShuffleErrorTriFR)/2;
BaselineMean = mean(MeanFR(:,BaselinePer),2);
BaselineSTD = std(MeanFR(:,BaselinePer),[],2);

BinNum = size(ShuffleCorrTriFR,2);
NormShuffleCorrTriFR = (ShuffleCorrTriFR-repmat(BaselineMean,1,BinNum))./repmat(BaselineSTD,1,BinNum);
NormShuffleErrorTriFR = (ShuffleErrorTriFR-repmat(BaselineMean,1,BinNum))./repmat(BaselineSTD,1,BinNum);
