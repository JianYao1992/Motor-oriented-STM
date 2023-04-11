function [CorrTriFR,ErrorTriFR,ShuffleCorrTriFR,ShuffleErrorTriFR] = ReSampCorrErrorTrialsForComparison(CorrTrialsFR,ErrorTrialsFR,BootsStrapNum,IsClusterBasedPermuTest,SPlen)

if abs(size(CorrTrialsFR,1)-size(ErrorTrialsFR,1)) < 10
    RealBootsStrapNum = 1;
else
    RealBootsStrapNum = BootsStrapNum;
end
MinTriNum = min([size(CorrTrialsFR,1) size(ErrorTrialsFR,1)]); 

%% Averaged FR in correct trials
CorrTriFR = zeros(RealBootsStrapNum,SPlen-1);
for j = 1:RealBootsStrapNum
    RandPickedTrialNum = randperm(size(CorrTrialsFR,1));
    RandomPickedTrialID = RandPickedTrialNum(1:MinTriNum);
    Sample1CorrTrialFR = mean(CorrTrialsFR(RandomPickedTrialID,1:end-1));
    CorrTriFR(j,:) = Sample1CorrTrialFR;
end
if RealBootsStrapNum > 1
    CorrTriFR = mean(CorrTriFR);
end

%% Averaged FR in error trials
ErrorTriFR = zeros(RealBootsStrapNum,SPlen-1);
for j = 1:RealBootsStrapNum
    RandPickedTrialNum = randperm(size(ErrorTrialsFR,1));
    RandomPickedTrialID = RandPickedTrialNum(1:MinTriNum);
    Sample1ErrorTrialFR = mean(ErrorTrialsFR(RandomPickedTrialID,1:end-1));
    ErrorTriFR(j,:) = Sample1ErrorTrialFR;
end
if RealBootsStrapNum > 1
    ErrorTriFR = mean(ErrorTriFR);
end

%% Shuffle correct and error trials
ShuffleCorrTriFR = zeros(1000,SPlen-1);
ShuffleErrorTriFR = zeros(1000,SPlen-1);
if IsClusterBasedPermuTest == 1        
    AllCorrErrTrialsFR = [CorrTrialsFR; ErrorTrialsFR];    
    for i = 1:1000 % 1000 times of shuffling
        f(i) = parfeval(@ShuffleCorrectErrorTrials,2,AllCorrErrTrialsFR,MinTriNum);
    end
    for i = 1:1000
        [~,ShuffleCorrTrialFR,ShuffleErrorTrialFR] = fetchNext(f);
        ShuffleCorrTriFR(i,:) = ShuffleCorrTrialFR;
        ShuffleErrorTriFR(i,:) = ShuffleErrorTrialFR;
    end
end
