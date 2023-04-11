function [ShuffleCorrTrialFR,ShuffleErrorTrialFR] = ShuffleCorrectErrorTrials(AllCorrErrTrials,MinTriNum)

TotalTrialNum = size(AllCorrErrTrials,1);
RandomTrialID = randperm(TotalTrialNum);
RandomPickedCorrTrialID = RandomTrialID(1:MinTriNum);
RandomPickedErrorTrialID = RandomTrialID(MinTriNum+1:MinTriNum*2);
if size(AllCorrErrTrials,2) > 1
    ShuffleCorrTrialFR = mean(AllCorrErrTrials(RandomPickedCorrTrialID,1:end-1));
    ShuffleErrorTrialFR = mean(AllCorrErrTrials(RandomPickedErrorTrialID,1:end-1));
else
    ShuffleCorrTrialFR = mean(AllCorrErrTrials(RandomPickedCorrTrialID));
    ShuffleErrorTrialFR = mean(AllCorrErrTrials(RandomPickedErrorTrialID));
end

