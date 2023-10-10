function [AverCorrTriFR,AverErrorTriFR] = CalculateResampledAverCorrErrorFR(CorrTrialsFR,ErrorTrialsFR,BootsStrapNum,SPlen)

MinTriNum = min([size(CorrTrialsFR,1) size(ErrorTrialsFR,1)]);
if MinTriNum <= 20
    MinTriNum = ceil(MinTriNum/2);
else
    MinTriNum = 10*(floor(MinTriNum/10)-1);
end

%% Averaged FR in correct trials
AverCorrTriFR = zeros(BootsStrapNum,SPlen);
for j = 1:BootsStrapNum
    RandPickedTrialNum = randperm(size(CorrTrialsFR,1));
    RandomPickedTrialID = RandPickedTrialNum(1:MinTriNum);
    AverCorrTriFR(j,:) = mean(CorrTrialsFR(RandomPickedTrialID,:),1);
end

%% Averaged FR in error trials
AverErrorTriFR = zeros(BootsStrapNum,SPlen);
for j = 1:BootsStrapNum
    RandPickedTrialNum = randperm(size(ErrorTrialsFR,1));
    RandomPickedTrialID = RandPickedTrialNum(1:MinTriNum);
    AverErrorTriFR(j,:) = mean(ErrorTrialsFR(RandomPickedTrialID,:),1);
end
