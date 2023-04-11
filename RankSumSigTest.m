function [p,h,MeanFR1,MeanFR2,SelectivityIndex] = RankSumSigTest(Trials1,Trials2,EventPeriod,EventPeriod2)

TargetFR1 = mean(Trials1(:,EventPeriod),2);
TargetFR2 = mean(Trials2(:,EventPeriod2),2);

[p,h] = ranksum(TargetFR1,TargetFR2);

MeanFR1 = mean(TargetFR1);
MeanFR2 = mean(TargetFR2);
SelectivityIndex = (MeanFR1-MeanFR2)/(MeanFR1+MeanFR2);