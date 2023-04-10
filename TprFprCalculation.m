function [TPR,FPR] = TprFprCalculation(DV1,DV2)

% true positive (TP), eqv. with hit
% true negative (TN), eqv. with correct rejection
% false positive (FP),eqv. with false alarm, Type I error
% false negative (FN),eqv. with miss, Type II error

% y axis (Hit rate=hit/(hit+miss))
% sensitivity or true positive rate (TPR), eqv. with hit rate, recall
% TPR = TP/P = TP/(TP+FN);

% x axis (False rate=false alarm/(false alarm+correct rejection))
% fall-out or false positive rate (FPR)(1-Specificity)
% FPR = FP/N = FP/(FP+TN);

ROCCriterion = sortrows([DV1';DV2']);
ROCCriterion1 = ROCCriterion*0.9;
ROCCriterion1 = interp1(1:length(ROCCriterion1),ROCCriterion1,1:0.5:length(ROCCriterion1));
ROCCriterion1 = [ROCCriterion(1)-1;ROCCriterion1';ROCCriterion(end)+1];
TPR = zeros(1,length(ROCCriterion1));
FPR = zeros(1,length(ROCCriterion1));
for i = 1:size(ROCCriterion1,1)
    TPR(1,i) = length(find(DV1>ROCCriterion1(i,1)))/length(DV1);
    FPR(1,i) = length(find(DV2>ROCCriterion1(i,1)))/length(DV2);
end
