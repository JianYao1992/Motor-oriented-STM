function AUC_RemoveNonFcsp = RemoveNonFcspToCalculateAUC(FcspNum_S1,FcspNum_S2,NonFcspNum_S1,NonFcspNum_S2,SpikeNum_S1,SpikeNum_S2)

FcspNum = vertcat(FcspNum_S1,FcspNum_S2);
NonFcspNum = vertcat(NonFcspNum_S1,NonFcspNum_S2);
% random non-FCSP events across trials
TotalFcspNum = sum(FcspNum);
RandSeqNonFcspNum = rand_seq(TotalFcspNum,numel(FcspNum),3*ceil(TotalFcspNum/numel(FcspNum))+2);
RandSeqNonFcspNum = RandSeqNonFcspNum(:);
% permutation
ID = randperm(numel(FcspNum));
PermNonFcspNum = NonFcspNum(ID);
% number of factual non-FCSP events
for i = 1:numel(ID)
    NonFcspNum(ID(i)) = min([RandSeqNonFcspNum(i) PermNonFcspNum(i)]);
end
NonFcspNum_S1 = NonFcspNum(1:numel(FcspNum_S1));
Num_woNonFcsp_S1 = SpikeNum_S1 - NonFcspNum_S1;
NonFcspNum_S2 = NonFcspNum(1+numel(FcspNum_S1):end);
Num_woNonFcsp_S2 = SpikeNum_S2 - NonFcspNum_S2;
% decison variable
[DV_S1,DV_S2] = DecisionVariableCalculation(Num_woNonFcsp_S1,Num_woNonFcsp_S2);
% TPR and FPR
[TPR,FPR] = TprFprCalculation(DV_S1,DV_S2);
% auROC of FCSP events
AUC_RemoveNonFcsp = AucAnalysis(FPR,TPR);

