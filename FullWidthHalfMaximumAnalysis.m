function [Result,LeftPoint,RightPoint] = FullWidthHalfMaximumAnalysis(Wave)
Threshold = min(Wave)/2;
LeftLeanCrossPoint = find([Wave 10000]<Threshold & [-10000 Wave]>Threshold);
RightLeanCrossPoint = find([Wave -10000]>Threshold & [10000 Wave]<Threshold);
[~,I] = min(Wave);
LeftCrossPoint = LeftLeanCrossPoint(LeftLeanCrossPoint<I);
RightCrossPoint = RightLeanCrossPoint(RightLeanCrossPoint>I);
if ~isempty(LeftCrossPoint) && ~isempty(RightCrossPoint)
    LeftPoint = max(LeftCrossPoint);
    RightPoint = min(RightCrossPoint);
    Result = (RightPoint-(Wave(RightPoint)-Threshold)/((Wave(RightPoint)-Wave(RightPoint-1)))) - (LeftPoint-((Wave(LeftPoint)-Threshold)/(Wave(LeftPoint)-Wave(LeftPoint-1))));
elseif ~isempty(RightCrossPoint)
    RightPoint = RightCrossPoint(1);
    Result = ((RightPoint-((Wave(RightPoint)-Threshold)/(Wave(RightPoint)-Wave(RightPoint-1))))-I)*2;
elseif ~isempty(LeftCrossPoint)
    LeftPoint = LeftCrossPoint(length(LeftCrossPoint));
    Result = (I-(LeftPoint-((Wave(LeftPoint)-Threshold)/(Wave(LeftPoint)-Wave(LeftPoint-1)))))*2;
else
    Result=length(Wave);
end