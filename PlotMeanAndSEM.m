function [YMin,YMax] = PlotMeanAndSEM(Data,X,Color,PlotRange,LineStyle,Area,STD)

bin_width = 5;
if  nargin == 7
    StdData = STD;
    IsArea = Area;
    bin_width = 2;
elseif nargin == 6
    StdData = std(Data)/sqrt(size(Data,1));
    IsArea = Area;
    bin_width = 2;
elseif nargin == 5
    StdData = std(Data)/sqrt(size(Data,1));
    IsArea = 1;    
    bin_width = 3;
elseif nargin == 4
    StdData = std(Data);
    LineStyle = '-';
    IsArea = 1;    
    bin_width = 2;
elseif nargin == 3
    StdData = std(Data);
    PlotRange = [min(X) max(X)];
    LineStyle = '-';
    IsArea = 1;
end
if size(Data,1) > 1
    AveragedData = smooth(mean(Data),bin_width)';
else
    AveragedData = Data;
end
StdData = StdData(:)';
PlotIndex = find(X>=PlotRange(1) & X<=PlotRange(2));
AveragedData = AveragedData(PlotIndex);
StdData = StdData(:,PlotIndex);
X = X(PlotIndex);
X2 = fliplr(X);
if size(StdData,1)==1 % STD or SEM
    y1 = AveragedData+StdData;
    y2 = AveragedData-StdData;
elseif size(StdData,1) == 2 % 95% CI
    y1 = StdData(1,:);
    y2 = StdData(2,:);
end
x0 = [X X2];
y0 = [y1 fliplr(y2)];
if IsArea == 1
    fill(x0,y0,Color,'edgecolor',[1 1 1],'FaceAlpha',0.5,'linestyle','none');
end
hold on
plot(X,AveragedData,'color',Color,'linestyle',LineStyle,'linewidth',1);
YMin = min(y2);
YMax = max(y1);
