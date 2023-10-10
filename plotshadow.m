function plotshadow(data,C,IsSTDorSEM,SmoothWindow,Xshift,TimeGain)

x = (1:size(data,2))/TimeGain+Xshift;
average = mean(data,1);
plot(x,smooth(average,SmoothWindow)','color',C,'linewidth',2);
hold on
if IsSTDorSEM == 2  % 1:STD; 2:SEM; 3:95% CI
    Highervalue = (smooth(average,SmoothWindow))' + std(data,0,1)/sqrt(size(data,1));
    Lowervalue = (smooth(average,SmoothWindow))' - std(data,0,1)/sqrt(size(data,1));
elseif IsSTDorSEM == 1
    Highervalue = (smooth(average,SmoothWindow))' + std(data,0,1);
    Lowervalue = (smooth(average,SmoothWindow))' - std(data,0,1);
elseif IsSTDorSEM == 3
    CI = prctile(data,[2.5 97.5],1);
    Highervalue = CI(2,:);
    Highervalue = smooth(Highervalue,SmoothWindow);
    Lowervalue = CI(1,:);
    Lowervalue = smooth(Lowervalue,SmoothWindow);
end
Time = [x, fliplr(x)];
value = [Highervalue(:)',fliplr(Lowervalue(:)')];
a = fill(Time,value,C,'edgecolor','none');
alpha(a,0.2);
box off;
hold on

