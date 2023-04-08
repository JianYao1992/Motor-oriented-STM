function SortedID = SortingByDiffBetweenTwoSampleTrials(FRinS1,FRinS2,ArraPeriod,IsArra)

meanFRinS1 = cellfun(@mean,FRinS1,'UniformOutput', 0);
meanFRinS2 = cellfun(@mean,FRinS2,'UniformOutput', 0);
MeanFRDiff = cellfun(@(x,y) x - y,meanFRinS1,meanFRinS2,'UniformOutput', 0);
MeanFRDiff = vertcat(MeanFRDiff{:});
Selec = cellfun(@(x,y) (x - y)./(x + y),meanFRinS1,meanFRinS2,'UniformOutput', 0);
Selec = vertcat(Selec{:});
if IsArra == 1
    ArraSelec = sortrows([-1*mean(MeanFRDiff(:,ArraPeriod),2) [1:size(MeanFRDiff,1)]' Selec],1);
    SortedID = ArraSelec(:,2);
    SubArraSelec = ArraSelec(:,3:end);
end
SmoothSelec = [];
for i = 1:size(SubArraSelec,1)
    SmoothSelec = [SmoothSelec; smooth(SubArraSelec(i,:),3)'];
end
figure;
colormap('jet');
imagesc(SmoothSelec,[-0.5 0.5]);
box('off');
unitnum = size(SmoothSelec,1);
hold on
plot([10 10],[0 unitnum],'color','k','linestyle','--','linewidth',3); % sample odor
plot([20 20],[0 unitnum],'color','k','linestyle','--','linewidth',3);
plot([60 60],[0 unitnum],'color','k','linestyle','--','linewidth',3); % test odor
plot([65 65],[0 unitnum],'color','k','linestyle','--','linewidth',3);
set(gca,'XTick',[10 20 60 65]);
set(gca,'XTickLabel',{'0','1','5','5.5'},'FontName','Arial','FontSize',16);
set(gca,'YTick',zeros(1,0));
set(gca,'YTickLabel',{''},'FontName','Arial','FontSize',16);
xlabel('Time (sec)','FontSize',18,'FontName','Arial');
box off;
