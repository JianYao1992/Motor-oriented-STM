function PlotWellTrainedResults(perfGroup1,perfGroup2,interval,popinterval,OptoGroup,IsPlotSingleMouse,IsControlBarLength,IsHitTrials)
if ~isempty(regexp(OptoGroup,'NpHR'))
    C = [0 125 0]/255;
else
    C = [67 106 178]/255;
end
if IsPlotSingleMouse == 1
    NumGroup1 = size(perfGroup1,1); NumGroup2 = size(perfGroup2,1);
    XPosGroup1Off = (1-interval*floor(NumGroup1/2)):interval:(1-interval*floor(NumGroup1/2)+interval*(NumGroup1-1));
    XPosGroup1On = (2-interval*floor(NumGroup1/2)):interval:(2-interval*floor(NumGroup1/2)+interval*(NumGroup1-1));
    XPosGroup2Off = (4-interval*floor(NumGroup2/2)):interval:(4-interval*floor(NumGroup2/2)+interval*(NumGroup2-1));
    XPosGroup2On = (5-interval*floor(NumGroup2/2)):interval:(5-interval*floor(NumGroup2/2)+interval*(NumGroup2-1));
    for k = 1:size(perfGroup1,1)
        plot([XPosGroup1Off(k) XPosGroup1On(k)],[perfGroup1(k,1) perfGroup1(k,2)],'k','LineWidth',2,'marker','o','markerfacecolor','k','markeredgecolor','none');
        hold on
    end
    for k = 1:size(perfGroup2,1)
        plot([XPosGroup2Off(k) XPosGroup2On(k)],[perfGroup2(k,1) perfGroup2(k,2)],'color',C,'LineWidth',2,'marker','o','markerfacecolor',C,'markeredgecolor','none');
        hold on
    end
    % population performance
    ci1 = bootci(1000,@mean,perfGroup1(:,1));
    ci2 = bootci(1000,@mean,perfGroup1(:,2));
    ci3 = bootci(1000,@mean,perfGroup2(:,1));
    ci4 = bootci(1000,@mean,perfGroup2(:,2));
    errorbar(1-interval*floor(NumGroup1/2)-popinterval,mean(perfGroup1(:,1)),abs(ci1(1)-mean(perfGroup1(:,1))),abs(ci1(2)-mean(perfGroup1(:,1))),'linestyle','none','color','k','marker','o','markerfacecolor','w','markeredgecolor','k','markersize',10);
    hold on
    errorbar(2-interval*floor(NumGroup1/2)+interval*(NumGroup1-1)+popinterval,mean(perfGroup1(:,2)),abs(ci2(1)-mean(perfGroup1(:,2))),abs(ci2(2)-mean(perfGroup1(:,2))),'linestyle','none','color','k','marker','o','markerfacecolor','k','markeredgecolor','k','markersize',10);
    hold on
    errorbar(4-interval*floor(NumGroup2/2)-popinterval,mean(perfGroup2(:,1)),abs(ci3(1)-mean(perfGroup2(:,1))),abs(ci3(2)-mean(perfGroup2(:,1))),'linestyle','none','color',C,'marker','o','markerfacecolor','w','markeredgecolor',C,'markersize',10);
    hold on
    errorbar(5-interval*floor(NumGroup2/2)+interval*(NumGroup2-1)+popinterval,mean(perfGroup2(:,2)),abs(ci4(1)-mean(perfGroup2(:,2))),abs(ci4(2)-mean(perfGroup2(:,2))),'linestyle','none','color',C,'marker','o','markerfacecolor',C,'markeredgecolor',C,'markersize',10);
    set(gca,'XTick',[1 2 4 5],'XTickLabel',{'','','',''},'xlim',[0 6]);
    set(gca,'YTick',0:1:4,'YTickLabel',{'0','1','2','3','4'},'FontName','Arial','FontSize',16,'ylim',[0 4]);
else
    NumGroup1 = size(perfGroup1,1); NumGroup2 = size(perfGroup2,1);
    % population performance
    ci1 = bootci(1000,@mean,perfGroup1(:,1));
    ci2 = bootci(1000,@mean,perfGroup1(:,2));
    ci3 = bootci(1000,@mean,perfGroup2(:,1));
    ci4 = bootci(1000,@mean,perfGroup2(:,2));
    if IsControlBarLength == 1
        errorbar(0.2,mean(perfGroup1(:,1)),abs(ci1(1)-mean(perfGroup1(:,1))),abs(ci1(2)-mean(perfGroup1(:,1))),'linestyle','none','color','k','marker','o','markerfacecolor','w','markeredgecolor','k','markersize',10); hold on
        errorbar(2.2,mean(perfGroup1(:,2)),abs(ci2(1)-mean(perfGroup1(:,2))),abs(ci2(2)-mean(perfGroup1(:,2))),'linestyle','none','color','k','marker','o','markerfacecolor','k','markeredgecolor','k','markersize',10); hold on
        errorbar(0.8,mean(perfGroup2(:,1)),abs(ci3(1)-mean(perfGroup2(:,1))),abs(ci3(2)-mean(perfGroup2(:,1))),'linestyle','none','color',C,'marker','o','markerfacecolor','w','markeredgecolor',C,'markersize',10); hold on
        errorbar(2.8,mean(perfGroup2(:,2)),abs(ci4(1)-mean(perfGroup2(:,2))),abs(ci4(2)-mean(perfGroup2(:,2))),'linestyle','none','color',C,'marker','o','markerfacecolor',C,'markeredgecolor',C,'markersize',10); hold on
        if IsHitTrials == 1
            plot([0.2 2.2],[mean(perfGroup1(:,1)) mean(perfGroup1(:,2))],'--k'); hold on
            plot([0.8 2.8],[mean(perfGroup2(:,1)) mean(perfGroup2(:,2))],'linestyle','--','color',C); hold on
        else
            plot([0.2 2.2],[mean(perfGroup1(:,1)) mean(perfGroup1(:,2))],'k'); hold on
            plot([0.8 2.8],[mean(perfGroup2(:,1)) mean(perfGroup2(:,2))],'color',C); hold on
        end
    else
        errorbar(1:2,mean(perfGroup1),std(perfGroup1)/sqrt(size(perfGroup1,1)),'linestyle','none','color','k','marker','o','markerfacecolor','k','markeredgecolor','k','markersize',10); hold on
        errorbar(1:2,mean(perfGroup2),std(perfGroup2)/sqrt(size(perfGroup2,1)),'linestyle','none','color',C,'marker','o','markerfacecolor',C,'markeredgecolor',C,'markersize',10); hold on
        if IsHitTrials == 1
            plot([1 2],mean(perfGroup1),'--k'); hold on
            plot([1 2],mean(perfGroup2),'linestyle','--','color',C); hold on
        else
            plot([1 2],mean(perfGroup1),'k'); hold on
            plot([1 2],mean(perfGroup2),'color',C); hold on
        end
    end
    set(gca,'XTick',[1 2],'xlim',[0.8 2.2]);
end
