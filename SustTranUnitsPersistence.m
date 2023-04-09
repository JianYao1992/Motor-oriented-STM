%% Persistence of cross-temporal decoding between sustained and transient neurons

clear; clc; close all;

%% Assignment
TargetBrain = {'mPFC','aAIC'};
TimeGain = 10;


p = [];
ClusBasedPermuTestIsSigBinClusSize = cell(numel(TargetBrain),1);
for ireg = 1:numel(TargetBrain)
    if strcmp(TargetBrain{ireg},'mPFC')
        PersistenceResults = dir('*mPFC*ClusterBased*.mat');
    elseif strcmp(TargetBrain{ireg},'aAIC')
        PersistenceResults = dir('*aAIC*ClusterBased*.mat');
    end
    for iFile=1:size(PersistenceResults,1)
        load(PersistenceResults(iFile).name);
        if ~isempty(regexp(PersistenceResults(iFile).name,'sustained'))
            ClusterBasedPermuTestIsSig_Sust = ClusterBasedPermuTestIsSignificant;
        else
            ClusterBasedPermuTestIsSig_Tran = ClusterBasedPermuTestIsSignificant;
        end
    end
    clearvars ClusterBasedPermuTestIsSignificant ClusterBasedPermuTestIsSigLessThanChance
    X = 6:60;
    ClusBasedPermuTestIsSig{1,1} = ClusterBasedPermuTestIsSig_Sust(X,X);
    ClusBasedPermuTestIsSig{1,2} = ClusterBasedPermuTestIsSig_Tran(X,X);
    
    %% Number of consecutive bins with significant decoding accuracy for each train time
    for iGroup = 1:numel(ClusBasedPermuTestIsSig)
        MaxSigBinNum = [];
        for iBin = 1:size(ClusBasedPermuTestIsSig{1,iGroup},1) % train time
            tempSigBinNum = 0;
            SigBinNum = [];
            for j = 1:size(ClusBasedPermuTestIsSig{1,iGroup},2) % test time
                if ClusBasedPermuTestIsSig{1,iGroup}(iBin,j) == 1
                    tempSigBinNum = tempSigBinNum + 1;
                    if j == size(ClusBasedPermuTestIsSig{1,iGroup},2)
                        SigBinNum = [SigBinNum tempSigBinNum];
                    end
                else
                    if ClusBasedPermuTestIsSig{1,iGroup}(iBin,j) == 0
                        SigBinNum = [SigBinNum tempSigBinNum];
                    end
                    tempSigBinNum = 0;
                end
            end
            if isempty(SigBinNum)
                MaxSigBinNum(iBin) = 0;
            else
                MaxSigBinNum(iBin) = max(SigBinNum);
            end
        end
        ClusBasedPermuTestIsSigBinClusSize{ireg}{1,iGroup} = MaxSigBinNum;
    end
    
    %% Wilcoxon signed-rank test
    p(1,end+1) = signrank(ClusBasedPermuTestIsSigBinClusSize{ireg}{1},ClusBasedPermuTestIsSigBinClusSize{ireg}{2});
end

%% Bar plot
figure('OuterPosition',[219 303 420 534]);
for ireg = 1:numel(ClusBasedPermuTestIsSigBinClusSize)
    tempPersistence = cellfun(@(x) x/10,ClusBasedPermuTestIsSigBinClusSize{ireg},'UniformOutput', 0);
    tempMeanPersistence = cellfun(@mean,tempPersistence);
    bar([1.1+1.6*(ireg-1) 1.9+1.6*(ireg-1)],[tempMeanPersistence(1) tempMeanPersistence(2)],0.6,'FaceColor','none','EdgeColor','k');
    hold on
    % individual persistence in each train bin
    for iBin = 1:length(tempPersistence{1})
        plot([1.1+1.6*(ireg-1) 1.9+1.6*(ireg-1)],[tempPersistence{1}(iBin) tempPersistence{2}(iBin)],'color','k','LineWidth',2,'marker','o','markerfacecolor','none','markeredgecolor','k','markersize',10);
        hold on
    end
end
set(gca,'XTick',zeros(1,0),'xlim',[0.6 1.9+1.6*(numel(ClusBasedPermuTestIsSigBinClusSize)-1)+0.5]);
set(gca,'YTick',0:2:6,'YTickLabel',{'0','2','4','6'},'FontName','Arial','FontSize',16,'ylim',[0 6]);
ylabel('Persistence (s)','FontSize',18,'FontName','Arial');
box off;
title(['p_mPFC=' num2str(p(1)) '; p_aAIC=' num2str(p(2))]);
set(gcf,'Renderer','Painter'); saveas(gcf,'Persistence-sustained and transient neurons','fig'); close;



