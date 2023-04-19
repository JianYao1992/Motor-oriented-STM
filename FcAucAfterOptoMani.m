%% auROC of FCSP events of control and optogenetic manipulation group

clear; clc; close all;

%% Assignment
Target = 'Within aAIC';
PairsType = {'Con Act','Con Inact','Incon Act','Incon Inact','Non-memory'}; 
PairsID = 1; 

%% auROC of FCSP events
% inhibition group
load(sprintf('Functional Coupling Coding_%dpairs_%s_NpHRGroup.mat',length(PairsType),Target));
data_inhibition = AllTypesAUC{PairsID};
PairsNuminEachBin_inhibition = cellfun(@(x) size(x,1),data_inhibition,'UniformOutput',1);
if all(PairsNuminEachBin_inhibition >= 2)
    FcspAucValue_inhibition = cellfun(@(x) x(:,3),data_inhibition,'UniformOutput',0);
end
% control group
load(sprintf('Functional Coupling Coding_%dpairs_%s_CtrlGroup.mat',length(PairsType),Target));
data_Ctrl = AllTypesAUC{PairsID};
PairsNuminEachBin_Ctrl = cellfun(@(x) size(x,1),data_Ctrl,'UniformOutput',1);
if all(PairsNuminEachBin_Ctrl >= 2)
    FcspAucValue_Ctrl = cellfun(@(x) x(:,3),data_Ctrl,'UniformOutput',0);
end
% activation group
load(sprintf('Functional Coupling Coding_%dpairs_%s_ChR2Group.mat',length(PairsType),Target));
data_activation = AllTypesAUC{PairsID};
PairsNuminEachBin_activation = cellfun(@(x) size(x,1),data_activation,'UniformOutput',1);
if all(PairsNuminEachBin_activation >= 2)
    FcspAucValue_activation = cellfun(@(x) x(:,3),data_activation,'UniformOutput',0);
end

%% Figures
if exist('FcspAucValue_inhibition') && exist('FcspAucValue_Ctrl') && exist('FcspAucValue_activation')
    figure('OuterPosition',[219 150 420 350]);
    % inhibition group
    meanAUC_inhibition = cellfun(@(x) mean(x),FcspAucValue_inhibition,'UniformOutput',true);
    stdAUC_inhibition = cellfun(@(x) std(x),FcspAucValue_inhibition,'UniformOutput',true);
    n_inhibition = cellfun(@(x) size(x,1),FcspAucValue_inhibition,'UniformOutput',true);
    % control group
    meanAUC_Ctrl = cellfun(@(x) mean(x),FcspAucValue_Ctrl,'UniformOutput',true);
    stdAUC_Ctrl = cellfun(@(x) std(x),FcspAucValue_Ctrl,'UniformOutput',true);
    n_Ctrl = cellfun(@(x) size(x,1),FcspAucValue_Ctrl,'UniformOutput',true);
    % activation group
    meanAUC_activation = cellfun(@(x) mean(x),FcspAucValue_activation,'UniformOutput',true);
    stdAUC_activation = cellfun(@(x) std(x),FcspAucValue_activation,'UniformOutput',true);
    n_activation = cellfun(@(x) size(x,1),FcspAucValue_activation,'UniformOutput',true);
    errorbar(1.5:1:4.5,meanAUC_inhibition,stdAUC_inhibition./sqrt(n_inhibition),'color',[0 125 0]/255,'linewidth',3,'marker','o','markerfacecolor',[0 125 0]/255,'markeredgecolor','none'); hold on
    errorbar(1.5:1:4.5,meanAUC_Ctrl,stdAUC_Ctrl./sqrt(n_Ctrl),'color','k','linewidth',3,'marker','o','markerfacecolor','k','markeredgecolor','none'); hold on
    errorbar(1.5:1:4.5,meanAUC_activation,stdAUC_activation./sqrt(n_activation),'color',[67 106 178]/255,'linewidth',3,'marker','o','markerfacecolor',[67 106 178]/255,'markeredgecolor','none'); hold on
    % unbalanced ANOVA test
    [p,tb1] = unbalanced_anova_test(FcspAucValue_inhibition,FcspAucValue_Ctrl,FcspAucValue_activation);
    title(['F(' num2str(tb1{2,3}) ')=' num2str(tb1{2,6}) 'P(1) = ' num2str(tb1{2,7}) 'F(' num2str(tb1{3,3}) ')=' num2str(tb1{3,6}) 'P(2) = ' num2str(tb1{3,7})]);
    set(gca,'XTick',0:2:6,'YTickLabel',{'0','2','4','6'},'xlim',[1 5],'FontName','Arial','FontSize',16);
    set(gca,'YTick',0.5:0.1:1,'YTickLabel',{'0.5','0.6','0.7','0.8','0.9','1'},'FontName','Arial','FontSize',16,'ylim',[0.45 0.8]);
    xlabel('Time from sample onset (s)','FontSize',18,'FontName','Arial');
    ylabel('Averaged delay auROC','FontSize',18,'FontName','Arial');
    box off;
    set(gcf,'Renderer','Painter'); saveas(gcf,['CompareFunctionalCouplingCoding_' Target '_' PairsType{PairsID} '_BetweenControlandOpto'],'fig'); close;
else
    disp('Not enough pairs');
end
