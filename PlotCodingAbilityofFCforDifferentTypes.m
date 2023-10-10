%% Plot FC coding ability during delay.

clear; clc; close all;

%% Assignment
Group = 'CtrlGroup';
TargetBrain = 'Within aAIC';
PairsType = {'CongruentActive','CongruentInactive','IncongruentActive','IncongruentInactive','Non-memory'};
load(sprintf('Functional Coupling Coding_%dpairs_%s_%s',length(PairsType),TargetBrain,Group));

%% Error bar plot
figure('OuterPosition',[219 150 420 350]);
t = 0.5:1:3.5;
C = [{[1 0 0]},{[0 0 1]},{[0 0 0]}];
TargetType = [1 2 4];
for i = 1:length(TargetType) % FC pairs types
    EmptyBinID = [];
    for bin = 1:length(TargetPairAUCofFC{TargetType(i)}) 
        if isempty(TargetPairAUCofFC{TargetType(i)}{bin})
            EmptyBinID = [EmptyBinID bin];
            data_Fcsp{1,i}{1,bin} = [];
            data_NonFcsp{1,i}{1,bin} = [];
        else
            data_Fcsp{1,i}{1,bin} = TargetPairAUCofFC{TargetType(i)}{bin}(:,3);
            data_NonFcsp{1,i}{1,bin} = TargetPairAUCofFC{TargetType(i)}{bin}(:,4);
        end
    end
    TargetPairAUCofFC{TargetType(i)}(:,EmptyBinID) = [];
    % FCSP
    tempResult = cellfun(@(x) x(:,3),TargetPairAUCofFC{TargetType(i)},'UniformOutput',false);
    tempAverage = cellfun(@mean,tempResult);
    tempSem = cellfun(@(x) std(x)/sqrt(size(x,1)),tempResult);
    tempT = t; 
    tempT(EmptyBinID) = [];
    errorbar(tempT,tempAverage,tempSem,'color',C{i},'linewidth',3,'marker','o','markerfacecolor',C{i},'markeredgecolor','none'); hold on;
    % non-FCSP
    tempResult = cellfun(@(x) x(:,4),TargetPairAUCofFC{TargetType(i)},'UniformOutput',false);
    tempAverage = cellfun(@mean,tempResult);
    tempSem = cellfun(@(x) std(x)/sqrt(size(x,1)),tempResult);
    tempT = t; 
    tempT(EmptyBinID) = [];
    errorbar(tempT,tempAverage,tempSem,'linestyle','--','color',C{i},'linewidth',3,'marker','o','markerfacecolor',C{i},'markeredgecolor','none'); hold on;
end
set(gca,'XTick',0:1:4,'XTickLabel',{'0','1','2','3','4'},'xlim',[0 4],'FontName','Arial','FontSize',16);
set(gca,'YTick',0.5:0.1:0.8,'YTickLabel',{'0.5','0.6','0.7','0.8'},'FontName','Arial','FontSize',16,'ylim',[0.47 0.8]);
xlabel('Time in delay (s)','FontSize',18,'FontName','Arial');
ylabel('Population averaged auROC','FontSize',18,'FontName','Arial');
box off;
plot([0 4],[0.5 0.5],'k--');
set(gcf,'Renderer','Painter'); saveas(gcf,['Functional Coupling Coding for 5 types in ' TargetBrain],'fig'); close;

%% Significance test
% between congruent active, congruent inactive, and incongruent inactive pairs for FCSP //////Two-way ANOVA//////
[~,tb1_BetwTypes_FCSP] = unbalanced_anova_test(data_Fcsp{1},data_Fcsp{2},data_Fcsp{3});
save(sprintf('P value_FC coding between Types_FCSP_%s_%s',TargetBrain,Group),'tb1_BetwTypes_FCSP','-v7.3');
% between FCSP and non-FCSP for congruent active pairs //////Tw-ANOVA-md//////
Frame = ConstructTwANOVAmdTestMatrix(data_Fcsp{1},data_NonFcsp{1});
[SSQs, DFs, MSQs, Fs, Ps] = mixed_between_within_anova(Frame,0,sprintf('P values_FC coding_between FCSP and NonFCSP_ConAct_%s_%s',TargetBrain,Group));
% between FCSP and non-FCSP for congruent inactive pairs //////Tw-ANOVA-md//////
Frame = ConstructTwANOVAmdTestMatrix(data_Fcsp{2},data_NonFcsp{2});
[SSQs, DFs, MSQs, Fs, Ps] = mixed_between_within_anova(Frame,0,sprintf('P values_FC coding_between FCSP and NonFCSP_ConInact_%s_%s',TargetBrain,Group));
% between FCSP and non-FCSP for incongruent inactive pairs //////Tw-ANOVA-md//////
Frame = ConstructTwANOVAmdTestMatrix(data_Fcsp{3},data_NonFcsp{3});
[SSQs, DFs, MSQs, Fs, Ps] = mixed_between_within_anova(Frame,0,sprintf('P values_FC coding_between FCSP and NonFCSP_InconInact_%s_%s',TargetBrain,Group));
% between congruent active, congruent inactive, and incongruent inactive pairs for non-FCSP //////Two-way ANOVA//////
[~,tb1_BetwTypes_nonFCSP] = unbalanced_anova_test(data_NonFcsp{1},data_NonFcsp{2},data_NonFcsp{3});
save(sprintf('P value_FC coding between Types_NonFCSP_%s_%s',TargetBrain,Group),'tb1_BetwTypes_nonFCSP','-v7.3');

