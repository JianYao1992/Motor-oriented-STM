%% Plot shift-predictor subtracted cross correlogram (SSCC).

clear; clc; close all;

%% Assignment
Group = 'CtrlGroup';
Target = 'mPFC-aAIC';
Period = {'Baseline','Sample','1stSecondOfDelay','2ndSecondOfDelay','3rdSecondOfDelay','4thSecondOfDelay','Test'};

%% ID of mPFC and aAIC neurons
load(['UnitsSummary_' Group '.mat']);
ID_mPFC = []; 
ID_aAIC = [];
% mPFC
for iUnit = 1:size(UnitsSummary.mPFC,1) 
    ID_mPFC = [ID_mPFC UnitsSummary.mPFC{iUnit,1}(1)];
end
% aAIC
for iUnit = 1:size(UnitsSummary.aAIC,1) 
    ID_aAIC = [ID_aAIC UnitsSummary.aAIC{iUnit,1}(1)];
end

%% Load result of FC
load(['FCofNeurons_PropertyIndividual_' Target '_' Group '.mat']);

%% SSCC plot of FC pair at each second
AIData = [];
for iSecond = 1:length(FcPairOfEachSecond) % second
    fprintf('Analyze FC of %sth bin over the time course of the trial\n',num2str(iSecond-2));
    % hit trials
    for iPair = 1:length(FcPairOfEachSecond{iSecond}.s1)
        figure('OuterPosition',[219 303 420 534]);
        tempNormBinCounts_hit = (FcPairOfEachSecond{iSecond}.s1{iPair}.hists1-smooth(FcPairOfEachSecond{iSecond}.s1{iPair}.shufs1))./std(FcPairOfEachSecond{iSecond}.s1{iPair}.shufs1);
        plot(tempNormBinCounts_hit,'r'); hold on
        tempAIPeak_hit = max(tempNormBinCounts_hit([46 47 48 49 52 53 54 55]));
        AI = FcPairOfEachSecond{iSecond}.s1{iPair}.AIs1;
        if AI > 0
            AI = 'P';
        else
            AI = 'N';
        end
        if strcmp(FcPairOfEachSecond{iSecond}.s1{iPair}.reg_su1,'mPFC') % unit 1
            tempsu1_ID = find(ID_mPFC==FcPairOfEachSecond{iSecond}.s1{iPair}.su1_clusterid);
            su1_ID = strcat('mPFC#',num2str(tempsu1_ID));
        else
            tempsu1_ID = find(ID_aAIC==FcPairOfEachSecond{iSecond}.s1{iPair}.su1_clusterid);
            su1_ID = strcat('aAIC#',num2str(tempsu1_ID));
        end
        if strcmp(FcPairOfEachSecond{iSecond}.s1{iPair}.reg_su2,'mPFC') % unit 2
            tempsu2_ID = find(ID_mPFC==FcPairOfEachSecond{iSecond}.s1{iPair}.su2_clusterid);
            su2_ID = strcat('mPFC#',num2str(tempsu2_ID));
        else
            tempsu2_ID = find(ID_aAIC==FcPairOfEachSecond{iSecond}.s1{iPair}.su2_clusterid);
            su2_ID = strcat('aAIC#',num2str(tempsu2_ID));
        end
        AIData(end+1,:) = [1 iSecond-2 tempsu1_ID tempsu2_ID tempAIPeak_hit]; % 1 column:go no-go; 2 column: second; 3 column: su1; 4 column: su2; 5 column: AI peak
        plot([45.5 45.5],[-20 20],'k--'); hold on
        plot([50.5 50.5],[0 20],'k--'); hold on
        plot([55.5 55.5],[-20 20],'k--'); hold on
        box off;
        title([Period{iSecond} '_hit_' su1_ID '_' su2_ID '_AI' num2str(AI) '_' Group]);
        set(gca,'XAxisLocation','origin','XTick',[0.5 50.5 100.5],'XTickLabel',{'-100','0','100'},'FontName','Arial','FontSize',16,'xlim',[0.5 100.5]);
        ylim([-1.2*tempAIPeak_hit 1.2*tempAIPeak_hit]);
        set(gcf,'Renderer','Painter'); saveas(gcf,[Period{iSecond} '_hit_' su1_ID '_' su2_ID '_AI' AI '_' Group],'fig'); close;
    end
    % CR trials
    for iPair = 1:length(FcPairOfEachSecond{iSecond}.s2)
        figure('OuterPosition',[219 303 420 534]);
        tempNormBinCounts_CR = (FcPairOfEachSecond{iSecond}.s2{iPair}.hists2-smooth(FcPairOfEachSecond{iSecond}.s2{iPair}.shufs2))./std(FcPairOfEachSecond{iSecond}.s2{iPair}.shufs2);
        plot(tempNormBinCounts_CR,'r');
        tempAIPeak_CR = max(tempNormBinCounts_CR([46 47 48 49 52 53 54 55]));
        AI = FcPairOfEachSecond{iSecond}.s2{iPair}.AIs2;
        if AI > 0
            AI = 'P';
        else
            AI = 'N';
        end
        if strcmp(FcPairOfEachSecond{iSecond}.s2{iPair}.reg_su1,'mPFC') % unit 1
            tempsu1_ID = find(ID_mPFC==FcPairOfEachSecond{iSecond}.s2{iPair}.su1_clusterid);
            su1_ID = strcat('mPFC Unit#',num2str(tempsu1_ID));
        else
            tempsu1_ID = find(ID_aAIC==FcPairOfEachSecond{iSecond}.s2{iPair}.su1_clusterid);
            su1_ID = strcat('aAIC Unit#',num2str(tempsu1_ID));
        end
        if strcmp(FcPairOfEachSecond{iSecond}.s2{iPair}.reg_su2,'mPFC') % unit 2
            tempsu2_ID = find(ID_mPFC==FcPairOfEachSecond{iSecond}.s2{iPair}.su2_clusterid);
            su2_ID = strcat('mPFC Unit#',num2str(tempsu2_ID));
        else
            tempsu2_ID = find(ID_aAIC==FcPairOfEachSecond{iSecond}.s2{iPair}.su2_clusterid);
            su2_ID = strcat('aAIC Unit#',num2str(tempsu2_ID));
        end
        AIData(end+1,:) = [2 iSecond-2 tempsu1_ID tempsu2_ID tempAIPeak_CR]; 
        hold on
        plot([45.5 45.5],[-20 20],'k--');
        hold on
        plot([50.5 50.5],[0 20],'k--');
        hold on
        plot([55.5 55.5],[-20 20],'k--');
        box off;
        title([Period{iSecond} '_CR_' su1_ID '_' su2_ID '_AI' num2str(AI) '_' Group]);
        set(gca,'XAxisLocation','origin','XTick',[0.5 50.5 100.5],'XTickLabel',{'-100','0','100'},'FontName','Arial','FontSize',16,'xlim',[0.5 100.5]);
        ylim([-1.2*tempAIPeak_CR 1.2*tempAIPeak_CR]);
        set(gcf,'Renderer','Painter'); saveas(gcf,[Period{iSecond} '_CR_' su1_ID '_' su2_ID '_AI' AI '_' Group],'fig'); close;
    end
end
