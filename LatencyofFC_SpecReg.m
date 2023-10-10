%% Significant lagging time of spike correlation in the specific region

clear; clc; close all;

%% Assignment
Group = 'CtrlGroup';
Target = 'aAIC-mPFC';
DelayPeriod = 3:6;
thresh = norminv(0.995); % Bonferroni correction for multiple comparison

%% Load result of FC
load(['FCofNeurons_PropertyIndividual_' Target '_' Group '.mat']);

%% Spike lag of SSCC plot
FcPairOfEachSecond = FcPairOfEachSecond(:,DelayPeriod);
lag_Go = [];
lag_NoGo = [];
for iSec = 1:numel(FcPairOfEachSecond) % second
    fprintf('Analyze FC of %dth bin of delay period\n',iSec);
    % hit trials
    for iPair = 1:length(FcPairOfEachSecond{iSec}.s1)
        tempNormBinCounts_hit = (FcPairOfEachSecond{iSec}.s1{iPair}.hists1-smooth(FcPairOfEachSecond{iSec}.s1{iPair}.shufs1))./std(FcPairOfEachSecond{iSec}.s1{iPair}.shufs1);
        tempNormBinCounts_hit = tempNormBinCounts_hit(46:55);
        tempNormBinCounts_hit(5:6) = 0;
        pos = find(tempNormBinCounts_hit>thresh);
        AI = FcPairOfEachSecond{iSec}.s1{iPair}.AIs1;
        if AI > 0
            pos(pos>=6) = [];
            pos = 5-max(pos);
        else
            pos(pos<=5) = [];
            pos = min(pos) - 6;
        end
        lag_Go = [lag_Go; pos];
    end
    % CR trials
    for iPair = 1:length(FcPairOfEachSecond{iSec}.s2)
        tempNormBinCounts_CR = (FcPairOfEachSecond{iSec}.s2{iPair}.hists2-smooth(FcPairOfEachSecond{iSec}.s2{iPair}.shufs2))./std(FcPairOfEachSecond{iSec}.s2{iPair}.shufs2);
        tempNormBinCounts_CR = tempNormBinCounts_CR(46:55);
        tempNormBinCounts_CR(5:6) = 0;
        pos = find(tempNormBinCounts_CR>thresh);
        AI = FcPairOfEachSecond{iSec}.s2{iPair}.AIs2;
        if AI > 0
            pos(pos>=6) = [];
            pos = 5-max(pos);
        else
            pos(pos<=5) = [];
            pos = min(pos) - 6;
        end
        lag_NoGo = [lag_NoGo; pos];
    end
end

%% Save result
save(sprintf('Spike lag of FC-%s-%s',Target,Group),'lag_Go','lag_NoGo','-v7.3');

%% Bar plot
lag = vertcat(lag_Go,lag_NoGo);
[Counts,Centers] = hist(lag,0:1:max(lag));
figure;
bar(Centers,Counts/sum(Counts),1,'FaceColor',[0 0 0],'FaceAlpha',1,'EdgeColor','none'); hold on;
set(gca,'XTick',-0.5:1:max(lag)+0.5,'XTickLabel',num2cell(0:2:10),'XLim',[-0.5 max(lag)+0.5]);
box off;
set(gcf,'Render','Painter'); saveas(gcf,sprintf('Spike lag of FC-%s-%s',Target,Group),'fig'); close all;
