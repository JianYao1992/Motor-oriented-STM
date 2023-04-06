%% Plot OGNG behavioral performance in learning phase for population of mice

clear; clc; close all;

%% Enter pathway
CurrPath = uigetdir;
cd(CurrPath);

%% Load results
% control group
[ C_Filename Pathway ] = uigetfile({'*.mat','Matlab files(*.mat)';},...
    'Pick some file','MultiSelect','on');
for iWindow = 1:size(C_Filename,2)
    C_Dataofmice{iWindow} = load(C_Filename{1,iWindow});
end
% experimental group
[ E_Filename Pathway ] = uigetfile({'*.mat','Matlab files(*.mat)';},...
    'Pick some file','MultiSelect','on');
for iWindow = 1:size(E_Filename,2)
    E_Dataofmice{iWindow} = load(E_Filename{1,iWindow});
end

%% Assignment
task = 'OGNG';
Manip = 'activate mPFC-aAIC';
C1 = [0 0 0];
if strcmp(Manip,'suppress mPFC-aAIC')
    C2 = [0 125 0]/255;
else
    C2 = [67 106 178]/255;
end
BlockNum = 4;

%% Hit rate and CR rate
HitRate_Ctrl = cell(1,BlockNum);
CRRate_Ctrl = cell(1,BlockNum);
HitRate_Exp = cell(1,BlockNum);
CRRate_Exp = cell(1,BlockNum);
% control group
for iWindow = 1:BlockNum % block
    for j = 1:size(C_Dataofmice,2) % mice
        if iWindow <= size(C_Dataofmice{1,j}.WindowHitRate{1,1},2) % task blocks
            HitRate_Ctrl{iWindow} = [HitRate_Ctrl{iWindow}; C_Dataofmice{1,j}.WindowHitRate{1,1}(iWindow)];
            CRRate_Ctrl{iWindow} = [CRRate_Ctrl{iWindow}; C_Dataofmice{1,j}.WindowCRRate{1,1}(iWindow)];
        end
    end
end
% experimental group
for iWindow = 1:BlockNum % block
    for j = 1:size(E_Dataofmice,2) % mice
        if iWindow <= size(E_Dataofmice{1,j}.WindowHitRate{1,1},2) % task blocks 
            HitRate_Exp{iWindow} = [HitRate_Exp{iWindow}; E_Dataofmice{1,j}.WindowHitRate{1,1}(iWindow)];
            CRRate_Exp{iWindow} = [CRRate_Exp{iWindow}; E_Dataofmice{1,j}.WindowCRRate{1,1}(iWindow)];
        end
    end
end

%% Plot hit rate and CR rate
figure('OuterPosition',[219 203 750 600]);
% control group
MeanHitRate_Ctrl = cellfun(@mean, HitRate_Ctrl); % hit rate
StdHitRate_Ctrl = cellfun(@std, HitRate_Ctrl);
[nrows,ncols] = cellfun(@size, HitRate_Ctrl);
errorbar(1:BlockNum,MeanHitRate_Ctrl,StdHitRate_Ctrl./sqrt(nrows),'color',C1,'linestyle','--','LineWidth',2,'marker','o','markerfacecolor',[1 1 1],'markeredgecolor',C1,'markersize',10);
hold on
MeanCRRate_Ctrl = cellfun(@mean, CRRate_Ctrl); % CR rate
StdCRRate_Ctrl = cellfun(@std, CRRate_Ctrl);
[nrows,ncols] = cellfun(@size, CRRate_Ctrl);
errorbar(1:BlockNum,MeanCRRate_Ctrl,StdCRRate_Ctrl./sqrt(nrows),'color',C1,'LineWidth',2,'marker','o','markerfacecolor',C1,'markeredgecolor',C1,'markersize',10);
hold on
% experimental group
MeanHitRate_Exp = cellfun(@mean, HitRate_Exp); % hit rate
StdHitRate_Exp = cellfun(@std, HitRate_Exp);
[nrows,ncols] = cellfun(@size, HitRate_Exp);
errorbar(1:BlockNum,MeanHitRate_Exp,StdHitRate_Exp./sqrt(nrows),'color',C2,'linestyle','--','LineWidth',2,'marker','o','markerfacecolor',[1 1 1],'markeredgecolor',C2,'markersize',10);
MeanCRRate_Exp = cellfun(@mean, CRRate_Exp); % CR rate
StdCRRate_Exp = cellfun(@std, CRRate_Exp);
[nrows,ncols] = cellfun(@size, CRRate_Exp);
errorbar(1:BlockNum,MeanCRRate_Exp,StdCRRate_Exp./sqrt(nrows),'color',C2,'LineWidth',2,'marker','o','markerfacecolor',C2,'markeredgecolor',C2,'markersize',10);
legend('boxoff');
set(gca,'XTick',0:1:BlockNum,'XTickLabel',cell2str(num2cell(0:1:BlockNum)),'FontName','Arial','FontSize',16,'xlim',[0 BlockNum+1]);
set(gca,'YTick',0:0.2:1,'YTickLabel',{'0','20','40','60','80','100'},'FontName','Arial','FontSize',16,'ylim',[0 1]);
xlabel('Learning day','FontSize',18,'FontName','Arial');
ylabel('Hit and CR rates (%)','FontSize',18,'FontName','Arial');
box off;

%% Tw-ANOVA-md
p_hit = GetDataMatrixForMixedRepeatedAnova(HitRate_Ctrl, HitRate_Exp, 'HitRate');
p_hit = p_hit{1,1};
p_cr = GetDataMatrixForMixedRepeatedAnova(CRRate_Ctrl, CRRate_Exp, 'CRRate');
p_cr = p_cr{1,1};
title(['p_hit=' num2str(p_hit) '; p_cr=' num2str(p_cr)]);
set(gcf, 'Renderer', 'Painter'); saveas(gcf,['1Hit and CR rates-' Manip '-Learning-' task],'fig'); close all;

%% Save groups of mice
save('TwoGroupsMiceIDandPValue','C_Filename','E_Filename');
