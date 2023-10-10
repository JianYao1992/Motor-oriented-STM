%% Compare proportion of memory following neurons in FC between mPFC-aAIC and aAIC-mPFC neurons

clear; clc; close all;

%% Load proportion data of mPFC-aAIC pairs
load('Proportion of following memory neurons in FC in Combine trials_mPFC-aAIC_CtrlGroup.mat');
PropMemUnit_mPFCtoaAIC = PropOfMemoryNeuron;

%% Load proportion data of aAIC-mPFC pairs
load('Proportion of following memory neurons in FC in Combine trials_aAIC-mPFC_CtrlGroup.mat');
PropMemUnit_aAICtomPFC = PropOfMemoryNeuron;

%% Bar plot
figure('Color','w','Position',[100,100,235,260])
hold on
bar(1.1,PropMemUnit_mPFCtoaAIC(1)/PropMemUnit_mPFCtoaAIC(2),0.6,'EdgeColor','k','FaceColor','k'); hold on
bar(1.9,PropMemUnit_aAICtomPFC(1)/PropMemUnit_aAICtomPFC(2),0.6,'EdgeColor','k','FaceColor','k'); hold on
bar(2.9,PropMemUnit_mPFCtoaAIC(3)/PropMemUnit_mPFCtoaAIC(4),0.6,'EdgeColor','k','FaceColor','w'); hold on
bar(3.7,PropMemUnit_aAICtomPFC(3)/PropMemUnit_aAICtomPFC(4),0.6,'EdgeColor','k','FaceColor','w'); hold on
[~,p_mem,~,~] = prop_test([PropMemUnit_mPFCtoaAIC(1) PropMemUnit_aAICtomPFC(1)],[PropMemUnit_mPFCtoaAIC(2) PropMemUnit_aAICtomPFC(2)],false);
[~,p_nonmem,~,~] = prop_test([PropMemUnit_mPFCtoaAIC(3) PropMemUnit_aAICtomPFC(3)],[PropMemUnit_mPFCtoaAIC(4) PropMemUnit_aAICtomPFC(4)],false);
[~,p_aAICmPFC,~,~] = prop_test([PropMemUnit_aAICtomPFC(1) PropMemUnit_aAICtomPFC(3)],[PropMemUnit_aAICtomPFC(2) PropMemUnit_aAICtomPFC(4)],false);
p_mem = 3*p_mem;
p_nonmem = 3*p_nonmem;
p_aAICmPFC = 3*p_aAICmPFC;
title(['p_mem = ' num2str(p_mem) '; p_nonmem = ' num2str(p_nonmem) '; p_aAICmPFC = ' num2str(p_aAICmPFC)]);
xlim([0.6 4.2])
set(gca,'YTick',0:0.2:1,'ylim',[0 1],'XTick',zeros(1,0));
ylabel('Proportion of memory F.U. (%)');
set(gcf,'Renderer','Painter'); saveas(gcf,'Proportion of following memory neuron in FC between mPFC-aAIC and aAIC-mPFC pairs','fig'); close;