
%% Compare effect of removing FCSP on following neurons between control and optogenetic groups

clear; clc; close all;

%% Assignment
Reg = 'Within aAIC';
FRrange = [10 13];

%% Load result
% inhibition group
load(sprintf('Decreased FU AUC with removing FCSP-related spikes-FRrange-%d-%dHz-%s-NpHRGroup-1msbin.mat',FRrange(1),FRrange(2),Reg));
ChangedValue_MemtoMem_inhibition = ChangedValue_MemtoMem{1};
% ChangedValue_MemtoMem_inhibition(UniqFR_MemtoMem{1}>9 & UniqFR_MemtoMem{1}<10) = [];
% control group
load(sprintf('Decreased FU AUC with removing FCSP-related spikes-FRrange-%d-%dHz-%s-CtrlGroup-1msbin.mat',FRrange(1),FRrange(2),Reg));
ChangedValue_MemtoMem_CtrlGroup = ChangedValue_MemtoMem{1};
% ChangedValue_MemtoMem_CtrlGroup(UniqFR_MemtoMem{1}>9 & UniqFR_MemtoMem{1}<10) = [];
% activation group
load(sprintf('Decreased FU AUC with removing FCSP-related spikes-FRrange-%d-%dHz-%s-ChR2Group-1msbin.mat',FRrange(1),FRrange(2),Reg));
ChangedValue_MemtoMem_activation = ChangedValue_MemtoMem{1};
% ChangedValue_MemtoMem_activation(UniqFR_MemtoMem{1}>9 & UniqFR_MemtoMem{1}<10) = [];

BsData_Ctrl = [];
BsData_Inhibition = [];
BsData_Activation = [];
BStimes = 1000;
for itr = 1:BStimes
    % control group
    tempNum = numel(ChangedValue_MemtoMem_CtrlGroup);
    BsSize = tempNum-3;
    tempID = randperm(tempNum);
    tempData_Ctrl = ChangedValue_MemtoMem_CtrlGroup(tempID(1:BsSize));
    BsData_Ctrl = [BsData_Ctrl; mean(tempData_Ctrl)];
    % inhibition group
    tempNum = numel(ChangedValue_MemtoMem_inhibition);
    BsSize = tempNum-1;
    tempID = randperm(tempNum);
    tempData_inhibition = ChangedValue_MemtoMem_inhibition(tempID(1:BsSize));
    BsData_Inhibition = [BsData_Inhibition; mean(tempData_inhibition)];
    % activation group
    tempNum = numel(ChangedValue_MemtoMem_activation);
    BsSize = tempNum-5;
    tempID = randperm(tempNum);
    tempData_activation = ChangedValue_MemtoMem_activation(tempID(1:BsSize));
    BsData_Activation = [BsData_Activation; mean(tempData_activation)];
end
% bootstrap test between control and inhibition
tempDiff = BsData_Ctrl - BsData_Inhibition;
if mean(tempDiff) >= 0
    tempP_CtrlInhibition = nnz(tempDiff<=0)/BStimes;
else
    tempP_CtrlInhibition = nnz(tempDiff>=0)/BStimes;
end
% bootstrap test between control and activation
tempDiff = BsData_Ctrl - BsData_Activation;
if mean(tempDiff) >= 0
    tempP_CtrlActivation = nnz(tempDiff<=0)/BStimes;
else
    tempP_CtrlActivation = nnz(tempDiff>=0)/BStimes;
end
bar(1.1,mean(BsData_Inhibition),0.6,'facecolor',[0 125 0]/255,'edgecolor','none'); hold on
bar(1.9,mean(BsData_Ctrl),0.6,'facecolor',[0 0 0],'edgecolor','none'); hold on
bar(2.7,mean(BsData_Activation),0.6,'facecolor',[67 106 178]/255,'edgecolor','none'); hold on
text(1.5,mean(BsData_Ctrl)-0.05,['p=' num2str(tempP_CtrlInhibition)],'Color','red'); hold on
text(2.3,mean(BsData_Ctrl)-0.05,['p=' num2str(tempP_CtrlActivation)],'Color','red'); hold on

% confidence interval
tempCI_Inhibition = prctile(BsData_Inhibition,[2.5 97.5],1);
tempCI_Ctrl = prctile(BsData_Ctrl,[2.5 97.5],1);
tempCI_Activation = prctile(BsData_Activation,[2.5 97.5],1);
errorbar(1.1,mean(BsData_Inhibition),mean(BsData_Inhibition)-tempCI_Inhibition(1),tempCI_Inhibition(2)-mean(BsData_Inhibition),'color',[0 125 0]/255); hold on
errorbar(1.9,mean(BsData_Ctrl),mean(BsData_Ctrl)-tempCI_Ctrl(1),tempCI_Ctrl(2)-mean(BsData_Ctrl),'color',[0 0 0]); hold on
errorbar(2.7,mean(BsData_Activation),mean(BsData_Activation)-tempCI_Activation(1),tempCI_Activation(2)-mean(BsData_Activation),'color',[67 106 178]/255); hold on
title('Bootstrap of 1000 times without correction');
set(gca,'XTick',zeros(1,0),'xlim',[0.6 3.2]);
set(gca,'YTick',-0.15:0.05:0,'YTickLabel',num2cell(-15:5:0),'FontName','Arial','FontSize',16,'ylim',[-0.15 0]);
box on;
set(gcf,'Render','Painter'); saveas(gcf,sprintf('Compare decreased FU AUC with removing FCSP-MemtoMem-%s-FRrange-%d-%dHz-1msbin',Reg,FRrange(1),FRrange(2)),'fig'); close all;


