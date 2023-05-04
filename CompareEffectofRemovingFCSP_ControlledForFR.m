
%% Compare effect of removing FCSP on following neurons between control and optogenetic groups

clear; clc; close all;

%% Assignment
Reg = 'Within aAIC';
FRrange = [6 12];

%% Load result
% inhibition group
load(sprintf('Decreased FU selectivity with removing FCSP-related spikes-FRrange-%d-%dHz-%s-NpHRGroup.mat',FRrange(1),FRrange(2),Reg));
ChangedValue_MemtoMem_inhibition = ChangedValue_MemtoMem{1};
ChangedValue_NonmemtoMem_inhibition = ChangedValue_NonmemtoMem{1};
% control group
load(sprintf('Decreased FU selectivity with removing FCSP-related spikes-FRrange-%d-%dHz-%s-CtrlGroup.mat',FRrange(1),FRrange(2),Reg));
ChangedValue_MemtoMem_CtrlGroup = ChangedValue_MemtoMem{1};
ChangedValue_NonmemtoMem_CtrlGroup = ChangedValue_NonmemtoMem{1};
% activation group
load(sprintf('Decreased FU selectivity with removing FCSP-related spikes-FRrange-%d-%dHz-%s-ChR2Group.mat',FRrange(1),FRrange(2),Reg));
ChangedValue_MemtoMem_activation = ChangedValue_MemtoMem{1};
ChangedValue_NonmemtoMem_activation = ChangedValue_NonmemtoMem{1};

%% One-way ANOVA
% memory to memory
Data_MemtoMem = vertcat(ChangedValue_MemtoMem_inhibition,ChangedValue_MemtoMem_CtrlGroup,ChangedValue_MemtoMem_activation);
ID_MemtoMem = vertcat(ones(numel(ChangedValue_MemtoMem_inhibition),1),2*ones(numel(ChangedValue_MemtoMem_CtrlGroup),1),3*ones(numel(ChangedValue_MemtoMem_activation),1));
[~,tbl_MemtoMem] = anova1(Data_MemtoMem,ID_MemtoMem,'off');
% non-memory to memory
Data_NonmemtoMem = vertcat(ChangedValue_NonmemtoMem_inhibition,ChangedValue_NonmemtoMem_CtrlGroup,ChangedValue_NonmemtoMem_activation);
ID_NonmemtoMem = vertcat(ones(numel(ChangedValue_NonmemtoMem_inhibition),1),2*ones(numel(ChangedValue_NonmemtoMem_CtrlGroup),1),3*ones(numel(ChangedValue_NonmemtoMem_activation),1));
[~,tbl_NonmemtoMem] = anova1(Data_NonmemtoMem,ID_NonmemtoMem,'off');

%% Figure
% memory to memory
figure('position',[200 200 500 300]);
bar(1.1,mean(ChangedValue_MemtoMem_inhibition),0.6,'facecolor',[0 125 0]/255,'edgecolor','none'); hold on
errorbar(1.1,mean(ChangedValue_MemtoMem_inhibition),std(ChangedValue_MemtoMem_inhibition)/sqrt(numel(ChangedValue_MemtoMem_inhibition)),'color',[0 125 0]/255); hold on
bar(1.9,mean(ChangedValue_MemtoMem_CtrlGroup),0.6,'facecolor',[0 0 0]/255,'edgecolor','none'); hold on
errorbar(1.9,mean(ChangedValue_MemtoMem_CtrlGroup),std(ChangedValue_MemtoMem_CtrlGroup)/sqrt(numel(ChangedValue_MemtoMem_CtrlGroup)),'color',[0 0 0]/255); hold on
bar(2.7,mean(ChangedValue_MemtoMem_activation),0.6,'facecolor',[67 106 178]/255,'edgecolor','none'); hold on
errorbar(2.7,mean(ChangedValue_MemtoMem_activation),std(ChangedValue_MemtoMem_activation)/sqrt(numel(ChangedValue_MemtoMem_activation)),'color',[67 106 178]/255); hold on
xlim([0.6 3.2]);
set(gcf,'Render','Painter'); saveas(gcf,sprintf('Compare decreased FU selectivity with removing FCSP-MemtoMem-%s-FRrange-%d-%dHz',Reg,FRrange(1),FRrange(2)),'fig'); close all;
% non-memory to memory
figure('position',[200 200 500 300]);
bar(1.1,mean(ChangedValue_NonmemtoMem_inhibition),0.6,'facecolor',[0 125 0]/255,'edgecolor','none'); hold on
errorbar(1.1,mean(ChangedValue_NonmemtoMem_inhibition),std(ChangedValue_NonmemtoMem_inhibition)/sqrt(numel(ChangedValue_NonmemtoMem_inhibition)),'color',[0 125 0]/255); hold on
bar(1.9,mean(ChangedValue_NonmemtoMem_CtrlGroup),0.6,'facecolor',[0 0 0]/255,'edgecolor','none'); hold on
errorbar(1.9,mean(ChangedValue_NonmemtoMem_CtrlGroup),std(ChangedValue_NonmemtoMem_CtrlGroup)/sqrt(numel(ChangedValue_NonmemtoMem_CtrlGroup)),'color',[0 0 0]/255); hold on
bar(2.7,mean(ChangedValue_NonmemtoMem_activation),0.6,'facecolor',[67 106 178]/255,'edgecolor','none'); hold on
errorbar(2.7,mean(ChangedValue_NonmemtoMem_activation),std(ChangedValue_NonmemtoMem_activation)/sqrt(numel(ChangedValue_NonmemtoMem_activation)),'color',[67 106 178]/255); hold on
xlim([0.6 3.2]);
set(gcf,'Render','Painter'); saveas(gcf,sprintf('Compare decreased FU selectivity with removing FCSP-NonmemtoMem-%s-FRrange-%d-%dHz',Reg,FRrange(1),FRrange(2)),'fig'); close all;
save(sprintf('Compare decreased FU selectivity with removing FCSP-%s-FRrange-%d-%dHz',Reg,FRrange(1),FRrange(2)),'tbl_MemtoMem','tbl_NonmemtoMem','-v7.3');
