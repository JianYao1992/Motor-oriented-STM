%% Compare FR of leading neurons, following projection-specific manipulation

clear; clc; close all;

%% Assignment
Reg = 'Within aAIC';
IsMemLU = 1;
if IsMemLU == 1
    PairType = 'MemtoMem';
else
    PairType = 'NonmemtoMem';
end

%% Load result of optogenetic inhibition group
load(sprintf('Decreased FU selectivity with removing FCSP-related spikes-%s-NpHRGroup.mat',Reg));
FR_MemtoMem_inhibition = FR_MemtoMem;
FR_NonmemtoMem_inhibition = FR_NonmemtoMem;

%% Load result of control group
load(sprintf('Decreased FU selectivity with removing FCSP-related spikes-%s-CtrlGroup.mat',Reg));
FR_MemtoMem_Ctrl = FR_MemtoMem;
FR_NonmemtoMem_Ctrl = FR_NonmemtoMem;

%% Load result of optogenetic activation group
load(sprintf('Decreased FU selectivity with removing FCSP-related spikes-%s-ChR2Group.mat',Reg));
FR_MemtoMem_activation = FR_MemtoMem;
FR_NonmemtoMem_activation = FR_NonmemtoMem;

%% Only pairs with one leading neuron
if IsMemLU == 1
    FR_inhibition = FR_MemtoMem_inhibition(:,1);
    FR_inhibition = vertcat(FR_inhibition{:});
    FR_Ctrl = FR_MemtoMem_Ctrl(:,1);
    FR_Ctrl = vertcat(FR_Ctrl{:});
    FR_activation = FR_MemtoMem_activation(:,1);
    FR_activation = vertcat(FR_activation{:});
else
    FR_inhibition = FR_NonmemtoMem_inhibition(:,1);
    FR_inhibition = vertcat(FR_inhibition{:});
    FR_Ctrl = FR_NonmemtoMem_Ctrl(:,1);
    FR_Ctrl = vertcat(FR_Ctrl{:});
    FR_activation = FR_NonmemtoMem_activation(:,1);
    FR_activation = vertcat(FR_activation{:});
end

%% One-way ANOVA
Data = vertcat(FR_inhibition,FR_Ctrl,FR_activation);
ID = vertcat(ones(numel(FR_inhibition),1),2*ones(numel(FR_Ctrl),1),3*ones(numel(FR_activation),1));
[~,tb1] = anova1(Data,ID,'off');

%% Error bar plot
bar(1.1,mean(FR_inhibition),0.6,'facecolor',[0 125 0]/255); hold on % NpHR group
errorbar(1.1,mean(FR_inhibition),std(FR_inhibition)/sqrt(numel(FR_inhibition)),'color',[0 125 0]/255,'marker','none'); hold on
bar(1.9,mean(FR_Ctrl),0.6,'facecolor',[0 0 0]/255); hold on % Ctrl group
errorbar(1.9,mean(FR_Ctrl),std(FR_Ctrl)/sqrt(numel(FR_Ctrl)),'k','marker','none'); hold on
bar(2.7,mean(FR_activation),0.6,'facecolor',[67 106 178]/255); hold on % ChR2 group
errorbar(2.7,mean(FR_activation),std(FR_activation)/sqrt(numel(FR_activation)),'color',[67 106 178]/255,'marker','none'); hold on
box off
set(gca,'XLim',[0.6 3.2],'YTick',0:10:30,'YTickLabel',num2cell(0:10:30),'YLim',[0 30]);
set(gcf,'Render','Painter'); saveas(gcf,sprintf('Compare FR of LU-%s-%s',PairType,Reg),'fig'); close all;

%% Save statistics table
save(sprintf('Compare FR of LU-%s-%s',PairType,Reg),'tb1','-v7.3');













