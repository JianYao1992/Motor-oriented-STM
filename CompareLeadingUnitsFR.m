%% Compare FR of leading neurons

clear; clc; close all;

%% Assignment
Reg = 'Within aAIC';

%% Load result of control group
load(sprintf('Decreased FU selectivity with removing FCSP-related spikes-%s-CtrlGroup.mat',Reg));

%% Minimal number of leading neurons
MinLuNum = size(FCpair_MemtoMem,2);
FR_MemtoMem = FR_MemtoMem(:,1:MinLuNum);
FR_NonmemtoMem = FR_NonmemtoMem(:,1:MinLuNum);

%% Bootstrap test
for iLuNum = 1:size(FR_MemtoMem,2)
    temp_MemtoMem = FR_MemtoMem(:,iLuNum);
    temp_MemtoMem = vertcat(temp_MemtoMem{:});
    PairNum_MemtoMem = numel(temp_MemtoMem);
    temp_NonmemtoMem = FR_NonmemtoMem(:,iLuNum);
    temp_NonmemtoMem = vertcat(temp_NonmemtoMem{:});
    PairNum_NonmemtoMem = numel(temp_NonmemtoMem);
    MinPairNum = min(horzcat(PairNum_MemtoMem,PairNum_NonmemtoMem));
    if MinPairNum >= 200
        BsSize = 100*(floor(MinPairNum/100)-1);
    elseif MinPairNum >= 100 && MinPairNum < 200
        BsSize = 50*(floor(MinPairNum/50)-1);
    elseif MinPairNum < 100 
        BsSize = 5*floor(MinPairNum/10);
    end
    BsTimes = 1000;
    BsValue = cell(1,2);
    for k = 1:BsTimes
        RandID = randperm(numel(temp_MemtoMem));
        RandID = RandID(1:BsSize);
        BsValue{1,1} = [BsValue{1,1}; mean(temp_MemtoMem(RandID))];
        RandID = randperm(numel(temp_NonmemtoMem));
        RandID = RandID(1:BsSize);
        BsValue{1,2} = [BsValue{1,2}; mean(temp_NonmemtoMem(RandID))];
    end
    bar(1.1+1.8*(iLuNum-1),mean(BsValue{1,1}),0.6,'facecolor',[0 0 0]); hold on
    CI_MemtoMem = prctile(BsValue{1,1},[2.5 97.5],1);
    errorbar(1.1+1.8*(iLuNum-1),mean(BsValue{1,1}),mean(BsValue{1,1})-CI_MemtoMem(1),CI_MemtoMem(2)-mean(BsValue{1,1}),'k','marker','none'); hold on
    bar(1.9+1.8*(iLuNum-1),mean(BsValue{1,2}),0.6,'facecolor',[0 0 0]); hold on
    CI_NonmemtoMem = prctile(BsValue{1,2},[2.5 97.5],1);
    errorbar(1.9+1.8*(iLuNum-1),mean(BsValue{1,2}),mean(BsValue{1,2})-CI_NonmemtoMem(1),CI_NonmemtoMem(2)-mean(BsValue{1,2}),'k','marker','none'); hold on
    diffvalue = BsValue{1,1} - BsValue{1,2};
    if mean(diffvalue) <= 0
        tempP = nnz(diffvalue>=0)/BsTimes;
    else
        tempP = nnz(diffvalue<=0)/BsTimes;
    end
    text(1.1+1.8*(iLuNum-1),mean(BsValue{1,1})/2,['n = ' num2str(numel(temp_MemtoMem))],'color','r'); hold on
    text(1.5+1.8*(iLuNum-1),mean(BsValue{1,1}),['p = ' num2str(tempP)],'color','r'); hold on
    text(1.9+1.8*(iLuNum-1),mean(BsValue{1,1})/2,['n = ' num2str(numel(temp_NonmemtoMem))],'color','r'); hold on
end
box off
title('Bootstrap test without Bonferroni correction');
set(gca,'XLim',[0.6 2.4+1.8*(size(FR_MemtoMem,2)-1)],'YTick',0:5:20,'YTickLabel',num2cell(0:5:20),'YLim',[0 20]);
set(gcf,'Render','Painter'); saveas(gcf,sprintf('Compare FR of LU between memory and nonmemory-%s-CtrlGroup',Reg),'fig'); close all;













