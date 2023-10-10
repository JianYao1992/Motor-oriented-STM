%% Compare FC density of congruent active, congruent inactive, incongruent active, incongruent inactive and non-memory neuronal pairs,
%% following projection-specific manipulation
% /// within mPFC, within aAIC, mPFC-aAIC, and aAIC-mPFC ///
% /// low or high FR ///

clear; clc; close all;

%% Assignment
Group = 'NpHRGroup';
if strcmp(Group,'NpHRGroup')
    C = [{[0 0 0]} {[0 125 0]/255}]; % black, green
elseif strcmp(Group,'ChR2Group')
    C = [{[0 0 0]} {[67 106 178]/255}]; % black, blue
end
Reg = 'mPFC-aAIC';
BStimes = 1000;
Type = [1 2 4]; % {'ConAct','ConInact','InconAct','InconInact','Non-memory'};
FRrangeID = 2; % 1: 2-5 Hz; 2: 5-12 Hz

%% Load result of FC density
% control group
load('Region-specific FC density within FR 2 to 5 to 12-CtrlGroup_1msbin.mat');
FcDensity_Ctrl = [{IsSigFC_hit_ConAct},{IsSigFC_hit_ConInact},{IsSigFC_hit_InconAct},{IsSigFC_hit_InconInact},{IsSigFC_hit_Nonmem}];
% experimental group
load(sprintf('Region-specific FC density within FR 2 to 5 to 12-%s_1msbin.mat',Group));
FcDensity_Opto = [{IsSigFC_hit_ConAct},{IsSigFC_hit_ConInact},{IsSigFC_hit_InconAct},{IsSigFC_hit_InconInact},{IsSigFC_hit_Nonmem}];

%% FC density of target type of pairs, controlled for FR
FcDensity_Ctrl = FcDensity_Ctrl(:,Type);
FcDensity_Opto = FcDensity_Opto(:,Type);
Data_Ctrl = cell(1,numel(Type));
Data_Opto = cell(1,numel(Type));
for itype = 1:numel(Type)
    if strcmp(Reg,'Within mPFC')
        tempData_Ctrl = FcDensity_Ctrl{itype}{2,FRrangeID};
        tempData_Opto = FcDensity_Opto{itype}{2,FRrangeID};
    elseif strcmp(Reg,'Within aAIC')
        tempData_Ctrl = FcDensity_Ctrl{itype}{3,FRrangeID};
        tempData_Opto = FcDensity_Opto{itype}{3,FRrangeID};
    elseif strcmp(Reg,'mPFC-aAIC')
        tempData_Ctrl = FcDensity_Ctrl{itype}{4,FRrangeID};
        tempData_Opto = FcDensity_Opto{itype}{4,FRrangeID};
    elseif strcmp(Reg,'aAIC-mPFC')
        tempData_Ctrl = FcDensity_Ctrl{itype}{5,FRrangeID};
        tempData_Opto = FcDensity_Opto{itype}{5,FRrangeID};
    end
    Data_Ctrl{itype} = vertcat(tempData_Ctrl{:});
    Data_Opto{itype} = vertcat(tempData_Opto{:});
end

%% Resampling and bar plot
figure('position',[200 200 600 400]);
BsData_Ctrl = cell(1,numel(Type));
BsData_Opto = cell(1,numel(Type));
for itype = 1:numel(Type)
    for itr = 1:BStimes
        % control group
        tempNum = numel(Data_Ctrl{itype});
        if tempNum <= 10
            BsSize = tempNum-1;
        elseif tempNum > 10 & tempNum <= 40
            BsSize = tempNum-2*(floor(tempNum/40)+1);
        else
            BsSize = tempNum-2*(floor(tempNum/40)+2);
        end
        tempID = randperm(tempNum);
        tempData_Ctrl = Data_Ctrl{itype}(tempID(1:BsSize));
        BsData_Ctrl{itype} = [BsData_Ctrl{itype}; mean(tempData_Ctrl)];
        % optogenetic group
        tempNum = numel(Data_Opto{itype});
        if tempNum <= 10
            BsSize = tempNum-1;
        elseif tempNum > 10 & tempNum <= 40
            BsSize = tempNum-2*(floor(tempNum/40)+1);
        else
            BsSize = tempNum-2*(floor(tempNum/40)+2);
        end
        tempID = randperm(tempNum);
        tempData_Opto = Data_Opto{itype}(tempID(1:BsSize));
        BsData_Opto{itype} = [BsData_Opto{itype}; mean(tempData_Opto)];
    end
    % bootstrap test
    tempDiff = BsData_Ctrl{itype} - BsData_Opto{itype};
    if mean(tempDiff) >= 0
        tempP = nnz(tempDiff<=0)/BStimes;
    else
        tempP = nnz(tempDiff>=0)/BStimes;
    end
    bar(1.1+1.8*(itype-1),mean(BsData_Ctrl{itype}),0.6,'facecolor',C{1},'edgecolor','none'); hold on
    bar(1.9+1.8*(itype-1),mean(BsData_Opto{itype}),0.6,'facecolor',C{2},'edgecolor','none'); hold on
    text(1.1+1.6*(itype-1),mean(BsData_Ctrl{itype})+0.05,['p=' num2str(tempP)],'Color','red');
    % confidence interval
    tempCI_Ctrl = prctile(BsData_Ctrl{itype},[2.5 97.5],1);
    tempCI_Opto = prctile(BsData_Opto{itype},[2.5 97.5],1);
    errorbar(1.1+1.8*(itype-1),mean(BsData_Ctrl{itype}),mean(BsData_Ctrl{itype})-tempCI_Ctrl(1),tempCI_Ctrl(2)-mean(BsData_Ctrl{itype}),'color',C{1}); hold on
    errorbar(1.9+1.8*(itype-1),mean(BsData_Opto{itype}),mean(BsData_Opto{itype})-tempCI_Opto(1),tempCI_Opto(2)-mean(BsData_Opto{itype}),'color',C{2}); hold on
end
set(gca,'XTick',zeros(1,0),'xlim',[0.6 1.9+1.8*(numel(Type)-1)+0.5]);
set(gca,'YTick',0:0.1:1,'YTickLabel',num2cell(0:10:100),'FontName','Arial','FontSize',16,'ylim',[0 0.50001]);
ylabel('FC density (%)','FontSize',18,'FontName','Arial');
box off;
set(gcf,'Renderer','Painter'); saveas(gcf,['CompareFCdensity-PairsType' num2str(Type) '-' Reg '-BetweenCtrland' Group '-1msbin-2to7to14Hz'],'fig'); close;


