%% Compare FC density of congruent active, congruent inactive, incongruent active, incongruent inactive and non-memory neuronal pairs,
%% between mPFC-aAIC and aAIC-mPFC
% /// low or high FR ///

clear; clc; close all;

%% Assignment
Group = 'CtrlGroup';
Type = [1 2 4]; % {'ConAct','ConInact','InconAct','InconInact','Non-memory'};
FRrangeID = 2; % 1: 2-5 Hz; 2: 5-12 Hz

%% Load result of FC density
load(sprintf('Region-specific FC density within FR 2 to 5 to 12-%s.mat',Group));
FcDensity = [{IsSigFC_hit_ConAct},{IsSigFC_hit_ConInact},{IsSigFC_hit_InconAct},{IsSigFC_hit_InconInact},{IsSigFC_hit_Nonmem}];

%% FC density of target type of pairs, controlled for FR
FcDensity = FcDensity(:,Type);
FcDensity_mPFCaAIC = cell(1,numel(Type));
FcDensity_aAICmPFC = cell(1,numel(Type));
for itype = 1:numel(Type)
    FcDensity_mPFCaAIC{itype} = cellfun(@mean,FcDensity{itype}{4,FRrangeID});
    FcDensity_mPFCaAIC{itype} = FcDensity_mPFCaAIC{itype}(:);
    FcDensity_aAICmPFC{itype} = cellfun(@mean,FcDensity{itype}{5,FRrangeID});
    FcDensity_aAICmPFC{itype} = FcDensity_aAICmPFC{itype}(:);
end
% Tw-ANOVA-md
Frame = ConstructTwANOVAmdTestMatrix(FcDensity_mPFCaAIC,FcDensity_aAICmPFC);
[SSQs, DFs, MSQs, Fs, Ps] = mixed_between_within_anova(Frame,0,'P value_FC density bewteen mPFC-aAIC and aAIC-mPFC');

%% Bar plot
figure('position',[200 200 600 400]);
AverFcDensity_mPFCaAIC = cellfun(@mean,FcDensity_mPFCaAIC);
SemFcDensity_mPFCaAIC = cellfun(@(x) std(x)/sqrt(size(x,1)),FcDensity_mPFCaAIC);
AverFcDensity_aAICmPFC = cellfun(@mean,FcDensity_aAICmPFC);
SemFcDensity_aAICmPFC = cellfun(@(x) std(x)/sqrt(size(x,1)),FcDensity_aAICmPFC);
for itype = 1:numel(Type)
    bar(1.1+1.8*(itype-1),AverFcDensity_mPFCaAIC(itype),0.6,'r'); hold on
    errorbar(1.1+1.8*(itype-1),AverFcDensity_mPFCaAIC(itype),SemFcDensity_mPFCaAIC(itype),'color','r'); hold on
    bar(1.9+1.8*(itype-1),AverFcDensity_aAICmPFC(itype),0.6,'k'); hold on
    errorbar(1.9+1.8*(itype-1),AverFcDensity_aAICmPFC(itype),SemFcDensity_aAICmPFC(itype),'color','k'); hold on
end
set(gca,'XTick',zeros(1,0),'xlim',[0.6 1.9+1.8*(numel(Type)-1)+0.5]);
set(gca,'YTick',0:0.1:1,'YTickLabel',num2cell(0:10:100),'FontName','Arial','FontSize',16,'ylim',[0 0.2]);
ylabel('FC density (%)','FontSize',18,'FontName','Arial');
box off;
set(gcf,'Renderer','Painter'); saveas(gcf,['CompareFCdensity-PairsType' num2str(Type) '-Between mPFC-aAIC and aAIC-mPFC'],'fig'); close;


