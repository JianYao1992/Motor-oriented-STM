%% Neurons showing response to all trunk movement, nose movement, and pupil dilation.

clear; clc; close all;

%% Assignment
RegName = {'mPFC','aAIC'};

%% Neurons responsive to all movements
for iReg = 1:numel(RegName)
    load(sprintf('Units responseness to movement-%s.mat',RegName{iReg}));
    % Go-preferred units in analysis
    AnalysisGoPreferredUnitsID = intersect(AnalysisGoPreferredUnitsID_Trunk,AnalysisGoPreferredUnitsID_Nose);
    AnalysisGoPreferredUnitsID = intersect(AnalysisGoPreferredUnitsID,AnalysisGoPreferredUnitsID_Pupil);
    % NoGo-preferred units in analysis
    AnalysisNoGoPreferredUnitsID = intersect(AnalysisNoGoPreferredUnitsID_Trunk,AnalysisNoGoPreferredUnitsID_Nose);
    AnalysisNoGoPreferredUnitsID = intersect(AnalysisNoGoPreferredUnitsID,AnalysisNoGoPreferredUnitsID_Pupil);
    % Activated units in analysis
    ActivatedUnitsID = intersect(SortedActivatedUnitsID_Trunk,SortedActivatedUnitsID_Nose);
    ActivatedUnitsID = intersect(ActivatedUnitsID,SortedActivatedUnitsID_Pupil);
    ActivatedUnitsID_GoPrefUnits = intersect(ActivatedUnitsID,AnalysisGoPreferredUnitsID);
    % Inhibited units in analysis
    InhibitedUnitsID = intersect(SortedInhibitedUnitsID_Trunk,SortedInhibitedUnitsID_Nose);
    InhibitedUnitsID = intersect(InhibitedUnitsID,SortedInhibitedUnitsID_Pupil);
    InhibitedUnitsID_NoGoPrefUnits = intersect(InhibitedUnitsID,AnalysisNoGoPreferredUnitsID);
    % figure
    figure('position',[200 200 600 400]);
    y = [numel(ActivatedUnitsID_GoPrefUnits)/numel(AnalysisGoPreferredUnitsID) 1-numel(ActivatedUnitsID_GoPrefUnits)/numel(AnalysisGoPreferredUnitsID); numel(InhibitedUnitsID_NoGoPrefUnits)/numel(AnalysisNoGoPreferredUnitsID) 1-numel(InhibitedUnitsID_NoGoPrefUnits)/numel(AnalysisNoGoPreferredUnitsID)];
    bar(y,'stacked');
    set(gca,'XTick',zeros(1,0),'XLim',[0.5 2.5],'YTick',0:0.2:1,'YTickLabel',{'0','20','40','60','80','100'},'YLim',[0 1]);
    box off
    set(gcf,'Render','Painter'); saveas(gcf,sprintf('Proportion of memory neurons modulating FR for all movements in memory neurons-%s',RegName{iReg}),'fig'); close all;
end

%% Proportion of movement-modulated neurons
% mPFC
load('Units responseness to movement-mPFC.mat');
UnitsID_AllMovements_mPFC = intersect(AnaUnitsID_Trunk,AnaUnitsID_Nose);
UnitsID_AllMovements_mPFC = intersect(UnitsID_AllMovements_mPFC,AnaUnitsID_Pupil);
ActivatedUnitsID_mPFC = intersect(SortedActivatedUnitsID_Trunk,SortedActivatedUnitsID_Nose);
ActivatedUnitsID_mPFC = intersect(ActivatedUnitsID_mPFC,SortedActivatedUnitsID_Pupil);
InhibitedUnitsID_mPFC = intersect(SortedInhibitedUnitsID_Trunk,SortedInhibitedUnitsID_Nose);
InhibitedUnitsID_mPFC = intersect(InhibitedUnitsID_mPFC,SortedInhibitedUnitsID_Pupil);

UnitsNum_Torso_mPFC = numel(AnaUnitsID_Trunk); % torso movement
UnitsNum_TorsoActivated_mPFC = numel(SortedActivatedUnitsID_Trunk);
UnitsNum_TorsoInhibited_mPFC = numel(SortedInhibitedUnitsID_Trunk);
UnitsNum_Nose_mPFC = numel(AnaUnitsID_Nose); % nose movement
UnitsNum_NoseActivated_mPFC = numel(SortedActivatedUnitsID_Nose);
UnitsNum_NoseInhibited_mPFC = numel(SortedInhibitedUnitsID_Nose);
UnitsNum_Pupil_mPFC = numel(AnaUnitsID_Pupil); % pupil dilation
UnitsNum_PupilActivated_mPFC = numel(SortedActivatedUnitsID_Pupil);
UnitsNum_PupilInhibited_mPFC = numel(SortedInhibitedUnitsID_Pupil);
% aAIC
load('Units responseness to movement-aAIC.mat');
UnitsID_AllMovements_aAIC = intersect(AnaUnitsID_Trunk,AnaUnitsID_Nose);
UnitsID_AllMovements_aAIC = intersect(UnitsID_AllMovements_aAIC,AnaUnitsID_Pupil);
ActivatedUnitsID_aAIC = intersect(SortedActivatedUnitsID_Trunk,SortedActivatedUnitsID_Nose);
ActivatedUnitsID_aAIC = intersect(ActivatedUnitsID_aAIC,SortedActivatedUnitsID_Pupil);
InhibitedUnitsID_aAIC = intersect(SortedInhibitedUnitsID_Trunk,SortedInhibitedUnitsID_Nose);
InhibitedUnitsID_aAIC = intersect(InhibitedUnitsID_aAIC,SortedInhibitedUnitsID_Pupil);

UnitsNum_Torso_aAIC = numel(AnaUnitsID_Trunk); % torso movement
UnitsNum_TorsoActivated_aAIC = numel(SortedActivatedUnitsID_Trunk);
UnitsNum_TorsoInhibited_aAIC = numel(SortedInhibitedUnitsID_Trunk);
UnitsNum_Nose_aAIC = numel(AnaUnitsID_Nose); % nose movement
UnitsNum_NoseActivated_aAIC = numel(SortedActivatedUnitsID_Nose);
UnitsNum_NoseInhibited_aAIC = numel(SortedInhibitedUnitsID_Nose);
UnitsNum_Pupil_aAIC = numel(AnaUnitsID_Pupil); % pupil dilation
UnitsNum_PupilActivated_aAIC = numel(SortedActivatedUnitsID_Pupil);
UnitsNum_PupilInhibited_aAIC = numel(SortedInhibitedUnitsID_Pupil);
% figure //////torso movement//////
figure('position',[200 200 600 400]);
y = [1-UnitsNum_TorsoActivated_mPFC/UnitsNum_Torso_mPFC-UnitsNum_TorsoInhibited_mPFC/UnitsNum_Torso_mPFC UnitsNum_TorsoInhibited_mPFC/UnitsNum_Torso_mPFC UnitsNum_TorsoActivated_mPFC/UnitsNum_Torso_mPFC; 1-UnitsNum_TorsoActivated_aAIC/UnitsNum_Torso_aAIC-UnitsNum_TorsoInhibited_aAIC/UnitsNum_Torso_aAIC UnitsNum_TorsoInhibited_aAIC/UnitsNum_Torso_aAIC UnitsNum_TorsoActivated_aAIC/UnitsNum_Torso_aAIC];
bar(y,'stacked');
set(gca,'XTick',zeros(1,0),'XLim',[0.5 2.5],'YTick',0:0.2:1,'YTickLabel',{'0','20','40','60','80','100'},'YLim',[0 1]);
box off
set(gcf,'Render','Painter'); saveas(gcf,'Proportion of neurons modulating FR-torso movement','fig'); close all;
% figure //////nose movement//////
figure('position',[200 200 600 400]);
y = [1-UnitsNum_NoseActivated_mPFC/UnitsNum_Nose_mPFC-UnitsNum_NoseInhibited_mPFC/UnitsNum_Nose_mPFC UnitsNum_NoseInhibited_mPFC/UnitsNum_Nose_mPFC UnitsNum_NoseActivated_mPFC/UnitsNum_Nose_mPFC; 1-UnitsNum_NoseActivated_aAIC/UnitsNum_Nose_aAIC-UnitsNum_NoseInhibited_aAIC/UnitsNum_Nose_aAIC UnitsNum_NoseInhibited_aAIC/UnitsNum_Nose_aAIC UnitsNum_NoseActivated_aAIC/UnitsNum_Nose_aAIC];
bar(y,'stacked');
set(gca,'XTick',zeros(1,0),'XLim',[0.5 2.5],'YTick',0:0.2:1,'YTickLabel',{'0','20','40','60','80','100'},'YLim',[0 1]);
box off
set(gcf,'Render','Painter'); saveas(gcf,'Proportion of neurons modulating FR-nose movement','fig'); close all;
% figure //////pupil dilation//////
figure('position',[200 200 600 400]);
y = [1-UnitsNum_PupilActivated_mPFC/UnitsNum_Pupil_mPFC-UnitsNum_PupilInhibited_mPFC/UnitsNum_Pupil_mPFC UnitsNum_PupilInhibited_mPFC/UnitsNum_Pupil_mPFC UnitsNum_PupilActivated_mPFC/UnitsNum_Pupil_mPFC; 1-UnitsNum_PupilActivated_aAIC/UnitsNum_Pupil_aAIC-UnitsNum_PupilInhibited_aAIC/UnitsNum_Pupil_aAIC UnitsNum_PupilInhibited_aAIC/UnitsNum_Pupil_aAIC UnitsNum_PupilActivated_aAIC/UnitsNum_Pupil_aAIC];
bar(y,'stacked');
set(gca,'XTick',zeros(1,0),'XLim',[0.5 2.5],'YTick',0:0.2:1,'YTickLabel',{'0','20','40','60','80','100'},'YLim',[0 1]);
box off
set(gcf,'Render','Painter'); saveas(gcf,'Proportion of neurons modulating FR-Pupil dilation','fig'); close all;
% figure //////all movements//////
figure('position',[200 200 600 400]);
y = [1-numel(InhibitedUnitsID_mPFC)/numel(UnitsID_AllMovements_mPFC)-numel(ActivatedUnitsID_mPFC)/numel(UnitsID_AllMovements_mPFC) numel(InhibitedUnitsID_mPFC)/numel(UnitsID_AllMovements_mPFC) numel(ActivatedUnitsID_mPFC)/numel(UnitsID_AllMovements_mPFC); 1-numel(InhibitedUnitsID_aAIC)/numel(UnitsID_AllMovements_aAIC)-numel(ActivatedUnitsID_aAIC)/numel(UnitsID_AllMovements_aAIC) numel(InhibitedUnitsID_aAIC)/numel(UnitsID_AllMovements_aAIC) numel(ActivatedUnitsID_aAIC)/numel(UnitsID_AllMovements_aAIC)];
bar(y,'stacked');
set(gca,'XTick',zeros(1,0),'XLim',[0.5 2.5],'YTick',0:0.2:1,'YTickLabel',{'0','20','40','60','80','100'},'YLim',[0 1]);
box off
set(gcf,'Render','Painter'); saveas(gcf,'Proportion of neurons modulating FR-All movements','fig'); close all;

%% Venn diagram //////proportion of modulated neurons in memory neurons//////
clear; clc; close all;

RegName = {'mPFC','aAIC'};
Cvenn = {[1 0 0],[0 1 0],[0 0 1]};
VennLabel = {'1','2','3','12','13','23','123'};

for iReg = 1:numel(RegName)
    load(sprintf('Units responseness to movement-%s.mat',RegName{iReg}));
    % Go-preferred units in analysis
    AnalysisGoPreferredUnitsID = intersect(AnalysisGoPreferredUnitsID_Trunk,AnalysisGoPreferredUnitsID_Nose);
    AnalysisGoPreferredUnitsID = intersect(AnalysisGoPreferredUnitsID,AnalysisGoPreferredUnitsID_Pupil);
    % Go-preferred units which were activated by torso movement in analysis
    AnalysisGoPreferredUnitsID_TorsoActivated = intersect(AnalysisGoPreferredUnitsID,SortedActivatedUnitsID_Trunk);
    % Go-preferred units which were activated by nose movement in analysis
    AnalysisGoPreferredUnitsID_NoseActivated = intersect(AnalysisGoPreferredUnitsID,SortedActivatedUnitsID_Nose);
    % Go-preferred units which were activated by pupil dilation in analysis
    AnalysisGoPreferredUnitsID_PupilActivated = intersect(AnalysisGoPreferredUnitsID,SortedActivatedUnitsID_Pupil);
    % Proportion of Go-preferred neurons activated by both torso and nose
    % movement
    Prop_ActivatedinGoPrefUnits_TorsoNoseOverlap = numel(intersect(AnalysisGoPreferredUnitsID_TorsoActivated,AnalysisGoPreferredUnitsID_NoseActivated))/numel(AnalysisGoPreferredUnitsID);
    % Proportion of Go-preferred neurons activated by both torso movement and
    % pupil dilation
    Prop_ActivatedinGoPrefUnits_TorsoPupilOverlap = numel(intersect(AnalysisGoPreferredUnitsID_TorsoActivated,AnalysisGoPreferredUnitsID_PupilActivated))/numel(AnalysisGoPreferredUnitsID);
    % Proportion of Go-preferred neurons activated by both nose movement and
    % pupil dilation
    Prop_ActivatedinGoPrefUnits_NosePupilOverlap = numel(intersect(AnalysisGoPreferredUnitsID_NoseActivated,AnalysisGoPreferredUnitsID_PupilActivated))/numel(AnalysisGoPreferredUnitsID);
    % Proportion of Go-preferred neurons activated by all movements
    Prop_ActivatedinGoPrefUnits_TorsoNosePupilOverlap = numel(intersect(intersect(AnalysisGoPreferredUnitsID_TorsoActivated,AnalysisGoPreferredUnitsID_NoseActivated),AnalysisGoPreferredUnitsID_PupilActivated))/numel(AnalysisGoPreferredUnitsID);
    % plot
    figure, axis equal, axis off
    F = struct('Display', 'iter');
    [~,S] = venn([numel(AnalysisGoPreferredUnitsID_TorsoActivated) numel(AnalysisGoPreferredUnitsID_NoseActivated) numel(AnalysisGoPreferredUnitsID_PupilActivated)]./numel(AnalysisGoPreferredUnitsID),[Prop_ActivatedinGoPrefUnits_TorsoNoseOverlap Prop_ActivatedinGoPrefUnits_TorsoPupilOverlap Prop_ActivatedinGoPrefUnits_NosePupilOverlap Prop_ActivatedinGoPrefUnits_TorsoNosePupilOverlap],'ErrMinMode','ChowRodgers','FaceColor',Cvenn,...
        'FaceAlpha', 0.6,'EdgeColor',Cvenn,'EdgeAlpha',0.6);
    % label the area
    gnNo = horzcat([numel(AnalysisGoPreferredUnitsID_TorsoActivated) numel(AnalysisGoPreferredUnitsID_NoseActivated) numel(AnalysisGoPreferredUnitsID_PupilActivated)]./numel(AnalysisGoPreferredUnitsID),[Prop_ActivatedinGoPrefUnits_TorsoNoseOverlap Prop_ActivatedinGoPrefUnits_TorsoPupilOverlap Prop_ActivatedinGoPrefUnits_NosePupilOverlap Prop_ActivatedinGoPrefUnits_TorsoNosePupilOverlap]);
    for i = 1:length(gnNo)
        text(S.ZoneCentroid(i,1), S.ZoneCentroid(i,2), [VennLabel{i} '-' num2str(100*gnNo(i)) '%']);
    end
    set(gcf,'Renderer','Painter'); saveas(gcf,['Proportion of overlapped activated neurons in Go-preferred neurons-' RegName{iReg}],'fig'); close;
    
    
    
    
    % NoGo-preferred units in analysis
    AnalysisNoGoPreferredUnitsID = intersect(AnalysisNoGoPreferredUnitsID_Trunk,AnalysisNoGoPreferredUnitsID_Nose);
    AnalysisNoGoPreferredUnitsID = intersect(AnalysisNoGoPreferredUnitsID,AnalysisNoGoPreferredUnitsID_Pupil);
    % NoGo-preferred units which were inhibited by torso movement in analysis
    AnalysisNoGoPreferredUnitsID_TorsoInhibited = intersect(AnalysisNoGoPreferredUnitsID,SortedInhibitedUnitsID_Trunk);
    % NoGo-preferred units which were inhibited by nose movement in analysis
    AnalysisNoGoPreferredUnitsID_NoseInhibited = intersect(AnalysisNoGoPreferredUnitsID,SortedInhibitedUnitsID_Nose);
    % NoGo-preferred units which were inhibited by pupil dilation in analysis
    AnalysisNoGoPreferredUnitsID_PupilInhibited = intersect(AnalysisNoGoPreferredUnitsID,SortedInhibitedUnitsID_Pupil);
    % Proportion of NoGo-preferred neurons inhibited by both torso and nose
    % movement
    Prop_InhibitedinNoGoPrefUnits_TorsoNoseOverlap = numel(intersect(AnalysisNoGoPreferredUnitsID_TorsoInhibited,AnalysisNoGoPreferredUnitsID_NoseInhibited))/numel(AnalysisNoGoPreferredUnitsID);
    % Proportion of NoGo-preferred neurons inhibited by both torso movement and
    % pupil dilation
    Prop_InhibitedinNoGoPrefUnits_TorsoPupilOverlap = numel(intersect(AnalysisNoGoPreferredUnitsID_TorsoInhibited,AnalysisNoGoPreferredUnitsID_PupilInhibited))/numel(AnalysisNoGoPreferredUnitsID);
    % Proportion of NoGo-preferred neurons inhibited by both nose movement and
    % pupil dilation
    Prop_InhibitedinNoGoPrefUnits_NosePupilOverlap = numel(intersect(AnalysisNoGoPreferredUnitsID_NoseInhibited,AnalysisNoGoPreferredUnitsID_PupilInhibited))/numel(AnalysisNoGoPreferredUnitsID);
    % Proportion of NoGo-preferred neurons inhibited by all movements
    Prop_InhibitedinNoGoPrefUnits_TorsoNosePupilOverlap = numel(intersect(intersect(AnalysisNoGoPreferredUnitsID_TorsoInhibited,AnalysisNoGoPreferredUnitsID_NoseInhibited),AnalysisNoGoPreferredUnitsID_PupilInhibited))/numel(AnalysisNoGoPreferredUnitsID);
    % plot
    figure, axis equal, axis off
    F = struct('Display', 'iter');
    [~,S] = venn([numel(AnalysisNoGoPreferredUnitsID_TorsoInhibited) numel(AnalysisNoGoPreferredUnitsID_NoseInhibited) numel(AnalysisNoGoPreferredUnitsID_PupilInhibited)]./numel(AnalysisNoGoPreferredUnitsID),[Prop_InhibitedinNoGoPrefUnits_TorsoNoseOverlap Prop_InhibitedinNoGoPrefUnits_TorsoPupilOverlap Prop_InhibitedinNoGoPrefUnits_NosePupilOverlap Prop_InhibitedinNoGoPrefUnits_TorsoNosePupilOverlap],'ErrMinMode','ChowRodgers','FaceColor',Cvenn,...
        'FaceAlpha', 0.6,'EdgeColor',Cvenn,'EdgeAlpha',0.6);
    % label the area
    gnNo = horzcat([numel(AnalysisNoGoPreferredUnitsID_TorsoInhibited) numel(AnalysisNoGoPreferredUnitsID_NoseInhibited) numel(AnalysisNoGoPreferredUnitsID_PupilInhibited)]./numel(AnalysisNoGoPreferredUnitsID),[Prop_InhibitedinNoGoPrefUnits_TorsoNoseOverlap Prop_InhibitedinNoGoPrefUnits_TorsoPupilOverlap Prop_InhibitedinNoGoPrefUnits_NosePupilOverlap Prop_InhibitedinNoGoPrefUnits_TorsoNosePupilOverlap]);
    for i = 1:length(gnNo)
        text(S.ZoneCentroid(i,1), S.ZoneCentroid(i,2), [VennLabel{i} '-' num2str(100*gnNo(i)) '%']);
    end
    set(gcf,'Renderer','Painter'); saveas(gcf,['Proportion of overlapped inhibited neurons in NoGo-preferred neurons-' RegName{iReg}],'fig'); close;
    
    
end



