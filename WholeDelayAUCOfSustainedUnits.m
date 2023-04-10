%% auROC of delay-period FR of sustained neurons in mPFC and aAIC

clear; clc; close all;

%% Assignment 
DelayBin = 31:70;
SustainedUnitID = 0;
TargetID = [];
IsDVOrFRROC = 1;

%% Load results
load('UnitsSummary_CtrlGroup.mat'); % selectivity
UnitsSummary  = vertcat(UnitsSummary.mPFC,UnitsSummary.aAIC);
load('All Units FR and AUC_CtrlGroup.mat'); % ID of mPFC and aAIC neurons
load('mPFCFRinS1S2_CtrlGroup'); UnitsFRinS1S2_mPFC = vertcat(TargetBrainUnitsFRinS1,TargetBrainUnitsFRinS2); UnitsFRinS1S2_mPFC = UnitsFRinS1S2_mPFC'; % FR
load('aAICFRinS1S2_CtrlGroup'); UnitsFRinS1S2_aAIC = vertcat(TargetBrainUnitsFRinS1,TargetBrainUnitsFRinS2); UnitsFRinS1S2_aAIC = UnitsFRinS1S2_aAIC';
UnitsFRinS1S2 = vertcat(UnitsFRinS1S2_mPFC,UnitsFRinS1S2_aAIC);

%% Averaged delay-period FR for individual neuron
for iUnit = 1:size(UnitsSummary,1)
    if all(UnitsSummary{iUnit,3}(3:6) == 1) | all(UnitsSummary{iUnit,3}(3:6) == 2)
        SustainedUnitID = SustainedUnitID + 1;
        if strcmp(UnitsSummary{iUnit,6},'mPFC') RegionID = 1; UnitID = find(ID_mPFC == UnitsSummary{iUnit,1}(1)); % mPFC
        else RegionID = 2; UnitID = find(ID_aAIC == UnitsSummary{iUnit,1}(1)); % aAIC
        end
        for iTrialType = 1:2
            UnitsMeanDelayFR{SustainedUnitID,iTrialType} = mean(UnitsFRinS1S2{iUnit,iTrialType}(:,DelayBin),2);
        end
        if all(UnitsSummary{iUnit,3}(3:6) == 1) SelecIndex = 1;
        else SelecIndex = 2;
        end
        TargetID = [TargetID; horzcat(RegionID,UnitID,SelecIndex)];
    end
end

%% auROC of individual sustained neuron
for iUnit = 1:size(UnitsMeanDelayFR,1)
    if iUnit == 1 | mod(iUnit,10) == 0 | iUnit == size(UnitsMeanDelayFR,1)
        disp(['Current Unit ID = ' num2str(iUnit)]);
    end
    % step 1: delay-period FR in Go and NoGo trials
    FRofWholeDelay_Go = UnitsMeanDelayFR{iUnit,1};
    FRofWholeDelay_NoGo = UnitsMeanDelayFR{iUnit,2};
    if IsDVOrFRROC == 1
        % step 2: DV in Go and NoGo trials
        [DV_Go,DV_NoGo] = DecisionVariableCalculation(FRofWholeDelay_Go,FRofWholeDelay_NoGo);
        % step3: FPR and TPR
        [TPR,FPR] = TprFprCalculation(DV_Go,DV_NoGo);
    else
        % step3: FPR and TPR
        [TPR,FPR] = TprFprCalculation(FRofWholeDelay_Go',FRofWholeDelay_NoGo');
    end
    % step4: auROC
    TargetUnitsAUC(iUnit,1) = AucAnalysis(FPR,TPR);
end
TargetUnitsAUC = horzcat(TargetID,TargetUnitsAUC);

%% Histogram plot for distribution of sustained neurons with differernt auROC
MaxAUC = ceil(max(TargetUnitsAUC(:,4)));
MinAUC = 0.6;
figure;
[Counts,Centers] = hist(TargetUnitsAUC(:,4), MinAUC:(MaxAUC - MinAUC)/20:MaxAUC);
bar(Centers,Counts/length(TargetUnitsAUC(:,4)),'FaceColor',[1 0 0],'FaceAlpha',1,'EdgeColor','none'); hold on;
TargetUnitsAUC_NogoPreferred = TargetUnitsAUC(TargetUnitsAUC(:,3) == 2,4);
[Counts,Centers] = hist(TargetUnitsAUC_NogoPreferred, MinAUC:(MaxAUC - MinAUC)/20:MaxAUC);
bar(Centers,Counts/length(TargetUnitsAUC(:,4)),'FaceColor',[0 0 1],'FaceAlpha',1,'EdgeColor','none'); hold on;
MidValue  = median(TargetUnitsAUC(:,4));
plot([MidValue MidValue],[0 0.16],'k--');
xlim([0.65 1]); ylim([0 0.16]); box off;
set(gcf,'Renderer','Painter'); saveas(gcf,'Distribution of AUC value of whole delay period of sustained neurons','fig'); close;
    
SortedUnitsAUC = sortrows(TargetUnitsAUC,4);
save('SustainedNeuronsAUCValueofWholeDelay','SortedUnitsAUC','Centers');



