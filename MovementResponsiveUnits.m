%% ID of movement responsive (excited and suppressed) neurons in Go trials.
% //////mPFC and aAIC respectively//////

clear; clc; close all;

%% Assignment
RegName = {'mPFC','aAIC'};
BaseLen = 9;
SampOdorLen = 1;
DelayLen = 4;
TimeGain = 10;
MovementWindowNumber = 10;
BootsStrapNum = 1000;

%% Target directory
CurrPath = uigetdir;
AllPath = genpath(CurrPath);
SplitedPath = strsplit(AllPath,';');
SubPath = SplitedPath';
SubPath = SubPath(2:end-1);

%% Remove first-class pathway
FirstGradeFolder = dir(CurrPath);
ToRemoveFileID = [];
for i = 3:size(FirstGradeFolder,1)
    % pathway name to remove
    ToDeleteFolderName = strcat(CurrPath,'\',FirstGradeFolder(i).name);
    % pathway ID to remove
    for iFile = 1:length(SubPath)
        if strcmp(ToDeleteFolderName,SubPath(iFile)) == 1
            ToRemoveFileID = [ToRemoveFileID; iFile];
        end
    end
end
% remove pathway unrelated to data
SubPath(ToRemoveFileID)= [];

%% Align movement and neural correlates for all neurons
BM_Trunk = [];
BM_Nose = [];
BM_Pupil = [];
BMsig_Trunk = [];
BMsig_Nose = [];
BMsig_Pupil = [];
TrunkMoveWindow_delay = [];
NoseMoveWindow_delay = [];
PupilDilationWindow_delay = [];
UnitsCompletedTrialID = [];
UnitsTrialMark = [];
for iPath = 1:size(SubPath,1)
    Path = SubPath{iPath,1};
    cd(Path); % enter subpath
    DataFiles = dir('*SumExtraCellularAndVideo*.mat');
    for iFile = 1:size(DataFiles,1)
        % load extracellular data
        if size(DataFiles,1) == 1
            filename = DataFiles.name;
        else
            filename = DataFiles(iFile,1).name;
        end
        load(filename);
        BM_Trunk = [BM_Trunk repmat({TrialBasedBM_Trunk},1,size(SingleUnitList,1))];
        BMsig_Trunk = [BMsig_Trunk repmat({IsSigBM_Trunk},1,size(SingleUnitList,1))];
        TrunkMoveWindow_delay = [TrunkMoveWindow_delay repmat({IsTrunkMoveWindow_delay_Go},1,size(SingleUnitList,1))];
        BM_Nose = [BM_Nose repmat({TrialBasedBM_Nose},1,size(SingleUnitList,1))];
        BMsig_Nose = [BMsig_Nose repmat({IsSigBM_Nose},1,size(SingleUnitList,1))];
        NoseMoveWindow_delay = [NoseMoveWindow_delay repmat({IsNoseMoveWindow_delay_Go},1,size(SingleUnitList,1))];
        BM_Pupil = [BM_Pupil repmat({TrialBasedBM_Pupil},1,size(SingleUnitList,1))];
        BMsig_Pupil = [BMsig_Pupil repmat({IsSigBM_Pupil},1,size(SingleUnitList,1))];
        PupilDilationWindow_delay = [PupilDilationWindow_delay repmat({IsPupilDilationWindow_delay_Go},1,size(SingleUnitList,1))];
        UnitsCompletedTrialID = [UnitsCompletedTrialID repmat({CompletedTrialID},1,size(SingleUnitList,1))];
        UnitsTrialMark = [UnitsTrialMark repmat({TrialMark},1,size(SingleUnitList,1))];
    end
end

%% Movement responsive neurons of the target region
cd(CurrPath);
for iReg = 1:numel(RegName)
    AnaUnitsID_Trunk = [];
    AnaUnitsID_Nose = [];
    AnaUnitsID_Pupil = [];
    TarReg = RegName{iReg};
    load(sprintf('%sSelectivityData_CtrlGroup.mat',TarReg));
    FinishedTrialID = UnitsCompletedTrialID(:,TarUnitsID);
    TrialMark = UnitsTrialMark(:,TarUnitsID);
    % movement result of target neurons
    TrunkMovement = BM_Trunk(:,TarUnitsID);
    TrunkMovementSig = BMsig_Trunk(:,TarUnitsID);
    TrunkMoveWindow = TrunkMoveWindow_delay(:,TarUnitsID);
    NoseMovement = BM_Nose(:,TarUnitsID);
    NoseMovementSig = BMsig_Nose(:,TarUnitsID);
    NoseMoveWindow = NoseMoveWindow_delay(:,TarUnitsID);
    PupilDilation = BM_Pupil(:,TarUnitsID);
    PupilDilationSig = BMsig_Pupil(:,TarUnitsID);
    PupilDilationWindow = PupilDilationWindow_delay(:,TarUnitsID);
    % FR in delay period
    TargetBrainUnitsFRinS1 = cellfun(@(x) x(:,2+TimeGain*(BaseLen+SampOdorLen):1+TimeGain*(BaseLen+SampOdorLen+DelayLen)),TargetBrainUnitsFRinS1,'UniformOutput',0);
    % sum FR across delay bins
    MergedWindowNumber = (DelayLen/MovementWindowNumber)/(1/TimeGain);
    WindowMergedFRinS1 = cell(size(TargetBrainUnitsFRinS1));
    for iUnit = 1:numel(TargetBrainUnitsFRinS1)
        for iStepNum = 1:MovementWindowNumber
            WindowMergedFRinS1{iUnit}(:,iStepNum) = mean(TargetBrainUnitsFRinS1{iUnit}(:,1+MergedWindowNumber*(iStepNum-1):MergedWindowNumber*iStepNum),2);
        end
        if ~isempty(TrunkMoveWindow{iUnit})
            AnaUnitsID_Trunk = [AnaUnitsID_Trunk iUnit];
        end
        if ~isempty(NoseMoveWindow{iUnit})
            AnaUnitsID_Nose = [AnaUnitsID_Nose iUnit];
        end
        if ~isempty(PupilDilationWindow{iUnit})
            AnaUnitsID_Pupil = [AnaUnitsID_Pupil iUnit];
        end
    end
    % whether single unit showed modulation by movement
    RespID_Trunk = [];
    RespID_Nose = [];
    RespID_Pupil = [];
    AverFR_Trunk = cell(0);
    AverFR_Nose = cell(0);
    AverFR_Pupil = cell(0);
    AverFRdiff_Trunk = [];
    AverFRdiff_Nose = [];
    AverFRdiff_Pupil = [];
    for iUnit = 1:numel(WindowMergedFRinS1)
        tempFR = WindowMergedFRinS1{iUnit};
        if ~isempty(TrunkMoveWindow{iUnit})
            [tempRespID_Trunk,tempAverFR_Trunk,tempAverFRdiff_Trunk] = JudgeSingleUnitResponsiveToMovement(tempFR,TrunkMoveWindow{iUnit},BootsStrapNum);
            RespID_Trunk = [RespID_Trunk; horzcat(iUnit,tempRespID_Trunk)];
            AverFR_Trunk{end+1,1} = tempAverFR_Trunk;
            AverFRdiff_Trunk(end+1,1) = tempAverFRdiff_Trunk;
        end
        if ~isempty(NoseMoveWindow{iUnit})
            [tempRespID_Nose,tempAverFR_Nose,tempAverFRdiff_Nose] = JudgeSingleUnitResponsiveToMovement(tempFR,NoseMoveWindow{iUnit},BootsStrapNum);
            RespID_Nose = [RespID_Nose; horzcat(iUnit,tempRespID_Nose)];
            AverFR_Nose{end+1,1} = tempAverFR_Nose;
            AverFRdiff_Nose(end+1,1) = tempAverFRdiff_Nose;
        end
        if ~isempty(PupilDilationWindow{iUnit})
            [tempRespID_Pupil,tempAverFR_Pupil,tempAverFRdiff_Pupil] = JudgeSingleUnitResponsiveToMovement(tempFR,PupilDilationWindow{iUnit},BootsStrapNum);
            RespID_Pupil = [RespID_Pupil; horzcat(iUnit,tempRespID_Pupil)];
            AverFR_Pupil{end+1,1} = tempAverFR_Pupil;
            AverFRdiff_Pupil(end+1,1) = tempAverFRdiff_Pupil;
        end
    end
    % Go and NoGo-preferred neurons of target region
    load('Go-NoGo-PreferredUnitsID-CtrlGroup.mat');
    if strcmp(TarReg,'mPFC')
        GoPrefUnitsID = PreferCodingNeuronID.mPFC.GoPrefUnitsID;
        NoGoPrefUnitsID = PreferCodingNeuronID.mPFC.NoGoPrefUnitsID;
    elseif strcmp(TarReg,'aAIC')
        GoPrefUnitsID = PreferCodingNeuronID.aAIC.GoPrefUnitsID;
        NoGoPrefUnitsID = PreferCodingNeuronID.aAIC.NoGoPrefUnitsID;
    end
    
    %% Proportion of movement-activated neurons in Go-preferred neurons, and movement-inhibited neurons in NoGo-preferred neurons
    % trunk movement
    [SortedActivatedUnitsID_Trunk,SortedActivatedUnitsFR_Trunk,SortedInhibitedUnitsID_Trunk,SortedInhibitedUnitsFR_Trunk,AnalysisGoPreferredUnitsID_Trunk,AnalysisNoGoPreferredUnitsID_Trunk] = OverlappedMovmentResponsibleUnitsWithMemoryUnits(TarReg,'Trunk movement',RespID_Trunk,AverFR_Trunk,AverFRdiff_Trunk,GoPrefUnitsID,NoGoPrefUnitsID,AnaUnitsID_Trunk);
    % nose movement
    [SortedActivatedUnitsID_Nose,SortedActivatedUnitsFR_Nose,SortedInhibitedUnitsID_Nose,SortedInhibitedUnitsFR_Nose,AnalysisGoPreferredUnitsID_Nose,AnalysisNoGoPreferredUnitsID_Nose] = OverlappedMovmentResponsibleUnitsWithMemoryUnits(TarReg,'Nose movement',RespID_Nose,AverFR_Nose,AverFRdiff_Nose,GoPrefUnitsID,NoGoPrefUnitsID,AnaUnitsID_Nose);
    % pupil dilation
    [SortedActivatedUnitsID_Pupil,SortedActivatedUnitsFR_Pupil,SortedInhibitedUnitsID_Pupil,SortedInhibitedUnitsFR_Pupil,AnalysisGoPreferredUnitsID_Pupil,AnalysisNoGoPreferredUnitsID_Pupil] = OverlappedMovmentResponsibleUnitsWithMemoryUnits(TarReg,'Pupil dilation',RespID_Pupil,AverFR_Pupil,AverFRdiff_Pupil,GoPrefUnitsID,NoGoPrefUnitsID,AnaUnitsID_Pupil);
    
    %% Save ID of movement-responsible neurons
    save(sprintf('Units responseness to movement-%s',TarReg),'FinishedTrialID','TrialMark','TrunkMovement','NoseMovement','PupilDilation','TrunkMovementSig','NoseMovementSig','PupilDilationSig','SortedActivatedUnitsID_Trunk','SortedActivatedUnitsFR_Trunk','SortedInhibitedUnitsID_Trunk','SortedInhibitedUnitsFR_Trunk','AnalysisGoPreferredUnitsID_Trunk','AnalysisNoGoPreferredUnitsID_Trunk',...
        'SortedActivatedUnitsID_Nose','SortedActivatedUnitsFR_Nose','SortedInhibitedUnitsID_Nose','SortedInhibitedUnitsFR_Nose','AnalysisGoPreferredUnitsID_Nose','AnalysisNoGoPreferredUnitsID_Nose',...
        'SortedActivatedUnitsID_Pupil','SortedActivatedUnitsFR_Pupil','SortedInhibitedUnitsID_Pupil','SortedInhibitedUnitsFR_Pupil','AnalysisGoPreferredUnitsID_Pupil','AnalysisNoGoPreferredUnitsID_Pupil','AnaUnitsID_Trunk','AnaUnitsID_Nose','AnaUnitsID_Pupil','-v7.3');
end



function [UnitResponsenessID,AverFR,AverFRdiff] = JudgeSingleUnitResponsiveToMovement(FR,MovementBin,BootsStrapNum)

UnitResponsenessID = 0;
AverFR_Movement = [];
AverFR_Still = [];
for iBin = 1:size(FR,2)
    tempFR = FR(:,iBin);
    tempTrialID_Movement = find(MovementBin(:,iBin)==1);
    tempTrialID_Still = find(MovementBin(:,iBin)==0);
    if numel(tempTrialID_Movement) >= 1 && numel(tempTrialID_Still) >= 1
        CurrBinFR_Movement = tempFR(tempTrialID_Movement);
        CurrBinFR_Still = tempFR(tempTrialID_Still);
        if abs(size(CurrBinFR_Movement,1)-size(CurrBinFR_Still,1)) < 10
            RealBootsStrapNum = 1;
        else
            RealBootsStrapNum = BootsStrapNum;
        end
        MinTriNum = min([size(CurrBinFR_Movement,1) size(CurrBinFR_Still,1)]);
        
        %% Averaged FR in movement trials
        CurrBinAverFR_Movement = zeros(RealBootsStrapNum,1);
        for iBoot = 1:RealBootsStrapNum
            RandPickedTrialNum = randperm(size(CurrBinFR_Movement,1));
            RandomPickedTrialID = RandPickedTrialNum(1:MinTriNum);
            CurrBinAverFR_Movement(iBoot,:) = mean(CurrBinFR_Movement(RandomPickedTrialID,:));
        end
        if RealBootsStrapNum > 1
            CurrBinAverFR_Movement = mean(CurrBinAverFR_Movement);
        end
        
        %% Averaged FR in still trials
        CurrBinAverFR_Still = zeros(RealBootsStrapNum,1);
        for iBoot = 1:RealBootsStrapNum
            RandPickedTrialNum = randperm(size(CurrBinFR_Still,1));
            RandomPickedTrialID = RandPickedTrialNum(1:MinTriNum);
            CurrBinAverFR_Still(iBoot,:) = mean(CurrBinFR_Still(RandomPickedTrialID,:));
        end
        if RealBootsStrapNum > 1
            CurrBinAverFR_Still = mean(CurrBinAverFR_Still);
        end
        AverFR_Movement(end+1,1) = CurrBinAverFR_Movement;
        AverFR_Still(end+1,1) = CurrBinAverFR_Still;
    end
end
AverFR = horzcat(AverFR_Still,AverFR_Movement);
AverFRdiff = mean(AverFR_Movement) - mean(AverFR_Still);
% signrank test
p = signrank(AverFR_Movement,AverFR_Still);
if p <= 0.05
    if AverFRdiff >= 1
        UnitResponsenessID = 1; % activated
    elseif AverFRdiff <= -1
        UnitResponsenessID = 2; % inhibited
    end
end
end

function [SortedActivatedUnitsID,SortedActivatedUnitsFR,SortedInhibitedUnitsID,SortedInhibitedUnitsFR,AnalysisGoPreferredUnitsID,AnalysisNoGoPreferredUnitsID] = OverlappedMovmentResponsibleUnitsWithMemoryUnits(RegName,ParameterName,ResponsiveUnitsID,FR,FRdiffBetwMoveAndStill,GoPreferredUnitsID,NoGoPreferredUnitsID,AnalysisUnitsID)
AnalysisGoPreferredUnitsID = intersect(GoPreferredUnitsID,AnalysisUnitsID);
AnalysisNoGoPreferredUnitsID = intersect(NoGoPreferredUnitsID,AnalysisUnitsID);
ActivatedUnitsID = ResponsiveUnitsID(ResponsiveUnitsID(:,2)==1,1);
InhibitedUnitsID = ResponsiveUnitsID(ResponsiveUnitsID(:,2)==2,1);
ActivatedUnitsID_GoPref = intersect(ActivatedUnitsID,AnalysisGoPreferredUnitsID);
InhibitedUnitsID_NoGoPref = intersect(InhibitedUnitsID,AnalysisNoGoPreferredUnitsID);
PropActUnits_GoPref_Torso = numel(ActivatedUnitsID_GoPref) / numel(AnalysisGoPreferredUnitsID);
PropNonActUnits_GoPref_Torso = 1 - PropActUnits_GoPref_Torso;
PropInhUnits_NoGoPref_Torso = numel(InhibitedUnitsID_NoGoPref) / numel(AnalysisNoGoPreferredUnitsID);
PropNonInhUnits_NoGoPref_Torso = 1 - PropInhUnits_NoGoPref_Torso;
figure('position',[200 200 600 400]);
y = [PropActUnits_GoPref_Torso PropNonActUnits_GoPref_Torso; PropInhUnits_NoGoPref_Torso PropNonInhUnits_NoGoPref_Torso];
bar(y,'stacked');
set(gca,'XTick',zeros(1,0),'XLim',[0.5 2.5],'YTick',0:0.2:1,'YTickLabel',{'0','20','40','60','80','100'},'YLim',[0 1]);
box off
set(gcf,'Render','Painter'); saveas(gcf,sprintf('Proportion of memory neurons modulating FR for %s in memory neurons-%s',ParameterName,RegName),'fig'); close all;
% sorting movement-activated units by (FR_movement - FR_still)
FRdiff_activated = FRdiffBetwMoveAndStill(ResponsiveUnitsID(:,2)==1);
[~,SortID_Activated] = sortrows(FRdiff_activated,1,'descend');
SortedActivatedUnitsID = ActivatedUnitsID(SortID_Activated);
ActivatedUnitsFR = FR(ResponsiveUnitsID(:,2)==1,:);
SortedActivatedUnitsFR = ActivatedUnitsFR(SortID_Activated,:);
% sorting movement-inhibited units by (FR_still - FR_movement)
FRdiff_inhibited = FRdiffBetwMoveAndStill(ResponsiveUnitsID(:,2)==2);
[~,SortID_Inhibited] = sortrows(FRdiff_inhibited,1);
SortedInhibitedUnitsID = InhibitedUnitsID(SortID_Inhibited);
InhibitedUnitsFR = FR(ResponsiveUnitsID(:,2)==2,:);
SortedInhibitedUnitsFR = InhibitedUnitsFR(SortID_Inhibited,:);
end




