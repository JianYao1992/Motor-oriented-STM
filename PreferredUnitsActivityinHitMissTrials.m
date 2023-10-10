%% FR in correct and error trials for go-preferred and no-go-preferred neurons
% because of rare miss trials in the DRT, FR in hit, FA and CR without miss trials is analysed

clear; clc; close all;

%% Assignment
CortexName = 'mPFC';
PreferType = 2; % 1 and 2 for Go- and NoGo-preferred neurons, respectively
switch PreferType
    case 1
        CodingType = 'GoPreferredUnits';
    case 2
        CodingType = 'NoGoPreferredUnits';
end
TimeGain = 10;

%% Target directory
CurrPath = uigetdir;
cd(CurrPath);
AllPath = genpath(CurrPath); 
SplitedPath = strsplit(AllPath,';');
SubPath = SplitedPath';
SubPath = SubPath(2:end-1);

tic
%% FR of all neurons and ID of hit, miss, FA, and CR trials 
AllUnitsBinFR = [];
ID_mPFC = [];
ID_aAIC = [];
TrialsID_Hit = [];
TrialsID_Miss = [];
TrialsID_FA = [];
TrialsID_CR = [];
MiceID_BilatmPFC = [{'M19'} {'M20'}];
FileID = 1;
PrevUnitsNum = 0;
for iPath = 1:size(SubPath,1)
    Path = SubPath{iPath,1};
    cd(Path); % enter subpath
    JAVAFiles = dir('*.mat');
    for j = 1:size(JAVAFiles,1)
        Filename{1,FileID} = JAVAFiles(j,1).name;
        load(Filename{1,FileID});
        if ~isempty(SingleUnitList)
            % cross-trial FR
            for itru = 1:size(SingleUnitList,1) % neuron
                tempSuBinFR = [];
                for iTrial = 1:size(TrialMark,1) % trial
                    tempSuBinFR = [tempSuBinFR; Results{1,iTrial}(itru,:)];
                end
                AllUnitsBinFR = [AllUnitsBinFR {tempSuBinFR}];
            end
            % ID of hit, miss, FA, and CR trials
            tempTrialID_Hit = repmat({find(TrialMark(:,4)==1)},1,size(SingleUnitList,1)); 
            tempTrialID_Miss = repmat({find(TrialMark(:,4)==2)},1,size(SingleUnitList,1));
            tempTrialID_FA = repmat({find(TrialMark(:,4)==3)},1,size(SingleUnitList,1)); 
            tempTrialID_CR = repmat({find(TrialMark(:,4)==4)},1,size(SingleUnitList,1)); 
            TrialsID_Hit = [TrialsID_Hit tempTrialID_Hit];
            TrialsID_Miss = [TrialsID_Miss tempTrialID_Miss];
            TrialsID_FA = [TrialsID_FA tempTrialID_FA];
            TrialsID_CR = [TrialsID_CR tempTrialID_CR];
            % ID of mPFC neurons
            tempID_mPFC = [];
            Pos = [];
            for k = 1:length(MiceID_BilatmPFC)
                Pos = [Pos regexp(Filename{1,FileID},MiceID_BilatmPFC{1,k})]; % if non-emp, the mat file is from mouse of which electrodes were implanted into bilateral mPFC
            end
            if ~isempty(Pos)
                tempID_mPFC = find(SingleUnitList(:,1)>=1 & SingleUnitList(:,1)<=16);
            else
                tempID_mPFC = find((SingleUnitList(:,1)>=9 & SingleUnitList(:,1)<=16) | (SingleUnitList(:,1)>=25 & SingleUnitList(:,1)<=32));
            end
            if ~isempty(tempID_mPFC)
                ID_mPFC = [ID_mPFC tempID_mPFC'+PrevUnitsNum];
            end
            % ID of aAIC neurons
            tempID_aAIC = [];
            if ~isempty(Pos) 
                tempID_aAIC = [];
            else
                tempID_aAIC = find((SingleUnitList(:,1)>=1 & SingleUnitList(:,1)<=8) | (SingleUnitList(:,1)>=17 & SingleUnitList(:,1)<=24));
            end
            if ~isempty(tempID_aAIC)
                ID_aAIC = [ID_aAIC tempID_aAIC'+PrevUnitsNum];
            end
            PrevUnitsNum = PrevUnitsNum + size(SingleUnitList,1);
        end
        FileID = FileID + 1;
    end
end
TrialsID = [TrialsID_Hit; TrialsID_Miss; TrialsID_FA; TrialsID_CR]; % row: trial type; column: neuron

%% FRs of target neurons
if strcmp(CortexName,'mPFC')
    UnitsID = ID_mPFC;
else
    UnitsID = ID_aAIC;
end
TrialID_TarReg = TrialsID(:,UnitsID);
BinFR_TarReg = AllUnitsBinFR(1,UnitsID);
WindowSize = 0.1;
SlidingWindow = 0.1;
FR_TarReg = cell(1,size(BinFR_TarReg,2));
for i = 1:size(BinFR_TarReg,2) 
    for j = 1:size(BinFR_TarReg{1,i},1) % trial
        for k = 1:SlidingWindow*TimeGain:size(BinFR_TarReg{1,i},2)-(WindowSize*TimeGain-1)
            FR_TarReg{1,i}(j,k) = sum(BinFR_TarReg{1,i}(j,k:k+WindowSize*TimeGain-1))/WindowSize;
        end
    end
end
% allocate FRs into Go and NoGo trial types
[FR_Go,~,~,FR_NoGo] = ExtractSpecTrialNumForEachNeuron(FR_TarReg,TrialID_TarReg,[],0,0);
% allocate FRs into different trial types
[FR_Hit,FR_Miss,FR_FA,FR_CR] = ExtractSpecTrialNumForEachNeuron(FR_TarReg,TrialID_TarReg,[],0,1);
clearvars -except CortexName CurrPath PreferType CodingType TimeGain FR_TarReg FR_Go FR_NoGo FR_Hit FR_Miss FR_FA FR_CR

%% Event duration
BaseLen = 4;
ShownBaseLen = 2;
BasePeriod = (BaseLen-ShownBaseLen)*TimeGain+1:BaseLen*TimeGain;
SampOdorLen = 1;
SampPeriod = BaseLen*TimeGain+1:(BaseLen+SampOdorLen)*TimeGain;
DelayLen = 4;
TestOdorLen = 0.5;
RespLen = 1;
ITI = 10;
X = -1*BaseLen:0.1:SampOdorLen+DelayLen+TestOdorLen+RespLen+ITI;
X = X(1:size(FR_TarReg{1,1},2)-1)+0.05;
IsClusBasedPermuTest = 1; % permutation test for false and CR trials
WorkerNum = 10;

%% FR modulation /// compare FR between sample, delay period and baseline period (2 s before sample onset) ///
RampingAna = zeros(length(FR_TarReg),6);
RampingAna_Go = zeros(length(FR_TarReg),6);
RampingAna_NoGo = zeros(length(FR_TarReg),6);
IsSigModulation = cell(1,length(FR_TarReg));
for iUnit = 1:length(FR_TarReg)
    % modulation direction and ramping pattern in combined Go and NoGo trials
    [tempNeuIsSigModulation,~,IsRampingNeu,SecondMeanFR] = ModulationAndRampingAnalysis(FR_TarReg{iUnit},BasePeriod,BaseLen,SampOdorLen,DelayLen,TimeGain,SampPeriod);
    RampingAna(iUnit,:) = [IsRampingNeu,SecondMeanFR];
    % modulation direction and ramping pattern in Go trials
    [tempNeuIsSigModulation_Go,~,IsRampingNeu_Go,SecondMeanFR_Go] = ModulationAndRampingAnalysis(FR_Go{iUnit},BasePeriod,BaseLen,SampOdorLen,DelayLen,TimeGain,SampPeriod);
    RampingAna_Go(iUnit,:) = [IsRampingNeu_Go,SecondMeanFR_Go];
    % modulation direction and ramping pattern in NoGo trials
    [tempNeuIsSigModulation_NoGo,~,IsRampingNeu_NoGo,SecondMeanFR_NoGo] = ModulationAndRampingAnalysis(FR_NoGo{iUnit},BasePeriod,BaseLen,SampOdorLen,DelayLen,TimeGain,SampPeriod);
    RampingAna_NoGo(iUnit,:) = [IsRampingNeu_NoGo,SecondMeanFR_NoGo];
    % summarize modulation direction
    IsSigModulation{iUnit} = [tempNeuIsSigModulation; tempNeuIsSigModulation_Go; tempNeuIsSigModulation_NoGo];
end

%% Identify the neuron group /// based on modulation direction ///
% load selectivity result
cd(CurrPath);
load([CortexName 'SelectivityData_CtrlGroup']);
% group ID of all mPFC or aAIC neurons
UnitsGroupID = zeros(length(FR_TarReg),2);
MemUnitsID = vertcat(SustainedUnitID,TransientCodingNeuronID,SwitchedSustainedUnitID); % ID of sustained and transient neurons
SigSelectivity = horzcat(SortedID,TargetBrainSigSelectivity); % selectivity of all mPFC or aAIC neurons
SigSelectivity = sortrows(SigSelectivity); % normal order
SigSelectivity = SigSelectivity(:,2:end); % selectivity in normal order
for iUnit = 1:length(MemUnitsID) % assignment of group ID
    tempUnitID = MemUnitsID(iUnit);
    tempSelectivity = [];
    for iSecond = 1:floor(size(SigSelectivity,2)/TimeGain)
        tempSelectivity(1,iSecond) = mean(SigSelectivity(tempUnitID,TimeGain*(iSecond-1)+1:TimeGain*iSecond),2);
    end
    tempSelectivity = tempSelectivity(:,3:6);
    tempIsSigModulation = IsSigModulation{tempUnitID};
    tempUnitDelayFR = [RampingAna_Go(tempUnitID,3:end); RampingAna_NoGo(tempUnitID,3:end)];
    [SelecDirection,SelecUnitType,SelecUnitGroupID] = JudgeSelecOfExcInhPattern(tempSelectivity,tempIsSigModulation,tempUnitDelayFR);
    UnitsGroupID(tempUnitID,1) = SelecDirection;
    UnitsGroupID(tempUnitID,2) = SelecUnitGroupID;
end

%% Plot averaged performance-dependent FR for different neuron groups
% number of hit, FA, and CR trials
TrialsNum_Hit = cellfun(@(x) size(x,1),FR_Hit,'UniformOutput',1);
TrialsNum_Miss = cellfun(@(x) size(x,1),FR_Miss,'UniformOutput',1);
TrialsNum_FA = cellfun(@(x) size(x,1),FR_FA,'UniformOutput',1);
TrialsNum_CR = cellfun(@(x) size(x,1),FR_CR,'UniformOutput',1);
% compare mean FR in Hit and Miss trials using same trial number
TriNumThres_HitandMiss = [3 3];
BootstrapNum = 1000;
tempNeuID = find(UnitsGroupID(:,1) == PreferType & TrialsNum_Hit(:) >= TriNumThres_HitandMiss(1) & TrialsNum_Miss(:) >= TriNumThres_HitandMiss(2));
save(sprintf('ID-%s-%s_Hit-%d-Miss-%d',CortexName,CodingType,TriNumThres_HitandMiss(1),TriNumThres_HitandMiss(2)),'tempNeuID','-v7.3');
if numel(tempNeuID) > 1
    [Sam1CorrErrIsSig,Sam1CorrErrIsSig_BelowChance,FR_CorrErrTrials] = PlotSelecUnitsHitMissFR(tempNeuID,FR_Hit,FR_Miss,X,SampOdorLen,DelayLen,TestOdorLen,...
        RespLen,size(FR_TarReg{1,1},2),TimeGain,BasePeriod,1,BootstrapNum,[CodingType 'Hit>=' num2str(TriNumThres_HitandMiss(1)) '-Miss>=' num2str(TriNumThres_HitandMiss(2))],IsClusBasedPermuTest,WorkerNum);
    set(gcf,'Render','Painter'); saveas(gcf,[CortexName '-' CodingType '-PerDependFR-Hit-' num2str(TriNumThres_HitandMiss(1)) '-Miss-' num2str(TriNumThres_HitandMiss(2))],'fig'); 
end

% % sample-preferred and delay-activated units
% UnitsGroupMarker = [{'GoSelec-GoExc'};{'GoSelec-GoInh'};{'GoSelec-NogoInh'};{'NogoSelec-NogoExc'};{'NogoSelec-NogoInh'};{'NogoSelec-GoInh'}];
% for iComb = 1:numel(FalseCRTrialsNumCritComb)
%     if PreferType == 1
%         GroupMarkerID = 1;
%     elseif PreferType == 2
%         GroupMarkerID = 4;
%     end
%     tempNeuID = find(UnitsGroupID(:,2) == GroupMarkerID & TrialsNum_FA(:) >= FalseCRTrialsNumCritComb{iComb}(1) & TrialsNum_CR(:) >= FalseCRTrialsNumCritComb{iComb}(2));
%     if numel(tempNeuID) > 1
%         [Sam2CorrErrIsSig,Sam2CorrErrIsSig_BelowChance,NormalizedCorrErrTrialsFR] = PlotSelecUnitsPerfDependFR(tempNeuID,FR_Hit,FR_CR,FR_False,X,SampOdorLen,DelayLen,TestOdorLen,...
%             RespLen,size(FR_TarReg{1,1},2),TimeGain,BasePeriod,RaworNormFR,BootsStrapNum,[UnitsGroupMarker{GroupMarkerID} '-False>=' num2str(FalseCRTrialsNumCritComb{iComb}(1)) '-CR>=' num2str(FalseCRTrialsNumCritComb{iComb}(2))],IsClusBasedPermuTest,WorkerNum);
%         if RaworNormFR == 1
%             saveas(gcf,[CortexName '-' UnitsGroupMarker{GroupMarkerID} '-PerDependFR-False-' num2str(FalseCRTrialsNumCritComb{iComb}(1)) '-CR-' num2str(FalseCRTrialsNumCritComb{iComb}(2))],'fig')
%         elseif RaworNormFR == 2
%             saveas(gcf,[CortexName '-' UnitsGroupMarker{GroupMarkerID} '-PerDepend-NormFR-False-' num2str(FalseCRTrialsNumCritComb{iComb}(1)) '-CR-' num2str(FalseCRTrialsNumCritComb{iComb}(2))],'fig')
%         end
%         close all
%     end
% end
toc





