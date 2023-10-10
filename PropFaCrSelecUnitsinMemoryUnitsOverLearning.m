%% Proportion of performance-dependent neurons in Go-preferred or NoGo-preferred neurons
% //////mPFC and aAIC respectively//////

clear; clc; close all;

%% Assignment
Group = 'CtrlGroup';
CortexName = {'mPFC','aAIC'};
PreferType = [1 2]; % 1 and 2 for Go- and NoGo-preferred neurons, respectively
LearningDayID = 1:1:6;
TimeGain = 10;

%% Target directory
CurrPath = uigetdir;
AllPath = genpath(CurrPath);
SplitedPath = strsplit(AllPath,';');
SubPath = SplitedPath';
SubPath = SubPath(2:end-1);

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

for iReg = 1:numel(CortexName)
    %% FRs of target neurons
    if strcmp(CortexName{iReg},'mPFC')
        UnitsID = ID_mPFC;
        IsPerfRelatedUnit.mPFC.GoPreferredUnits = cell(1,numel(LearningDayID));
        IsPerfRelatedUnit.mPFC.NoGoPreferredUnits = cell(1,numel(LearningDayID));
    else
        UnitsID = ID_aAIC;
        IsPerfRelatedUnit.aAIC.GoPreferredUnits = cell(1,numel(LearningDayID));
        IsPerfRelatedUnit.aAIC.NoGoPreferredUnits = cell(1,numel(LearningDayID));
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
    
    %% FR in different tiral types of target neuorns in target learning days
    cd(CurrPath);
    % allocate FRs into Go and NoGo trial types
    [FR_Go,~,~,FR_NoGo] = ExtractSpecTrialNumForEachNeuron(FR_TarReg,TrialID_TarReg,[],0,0);
    % allocate FRs into different trial types
    [FR_Hit,~,FR_FA,FR_CR] = ExtractSpecTrialNumForEachNeuron(FR_TarReg,TrialID_TarReg,[],0,1);
    
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
    
    %% FR modulation ///compare FR between sample, delay period and baseline period (2 s before sample onset)///
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
    load([CortexName{iReg} 'SelectivityData_CtrlGroup']);
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
    
    %% Load origin information
    addpath('X:\20190811DRT_RQP_LaserAndNoLaser\SelectivityAndLaserFiring\Ctrl');
    load(['NeuronOriginandSpikeRasterInformationfor' Group '.mat']); % neuron origin information
    if strcmp(CortexName{iReg},'mPFC')
        OriginInfo = mPFCOriginInfo;
    else
        OriginInfo = aAICOriginInfo;
    end
    
    %% Plot averaged performance-dependent FR for different neuron groups
    % number of hit, FA, and CR trials
    TrialsNum_Hit = cellfun(@(x) size(x,1),FR_Hit,'UniformOutput',1);
    TrialsNum_FA = cellfun(@(x) size(x,1),FR_FA,'UniformOutput',1);
    TrialsNum_CR = cellfun(@(x) size(x,1),FR_CR,'UniformOutput',1);
    % compare mean FR in false and CR trials using same trial number
    TriNumThres_FAandCR = {[3 3],[4 4],[5 5],[6 6],[7 7],[8 8],[9 9],[10 10],[11 11],[12 12],[13 13],[14 14],[15 15],[16 16],[17 17],[18 18],[19 19],[20 20],...
        [21 21],[22 22],[23 23],[24 24],[25 25],[26 26],[27 27],[28 28],[29 29],[30 30],[31 31],[32 32],[33 33],[34 34],[35 35],[36 36],[37 37],...
        [38 38],[39 39],[40 40],[5 10],[5 15],[5 20],[5 25],[5 30],[5 35],[5 40],[5 45],[5 50],[5 55],[5 60],[5 65],[5 70],[5 75],...
        [10 5],[10 15],[10 20],[10 25],[10 30],[10 35],[10 40],[10 45],[10 50],[10 55],[10 60],[10 65],[10 70], [15 5],[15 10],[15 20],...
        [15 25],[15 30],[15 35],[15 40],[15 45],[15 50],[15 55],[15 60],[15 65],[20 5],[20 10],[20 15],[20 25],[20 30],[20 35],[20 40],[20 45],[20 50],...
        [20 55],[20 60],[25 5],[25 10],[25 15],[25 20],[25 30],[25 35],[25 40],[25 45],[25 50],[25 55],[30 5],[30 10],[30 15],[30 20],[30 25],[30 35],...
        [30 40],[30 45],[30 50],[35 5],[35 10],[35 15],[35 20],[35 25],[35 30],[35 40],[35 45],[40 40],[45 5],[45 10],[45 15],[45 20],[45 25],[45 30],[45 35],...
        [50 5],[50 10],[50 15],[50 20],[50 25],[50 30],[55 5],[55 10],[55 15],[55 20],[55 25],[60 5],[60 10],[60 15],[60 20],[65 5],[65 10],[65 15],[70 5],...
        [70 10],[75 5]};
    for iPrefType = 1:numel(PreferType)
        switch PreferType(iPrefType)
            case 1
                CodingType = 'GoPreferredUnits';
            case 2
                CodingType = 'NoGoPreferredUnits';
        end
        for iDay = 1:numel(LearningDayID)
            tempdata = cell(1,numel(TriNumThres_FAandCR));
            for iComb = 1:numel(TriNumThres_FAandCR)
                tempNeuID = find(ismember(OriginInfo(:,2),LearningDayID(iDay)) & UnitsGroupID(:,1)==PreferType(iPrefType) & TrialsNum_FA(:)>=TriNumThres_FAandCR{iComb}(1) & TrialsNum_CR(:)>=TriNumThres_FAandCR{iComb}(2));
                if numel(tempNeuID) >= 3
                    fprintf('LearningDay-%s-%s-FA>=%d-CR>=%d-%sUnitsNum-%d\n',num2str(LearningDayID(iDay)),CortexName{iReg},TriNumThres_FAandCR{iComb}(1),TriNumThres_FAandCR{iComb}(2),CodingType,numel(tempNeuID));
                    tempdata{iComb} = zeros(2,numel(tempNeuID));
                    for iUnit = 1:numel(tempNeuID)
                        fprintf('Analyzing %dth neuron for neural correlates with behavior\n',iUnit);
                        tempdata{iComb}(1,iUnit) = tempNeuID(iUnit);
                        tempdata{iComb}(2,iUnit) = JudgeSingleUnitSelecForFAandCR(tempNeuID(iUnit),FR_CR,FR_FA,BaseLen,SampOdorLen,DelayLen,TimeGain,size(FR_TarReg{1,1},2),80,PreferType(iPrefType));
                    end
                end
            end
            if strcmp(CortexName{iReg},'mPFC')
                if strcmp(CodingType,'GoPreferredUnits')
                    IsPerfRelatedUnit.mPFC.GoPreferredUnits{iDay} = tempdata;
                elseif strcmp(CodingType,'NoGoPreferredUnits')
                    IsPerfRelatedUnit.mPFC.NoGoPreferredUnits{iDay} = tempdata;
                end
            elseif strcmp(CortexName{iReg},'aAIC')
                if strcmp(CodingType,'GoPreferredUnits')
                    IsPerfRelatedUnit.aAIC.GoPreferredUnits{iDay} = tempdata;
                elseif strcmp(CodingType,'NoGoPreferredUnits')
                    IsPerfRelatedUnit.aAIC.NoGoPreferredUnits{iDay} = tempdata;
                end
            end
        end
    end
end
save(['Proportion of performance dependent neurons in memory neurons over learning_RanksumTest-' Group],'IsPerfRelatedUnit','-v7.3');




function IsPerfRelatedUnit = JudgeSingleUnitSelecForFAandCR(UnitID,Samp2CorrTriFR,Samp2ErrorTriFR,BaseLen,SampOdorLen,DelayLen,TimeGain,TriLen,BootsStrapNum,Preference)

IsPerfRelatedUnit = 0;
DelayBinID = (1+SampOdorLen+1):(1+SampOdorLen+DelayLen);

%% Averaged FR in FA and CR trials
UnitFR_CR = Samp2CorrTriFR{UnitID};
UnitFR_FA = Samp2ErrorTriFR{UnitID};

%% Test significance of selectivity
if size(UnitFR_CR,1) >= 2 && size(UnitFR_FA,1) >= 2
    [FR_CR,FR_FA] = CalculateResampledAverCorrErrorFR(UnitFR_CR,UnitFR_FA,BootsStrapNum,TriLen);
    FR_CR = FR_CR(:,(BaseLen-1)*TimeGain+1:end);
    FR_FA = FR_FA(:,(BaseLen-1)*TimeGain+1:end);
    p = [];
    for j = 1:floor(size(FR_CR,2)/TimeGain) % window
        p(1,j) = ranksum(mean(FR_CR(:,1+10*(j-1):10*j),2),mean(FR_FA(:,1+10*(j-1):10*j),2));
    end
    p = p * floor(size(FR_CR,2)/TimeGain);
    SigSelectivity = [];
    SigSelecBinNum = 0;
    for j = 1:floor(size(FR_CR,2)/TimeGain) % window
        if p(j) <= 0.05
            SigSelectivity(1,j) = (mean(mean(FR_FA(:,1+10*(j-1):10*j),2))-mean(mean(FR_CR(:,1+10*(j-1):10*j),2)))/(mean(mean(FR_FA(:,1+10*(j-1):10*j),2))+mean(mean(FR_CR(:,1+10*(j-1):10*j),2)));
            if j >= 1+SampOdorLen && j <= 1+SampOdorLen+DelayLen
                SigSelecBinNum = SigSelecBinNum + 1;
            end
        else
            SigSelectivity(1,j) = 0;
        end
    end
    if SigSelecBinNum == DelayLen
        if (Preference == 1 && nnz(SigSelectivity(DelayBinID)>0) > 0) || (Preference == 2 && nnz(SigSelectivity(DelayBinID)<0) > 0)
            IsPerfRelatedUnit = 1;
        end
    elseif SigSelecBinNum > 0 && SigSelecBinNum < DelayLen
        [~,SigBinID,~,~] = TransientCodingTest(UnitID,FR_FA,FR_CR,find(SigSelectivity~=0));
        SigBinID = SigBinID(SigBinID>=1+SampOdorLen & SigBinID<=1+SampOdorLen+DelayLen);
        SigSelectivity = SigSelectivity(SigBinID);
        if ~isempty(SigSelectivity) && ((Preference==1 && nnz(SigSelectivity>0)>0)||(Preference==2 && nnz(SigSelectivity<0)>0))
            IsPerfRelatedUnit = 1;
        end
    end
end
end


