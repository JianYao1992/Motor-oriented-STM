%% ID of Go-preferred and No-Go-preferred neurons in mPFC and aAIC.

clear; clc; close all;

%% Assignment
Group = 'CtrlGroup';
Reg = {'mPFC','aAIC'};
PreferCodingNeuronID = struct();
% task event duration
TimeGain = 10;
BaseLen = 2;
ShownBaseLen = 2;
BasePeriod = (BaseLen-ShownBaseLen)*TimeGain+1:BaseLen*TimeGain;
SampOdorLen = 1;
SampPeriod = BaseLen*TimeGain+1:(BaseLen+SampOdorLen)*TimeGain;
DelayLen = 4;
WindowSize = 0.1;
SlidingWindow = 0.1;

%% Enter path
CurrPath = uigetdir;
cd(CurrPath);

%% FR of all neurons in mPFC and aAIC
for iReg = 1:numel(Reg)
    load(sprintf('%sFRinS1S2_%s.mat',Reg{iReg},Group));
    FR_Go = TargetBrainUnitsFRinS1;
    FR_Nogo = TargetBrainUnitsFRinS2;
    FR_TarReg = cellfun(@(x,y) vertcat(x,y),FR_Go,FR_Nogo,'UniformOutput',0);
    
    %% FR modulation, comparing the activity between sample, delay period and baseline period
    RampingAna = zeros(length(FR_TarReg),6);
    RampingAna_Go = zeros(length(FR_TarReg),6);
    RampingAna_Nogo = zeros(length(FR_TarReg),6);
    IsSigModulation = cell(1,length(FR_TarReg));
    for iUnit = 1:numel(FR_TarReg)
        % modulation direction and ramping pattern in all trials
        [tempNeuIsSigModulation,~,IsRampingNeu,SecondMeanFR] = ModulationAndRampingAnalysis(FR_TarReg{iUnit},BasePeriod,BaseLen,SampOdorLen,DelayLen,TimeGain,SampPeriod);
        RampingAna(iUnit,:) = [IsRampingNeu,SecondMeanFR];
        % modulation direction and ramping pattern in Go trials
        [tempNeuIsSigModulation_Go,~,IsRampingNeu_Go,SecondMeanFR_Go] = ModulationAndRampingAnalysis(FR_Go{iUnit},BasePeriod,BaseLen,SampOdorLen,DelayLen,TimeGain,SampPeriod);
        RampingAna_Go(iUnit,:) = [IsRampingNeu_Go,SecondMeanFR_Go];
        % modulation direction and ramping pattern in NoGo trials
        [tempNeuIsSigModulation_Nogo,~,IsRampingNeu_Nogo,SecondMeanFR_Nogo] = ModulationAndRampingAnalysis(FR_Nogo{iUnit},BasePeriod,BaseLen,SampOdorLen,DelayLen,TimeGain,SampPeriod);
        RampingAna_Nogo(iUnit,:) = [IsRampingNeu_Nogo SecondMeanFR_Nogo];
        % summarize modulation direction
        IsSigModulation{iUnit}=[tempNeuIsSigModulation;tempNeuIsSigModulation_Go;tempNeuIsSigModulation_Nogo];
    end
    
    %% Identify neuron groups, based on modulation direction
    % load selectivity result
    load([Reg{iReg} 'SelectivityData_' Group]);
    CodingUnitsID = vertcat(SustainedUnitID,TransientCodingNeuronID,SwitchedSustainedUnitID); % ID of memory neurons
    SigSelectivity = horzcat(SortedID,TargetBrainSigSelectivity); % selectivity of all neurons in mPFC and aAIC
    SigSelectivity = sortrows(SigSelectivity);
    SigSelectivity = SigSelectivity(:,2:end);
    % whether memory neuron belongs to Go- or NoGo-preferred neuron 
    UnitsGroupID = struct();
    UnitsGroupID.GoPrefUnitsID = [];
    UnitsGroupID.NoGoPrefUnitsID = [];
    UnitsGroupID.totalNum = numel(FR_TarReg);
    for iUnit = 1:length(CodingUnitsID) % assignment of group ID
        tempUnitID = CodingUnitsID(iUnit);
        tempSelectivity = [];
        for iSecond = 1:floor(size(SigSelectivity,2)/TimeGain)
            tempSelectivity(1,iSecond) = mean(SigSelectivity(tempUnitID,TimeGain*(iSecond-1)+1:TimeGain*iSecond),2);
        end
        tempSelectivity = tempSelectivity(:,3:6);
        tempIsSigModulation = IsSigModulation{tempUnitID};
        tempUnitDelayFR = [RampingAna_Go(tempUnitID,3:end); RampingAna_Nogo(tempUnitID,3:end)];
        [SelecDirection,SelecUnitType,SelecUnitGroupID] = JudgeSelecOfExcInhPattern(tempSelectivity,tempIsSigModulation,tempUnitDelayFR);
        if SelecDirection == 1
            UnitsGroupID.GoPrefUnitsID = [UnitsGroupID.GoPrefUnitsID tempUnitID];
        elseif SelecDirection == 2
            UnitsGroupID.NoGoPrefUnitsID = [UnitsGroupID.NoGoPrefUnitsID tempUnitID];
        end
    end
    if strcmp(Reg{iReg},'mPFC')
        PreferCodingNeuronID.mPFC = UnitsGroupID;
    elseif strcmp(Reg{iReg},'aAIC')
        PreferCodingNeuronID.aAIC = UnitsGroupID;
    end
end
save(sprintf('Go-NoGo-PreferredUnitsID-%s',Group),'PreferCodingNeuronID','-v7.3');


