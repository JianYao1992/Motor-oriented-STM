function [IsSigModulation,SecondModuPValue,IsRampingNeu,SecondMeanFR] = ModulationAndRampingAnalysis(tempFR,BaselinePeriod,BaseLen,OdorLen,DelayLen,TimeGain,SamplePeriod)

IsSigModulation = zeros(1,OdorLen+DelayLen);
SecondModuPValue = zeros(1,OdorLen+DelayLen);

%% FR in baseline vs. sample period
[Psample,~,BaselineFR,SampleFR,~] = RankSumSigTest(tempFR,tempFR,BaselinePeriod,SamplePeriod);
SecondModuPValue(1,1) = Psample;
if Psample < 0.05
    SecondModulation(1,1) = ceil(SampleFR-BaselineFR);
    if SampleFR-BaselineFR > 0
        IsSigModulation(1) = 1;
    else
        IsSigModulation(1) = -1;
    end
end

%% FR in baseline vs. delay period 
AverFR = mean(tempFR);
SecondMeanFR = zeros(1,DelayLen);
BonferroniCorrection = round(DelayLen);
for iDelay = 1:DelayLen
    TarDelayPeriod = (BaseLen+OdorLen+iDelay-1)*TimeGain+1:(BaseLen+OdorLen+iDelay)*TimeGain;
    SecondMeanFR(1,iDelay) = mean(AverFR(TarDelayPeriod));
    [p,~,BaselineFR,TarDelayFR,~] = RankSumSigTest(tempFR,tempFR,BaselinePeriod,TarDelayPeriod);
    SecondModuPValue(1,iDelay+1) = p*BonferroniCorrection;
    if p*BonferroniCorrection < 0.05
        SecondModulation(1,iDelay) = ceil(TarDelayFR-BaselineFR);
        if TarDelayFR-BaselineFR > 0
            IsSigModulation(1,iDelay+1) = 1;
        else
            IsSigModulation(1,iDelay+1) = -1;
        end
    end
end

%% Modulation direction (excitatory/inhibitory)
DelayIsSigModulation = IsSigModulation(2:end);
SigModulationBinNum = nnz(DelayIsSigModulation~=0);
ExcBinID = find(DelayIsSigModulation==1);
InhBinID = find(DelayIsSigModulation==-1);
ModulationDirection = 0;
if SigModulationBinNum >= 2 % at least 2 bins with significant activity modulation
    if ~isempty(ExcBinID) && isempty(InhBinID)
        ModulationDirection = 1;
    elseif isempty(ExcBinID) && ~isempty(InhBinID)
        ModulationDirection = -1;
    elseif length(ExcBinID) > length(InhBinID)
        ModulationDirection = 1;
    elseif length(ExcBinID) < length(InhBinID)
        ModulationDirection = -1;
    elseif length(ExcBinID) == length(InhBinID)
        if min(InhBinID) < min(ExcBinID)
            ModulationDirection = -1;
        elseif min(InhBinID) > min(ExcBinID)
            ModulationDirection = 1;
        end
    end    
end
IsSigModulation = [IsSigModulation ModulationDirection];

%% Ramping direction
Diff = diff(SecondMeanFR);
IsRampingNeu=0;
if isempty(find(Diff<0,1))
    IsRampingNeu=1;
elseif isempty(find(Diff>0,1))
    IsRampingNeu=-1;
end

%% Averaged FR during baseline and delay period
SecondMeanFR = [BaselineFR SecondMeanFR];