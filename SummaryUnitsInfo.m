%% Summarize all units information including ID in specific brain, averaged FR, peak-to-valley duration, and sample preference in each sec

clear; clc; close all;

%% Target directory
CurrPath = uigetdir; 
AllPath = genpath(CurrPath);
SplitPath = strsplit(AllPath,';');
SubPath = SplitPath';
SubPath = SubPath(2:end-1);

%% Distinguish FR of mPFC and aAIC neurons
AllUnitsBinnedFR = [];
AllUnitsTrialMark = [];
ID_mPFC = [];
ID_aAIC = []; 
BilateralmPFCMiceID = [{'M19'} {'M20'}];
PrevUnitsNum = 0;
FileID = 1;
Waveform = []; 
for iPath = 1:size(SubPath,1) 
    Path = SubPath{iPath,1};
    cd(Path);
    JAVAFiles = dir('*.mat');
    for j = 1:size(JAVAFiles,1)
        Filename{1,FileID} = JAVAFiles(j,1).name;
        load(Filename{1,FileID});
        if ~isempty(SingleUnitList)
            for itru = 1:size(SingleUnitList,1) % unit
                tempsingleunitBinnedFR = [];
                for itrt = 1:size(TrialMark,1) % trial
                    tempsingleunitBinnedFR = [tempsingleunitBinnedFR; Results{1,itrt}(itru,:)];
                end
                AllUnitsBinnedFR = [AllUnitsBinnedFR {tempsingleunitBinnedFR}];
                AllUnitsTrialMark = [AllUnitsTrialMark {TrialMark}];
            end
            % mPFC neurons ID
            Position = [];
            tempID_mPFC = [];
            for k = 1:length(BilateralmPFCMiceID)
                Position = [Position regexp(Filename{1,FileID},BilateralmPFCMiceID{1,k})]
            end
            if ~isempty(Position)
                tempID_mPFC = find(SingleUnitList(:,1)>=1 & SingleUnitList(:,1)<=16);  % mice with bilateral mPFC recording
            else
                tempID_mPFC = find((SingleUnitList(:,1)>=9 & SingleUnitList(:,1)<=16) | (SingleUnitList(:,1)>=25 & SingleUnitList(:,1)<=32));
            end
            if ~isempty(tempID_mPFC)
                ID_mPFC = [ID_mPFC tempID_mPFC'+PrevUnitsNum];
            end
            % aAIC neurons ID
            tempID_aAIC = [];
            if ~isempty(Position)
                tempID_aAIC = [];
            else
                tempID_aAIC = find((SingleUnitList(:,1)>=1 & SingleUnitList(:,1)<=8) | (SingleUnitList(:,1)>=17 & SingleUnitList(:,1)<=24));
            end
            if ~isempty(tempID_aAIC)
                ID_aAIC = [ID_aAIC tempID_aAIC'+PrevUnitsNum];
            end
            PrevUnitsNum = PrevUnitsNum + size(SingleUnitList,1);
            % summarize waveform
            Waveform = [Waveform; NewWaveForm];
        end
        FileID = FileID + 1;
    end
end
cd(CurrPath);
clearvars -except AllUnitsBinnedFR AllUnitsTrialMark ID_mPFC ID_aAIC Waveform

%% FR from baseline (1 sec before sample onset) to 0.5 sec after test offset
Group = 'CtrlGroup';
windowsize = 0.1; % 100-msec analysis window
sliding = 0.1; % 100-msec sliding window
TimeGain = 10;
SelPeriod = 31:100;
% averaged FR of mPFC neurons
UnitBinnedFR_mPFC = AllUnitsBinnedFR(1,ID_mPFC);
TrialMark_mPFC = AllUnitsTrialMark(1,ID_mPFC);
NewUnitBinnedFR_mPFC = cell(1,size(UnitBinnedFR_mPFC,2));
for iUnit = 1:size(UnitBinnedFR_mPFC,2)
    for iTrial = 1:size(UnitBinnedFR_mPFC{1,iUnit},1) % trial
        for k = 1:sliding*TimeGain:size(UnitBinnedFR_mPFC{1,iUnit},2)-(windowsize*TimeGain-1) % bin
            NewUnitBinnedFR_mPFC{1,iUnit}(iTrial,k) = sum(UnitBinnedFR_mPFC{1,iUnit}(iTrial,k:k+windowsize*TimeGain-1))/windowsize;
        end
    end
end
MeanFR_mPFC = cellfun(@(x) mean(mean(x(:,SelPeriod))),NewUnitBinnedFR_mPFC,'UniformOutput',1);
% averaged FR of aAIC neurons
UnitBinnedFR_aAIC = AllUnitsBinnedFR(1,ID_aAIC);
TrialMark_aAIC = AllUnitsTrialMark(1,ID_aAIC);
NewUnitBinnedFR_aAIC = cell(1,size(UnitBinnedFR_aAIC,2));
for iUnit=1:size(UnitBinnedFR_aAIC,2) 
    for iTrial = 1:size(UnitBinnedFR_aAIC{1,iUnit},1) % trial
        for k = 1:sliding*TimeGain:size(UnitBinnedFR_aAIC{1,iUnit},2)-(windowsize*TimeGain-1) % bin
            NewUnitBinnedFR_aAIC{1,iUnit}(iTrial,k) = sum(UnitBinnedFR_aAIC{1,iUnit}(iTrial,k:k+windowsize*TimeGain-1))/windowsize;
        end
    end
end
MeanFR_aAIC = cellfun(@(x) mean(mean(x(:,SelPeriod))),NewUnitBinnedFR_aAIC,'UniformOutput',1);

%% Peak-to-valley duration
AllUnitsPeakTroughDuration = []; AllUnitsFWHM = [];
for iUnit = 1:size(Waveform,1)
    [PeakTroughDuration,LargestChannelWF,~,FWHM,~,~] = PeakTroughDurationAnalysis(Waveform(iUnit,:));
    AllUnitsPeakTroughDuration(1,iUnit) = PeakTroughDuration;
    AllUnitsFWHM(1,iUnit) = FWHM;
    AllUnitsLargestChannelWF{iUnit} = LargestChannelWF;
end
% mPFC 
UnitsWFPeaktoValleyDuration_mPFC = AllUnitsPeakTroughDuration(1,ID_mPFC);
UnitsWFFWHM_mPFC = AllUnitsFWHM(1,ID_mPFC);
UnitsLargetstChannelWF_mPFC = AllUnitsLargestChannelWF(:,ID_mPFC);
% aAIC
UnitsWFPeaktoValleyDuration_aAIC = AllUnitsPeakTroughDuration(1,ID_aAIC);
UnitsWFFWHM_aAIC = AllUnitsFWHM(1,ID_aAIC);
UnitsLargetstChannelWF_aAIC = AllUnitsLargestChannelWF(:,ID_aAIC);

%% Sample preference in each sec /// sample, delay, and test period ///
SecNum = 7;
% mPFC
load(['mPFCSelectivityData_' Group]);
load(['mPFCFRinS1S2_' Group]);
WhetherMemInEachSec_mPFC = cell(1,length(TargetBrainUnitsSigBinID));
for iUnit = 1:length(TargetBrainUnitsSigBinID)
    WhetherMemInEachSec_mPFC{iUnit} = zeros(1,SecNum);
    if ~isempty(TargetBrainUnitsSigBinID{iUnit})
        WhetherMemInEachSec_mPFC{iUnit}(TargetBrainUnitsSigBinID{iUnit}) = 1;
        % S2 preference
        for iSigBin = 1:length(TargetBrainUnitsSigBinID{iUnit})
            if mean(mean(TargetBrainUnitsFRinS1{iUnit}(:,(TimeGain*(TargetBrainUnitsSigBinID{iUnit}(iSigBin)-1)+11):(TimeGain*TargetBrainUnitsSigBinID{iUnit}(iSigBin)+10))))...
                    < mean(mean(TargetBrainUnitsFRinS2{iUnit}(:,(TimeGain*(TargetBrainUnitsSigBinID{iUnit}(iSigBin)-1)+11):(TimeGain*TargetBrainUnitsSigBinID{iUnit}(iSigBin)+10))))
                WhetherMemInEachSec_mPFC{iUnit}(TargetBrainUnitsSigBinID{iUnit}(iSigBin)) = 2;
            end
        end
    end
end
% aAIC
load(['aAICSelectivityData_' Group]); 
load(['aAICFRinS1S2_' Group]);
WhetherMemInEachSec_aAIC = cell(1,length(TargetBrainUnitsSigBinID));
for iUnit = 1:length(TargetBrainUnitsSigBinID)
    WhetherMemInEachSec_aAIC{iUnit} = zeros(1,SecNum);
    if ~isempty(TargetBrainUnitsSigBinID{iUnit})
        WhetherMemInEachSec_aAIC{iUnit}(TargetBrainUnitsSigBinID{iUnit}) = 1;
        for iSigBin = 1:length(TargetBrainUnitsSigBinID{iUnit})
            if mean(mean(TargetBrainUnitsFRinS1{iUnit}(:,(TimeGain*(TargetBrainUnitsSigBinID{iUnit}(iSigBin)-1)+11):(TimeGain*TargetBrainUnitsSigBinID{iUnit}(iSigBin)+10))))...
                    < mean(mean(TargetBrainUnitsFRinS2{iUnit}(:,(TimeGain*(TargetBrainUnitsSigBinID{iUnit}(iSigBin)-1)+11):(TimeGain*TargetBrainUnitsSigBinID{iUnit}(iSigBin)+10))))
                WhetherMemInEachSec_aAIC{iUnit}(TargetBrainUnitsSigBinID{iUnit}(iSigBin)) = 2;
            end
        end
    end
end

%% Averaged FR in each sec of hit or CR trials
NewUnitBinnedFR_mPFC = cellfun(@(x) x(:,SelPeriod),NewUnitBinnedFR_mPFC,'UniformOutput',0);
NewUnitBinnedFR_aAIC = cellfun(@(x) x(:,SelPeriod),NewUnitBinnedFR_aAIC,'UniformOutput',0);
FrHitTrials_mPFC = cell(1,length(ID_mPFC)); 
FrCrTrials_mPFC = cell(1,length(ID_mPFC)); 
FrHitTrials_aAIC = cell(1,length(ID_aAIC)); 
FrCrTrials_aAIC = cell(1,length(ID_aAIC));
% mPFC
for iUnit = 1:length(ID_mPFC)
    % FR in hit trials
    tempFRinHitTrials = NewUnitBinnedFR_mPFC{1,iUnit}(TrialMark_mPFC{1,iUnit}(:,4)==1,:);
    if ~isempty(tempFRinHitTrials)
        for iSecond = 1:floor(size(tempFRinHitTrials,2)/TimeGain)
            FrHitTrials_mPFC{iUnit}(1,iSecond) = mean(mean(tempFRinHitTrials(:,1+TimeGain*(iSecond-1):TimeGain*iSecond)));
        end
    else
        FrHitTrials_mPFC{iUnit} = [];
    end
    % FR in CR trials
    tempFRinCrTrials = NewUnitBinnedFR_mPFC{1,iUnit}(TrialMark_mPFC{1,iUnit}(:,4)==4,:);
    if ~isempty(tempFRinCrTrials)
        for iSecond = 1:floor(size(tempFRinCrTrials,2)/TimeGain)
            FrCrTrials_mPFC{iUnit}(1,iSecond) = mean(mean(tempFRinCrTrials(:,1+TimeGain*(iSecond-1):TimeGain*iSecond)));
        end
    else
        FrCrTrials_mPFC{iUnit} = [];
    end
end
% aAIC
for iUnit = 1:length(ID_aAIC)
    % FR in hit trials
    tempFRinHitTrials = NewUnitBinnedFR_aAIC{1,iUnit}(TrialMark_aAIC{1,iUnit}(:,4)==1,:);
    if ~isempty(tempFRinHitTrials)
        for iSecond = 1:floor(size(tempFRinHitTrials,2)/TimeGain)
            FrHitTrials_aAIC{iUnit}(1,iSecond) = mean(mean(tempFRinHitTrials(:,1+TimeGain*(iSecond-1):TimeGain*iSecond)));
        end
    else
        FrHitTrials_aAIC{iUnit} = [];
    end
    % FR in CR trials
    tempFRinCrTrials = NewUnitBinnedFR_aAIC{1,iUnit}(TrialMark_aAIC{1,iUnit}(:,4)==4,:);
    if ~isempty(tempFRinCrTrials)
        for iSecond = 1:floor(size(tempFRinCrTrials,2)/TimeGain)
            FrCrTrials_aAIC{iUnit}(1,iSecond) = mean(mean(tempFRinCrTrials(:,1+TimeGain*(iSecond-1):TimeGain*iSecond)));
        end
    else
        FrCrTrials_aAIC{iUnit} = [];
    end
end

%% Construct matrix /// first column: ID in all brain areas, FR, peak-to-valley duration, FWHM; second column: waveform shape; third column: selectivity in each sec ///
% /// forth colmn: brain area ///
UnitsSummary = struct();
% mPFC
UnitsSummary.mPFC = cell(length(ID_mPFC),4);
for iUnit = 1:length(ID_mPFC)
    UnitsSummary.mPFC{iUnit,1} = horzcat(ID_mPFC(iUnit),MeanFR_mPFC(iUnit),UnitsWFPeaktoValleyDuration_mPFC(iUnit),UnitsWFFWHM_mPFC(iUnit));
    UnitsSummary.mPFC{iUnit,2} = UnitsLargetstChannelWF_mPFC{iUnit};
    UnitsSummary.mPFC{iUnit,3} = WhetherMemInEachSec_mPFC{iUnit};
    UnitsSummary.mPFC{iUnit,4} = FrHitTrials_mPFC{iUnit};
    UnitsSummary.mPFC{iUnit,5} = FrCrTrials_mPFC{iUnit};
    UnitsSummary.mPFC{iUnit,6} = 'mPFC';
end
% aAIC
UnitsSummary.aAIC = cell(length(ID_aAIC),4);
for iUnit = 1:length(ID_aAIC)
    UnitsSummary.aAIC{iUnit,1} = horzcat(ID_aAIC(iUnit),MeanFR_aAIC(iUnit),UnitsWFPeaktoValleyDuration_aAIC(iUnit),UnitsWFFWHM_aAIC(iUnit));
    UnitsSummary.aAIC{iUnit,2} = UnitsLargetstChannelWF_aAIC{iUnit};
    UnitsSummary.aAIC{iUnit,3} = WhetherMemInEachSec_aAIC{iUnit};
    UnitsSummary.aAIC{iUnit,4} = FrHitTrials_aAIC{iUnit};
    UnitsSummary.aAIC{iUnit,5} = FrCrTrials_aAIC{iUnit};
    UnitsSummary.aAIC{iUnit,6} = 'aAIC';
end
save(['UnitsSummary_' Group],'UnitsSummary');