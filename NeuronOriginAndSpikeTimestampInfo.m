%% Information of neuron origion.

clear; clc; close all;

%% Assignment
group = 'CtrlGroup';
PrevUnitsNum = 0;
BilateralmPFCMiceID = [{'M19'} {'M20'}]; 
FileID = 1;

%% Target directory
CurrPath = uigetdir;
AllPath = genpath(CurrPath);
SplitPath = strsplit(AllPath,';');
SubPath = SplitPath';
SubPath = SubPath(2:end-1);

%% Spike time of each trial and ID of mPFC and aAIC neurons
ID_mPFC = [];
ID_aAIC = [];
AllUnitSpikeTime = [];
AllUnitSpikeTimeIndivTrial = [];
AllUnitTrialMark = [];
MiceID = []; 
LearningDayID = [];
for iPath = 1:size(SubPath,1)
    Path = SubPath{iPath,1};
    cd(Path);
    JAVAFiles = dir('*.mat');
    for j = 1:size(JAVAFiles,1)
        Filename{1,FileID} = JAVAFiles(j,1).name;
        load(Filename{1,FileID});
        if ~isempty(SingleUnitList)
            MiceID = [MiceID; repmat(str2num(Filename{1,FileID}(8:9)),size(SingleUnitList,1),1)];
            LearningDayID = [LearningDayID; repmat(str2num(Filename{1,FileID}(end-4)),size(SingleUnitList,1),1)];
            AllUnitSpikeTime = [AllUnitSpikeTime; SpikeTime];
            % spike time of each trial
            for iRGResults = 1:size(RGResults,1)
                AllUnitSpikeTimeIndivTrial{iRGResults+PrevUnitsNum,1} = RGResults(iRGResults,:);
                AllUnitTrialMark{iRGResults+PrevUnitsNum,1} = TrialMark;
            end
            % ID of mPFC neurons
            Position = [];
            for k = 1:length(BilateralmPFCMiceID)
                Position = [Position regexp(Filename{1,FileID},BilateralmPFCMiceID{1,k})];
            end
            if ~isempty(Position)
                tempmPFCID = find(SingleUnitList(:,1)>=1 & SingleUnitList(:,1)<=16);  % mice with bilateral mPFC recording
            else
                tempmPFCID = find((SingleUnitList(:,1)>=9 & SingleUnitList(:,1)<=16) | (SingleUnitList(:,1)>=25 & SingleUnitList(:,1)<=32));
            end
            if ~isempty(tempmPFCID)
                ID_mPFC = [ID_mPFC tempmPFCID'+PrevUnitsNum];
            end
            % ID of aAIC neurons
            if ~isempty(Position)
                tempAIID = [];
            else
                tempAIID = find((SingleUnitList(:,1)>=1 & SingleUnitList(:,1)<=8) | (SingleUnitList(:,1)>=17 & SingleUnitList(:,1)<=24));
            end
            if ~isempty(tempAIID)
                ID_aAIC = [ID_aAIC tempAIID'+PrevUnitsNum];
            end
            PrevUnitsNum = PrevUnitsNum + size(SingleUnitList,1);
        end
        FileID = FileID + 1;
    end
end

%% Origin information
cd(CurrPath);
OriginInfo = [MiceID LearningDayID];
mPFCOriginInfo = OriginInfo(ID_mPFC,:);
aAICOriginInfo = OriginInfo(ID_aAIC,:);

%% Determine selectivity pattern of each neuron (transient, sustained, or non-memory, transient or sustained), for spike-correlogram analysis
% mPFC
load(['mPFCSelectivityData_' group]);
mPFCUnitID_Sust = SustainedUnitID; 
mPFCUnitID_Tran = [TransientCodingNeuronID; SwitchedSustainedUnitID]; 
for iUnit = 1:size(mPFCOriginInfo,1)
    if ~isempty(find(mPFCUnitID_Sust==iUnit))
        mPFCUnitSelecPat(iUnit,1) = 2;
    elseif ~isempty(find(mPFCUnitID_Tran==iUnit))
        mPFCUnitSelecPat(iUnit,1) = 1;
    else
        mPFCUnitSelecPat(iUnit,1) = 0;
    end
end
mPFCOriginInfo = [mPFCOriginInfo mPFCUnitSelecPat];
% aAIC
load(['aAICSelectivityData_' group]);
aAICUnitID_Sust = SustainedUnitID;
aAICUnitID_Tran = [TransientCodingNeuronID; SwitchedSustainedUnitID]; 
for iUnit = 1:size(aAICOriginInfo,1)
    if ~isempty(find(aAICUnitID_Sust==iUnit))
        aAICUnitSelecPat(iUnit,1) = 2;
    elseif ~isempty(find(aAICUnitID_Tran==iUnit))
        aAICUnitSelecPat(iUnit,1) = 1;
    else
        aAICUnitSelecPat(iUnit,1) = 0;
    end
end
aAICOriginInfo = [aAICOriginInfo aAICUnitSelecPat];
% assign spike time and trial mark with ID of mPFC and aAIC neurons
mPFCSpikeTime = AllUnitSpikeTime(ID_mPFC,:);
aAICSpikeTime = AllUnitSpikeTime(ID_aAIC,:);
mPFCSpikeTimeinIndividualTrial = AllUnitSpikeTimeIndivTrial(ID_mPFC,:);
aAICSpikeTimeinIndividualTrial = AllUnitSpikeTimeIndivTrial(ID_aAIC,:);
mPFCUnitTrialMark = AllUnitTrialMark(ID_mPFC,:);
aAICUnitTrialMark = AllUnitTrialMark(ID_aAIC,:);
save(['NeuronOriginandSpikeRasterInformationfor' group],'mPFCOriginInfo','aAICOriginInfo','mPFCSpikeTime','aAICSpikeTime',...
    'mPFCSpikeTimeinIndividualTrial','aAICSpikeTimeinIndividualTrial','mPFCUnitTrialMark','aAICUnitTrialMark');









