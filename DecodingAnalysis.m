%% Decoding accuracy (template matching)
% addpath(genpath('/gpfsdata/home/jyao/ERCodes/ERProcessingcode')); % ION computing center

clear; clc; close all;

%% Assignment
Target = 'mPFC';
NeuronCodingType = 4;  % 0: no remove; 1: remove sustained neurons; 2: remove transient neurons 3: only sustained 4: only transient 5: only non-memory neurons
switch(NeuronCodingType)
    case 0
        RegionCoding = [Target,' Neurons_All'];
    case 1
        RegionCoding = [Target,' Neurons_Remove sustained'];
    case 2
        RegionCoding = [Target,' Neurons_Remove transient'];
    case 3
        RegionCoding = [Target,' Neurons_Only sustained'];
    case 4
        RegionCoding = [Target,' Neurons_Only transient'];
    case 5
        RegionCoding = [Target,' Neurons_Only nonmemory'];
end
Group = 'CtrlGroup';
WorkerNumber = 10;
NullDistributionNum = 10; 
ResampRunNum = 100; 
IsCTDecoding = 1; % decoding or cross-temporal decoding
beforeneuronnum = 0; 
UnitsNumUsedForDecoding = 150;
TrialNumForEachCondition = 60;
AllUnitBinnedFR = [];
AllUnitTrialID_S1 = [];
AllUnitTrialID_S2 = [];
mPFC_ID = []; aAIC_ID = [];
BilateralmPFCMiceID = [{'M19'} {'M20'}];
TimeGain = 10;

%% Load selectivity
if strcmp(Target,'mPFC')
    SelectivityFile = dir('*mPFCSelectivityData*.mat');
    load(SelectivityFile.name);
else
    SelectivityFile = dir('*aAICSelectivityData*.mat');
    load(SelectivityFile.name);
end

%% FR file
File = dir('*short*.mat');

%% FR and ID of trials for all neurons
disp('-----Step 1: start loading file-----');
for itrm = 1:size(File,1) 
    if size(File,1)== 1
        load(File.name);
    else
        load(File(itrm,1).name);
    end
    disp(['---load File' num2str(itrm)]);
    if ~isempty(SingleUnitList)
        for itru = 1:size(SingleUnitList,1) % neuron
            tempsingleunitBinnedFR = [];
            for itrt = 1:size(TrialMark,1) % trial
                tempsingleunitBinnedFR = [tempsingleunitBinnedFR; Results{1,itrt}(itru,:)];
            end
            AllUnitBinnedFR = [AllUnitBinnedFR {tempsingleunitBinnedFR}];
        end
        S1TrialID = repmat({find(TrialMark(:,2)==1)},1,size(SingleUnitList,1));
        S2TrialID = repmat({find(TrialMark(:,2)==2)},1,size(SingleUnitList,1));
        AllUnitTrialID_S1 = [AllUnitTrialID_S1 S1TrialID];
        AllUnitTrialID_S2 = [AllUnitTrialID_S2 S2TrialID];
        % ID of mPFC neurons
        tempmPFCID = [];
        Position = [];
        for k = 1:length(BilateralmPFCMiceID)
            Position = [Position regexp(File(itrm,1).name,BilateralmPFCMiceID{1,k})]
        end
        if ~isempty(Position)
            tempmPFCID = find(SingleUnitList(:,1)>=1 & SingleUnitList(:,1)<=16);  % mice with bilateral mPFC recording
        else
            tempmPFCID = find((SingleUnitList(:,1)>=9 & SingleUnitList(:,1)<=16) | (SingleUnitList(:,1)>=25 & SingleUnitList(:,1)<=32));
        end
        if ~isempty(tempmPFCID)
            mPFC_ID = [mPFC_ID tempmPFCID'+beforeneuronnum];
        end
        % ID of aAIC neurons
        if ~isempty(Position)
            tempAIID = [];
        else
            tempAIID = find((SingleUnitList(:,1)>=1 & SingleUnitList(:,1)<=8) | (SingleUnitList(:,1)>=17 & SingleUnitList(:,1)<=24));
        end
        if ~isempty(tempAIID)
            aAIC_ID = [aAIC_ID tempAIID'+beforeneuronnum];
        end
        beforeneuronnum = beforeneuronnum + size(SingleUnitList,1);
    end
end
AllUnitTrialID = [AllUnitTrialID_S1; AllUnitTrialID_S2];

%% FR and ID of trials for target neurons
windowsize = 0.1; sliding = 0.1; % 100-msec window, 100-msec sliding window
if strcmp(Target,'mPFC')
    TarUnitsID = mPFC_ID;
else
    TarUnitsID = aAIC_ID;
end
TarUnitBinnedFR = AllUnitBinnedFR(1,TarUnitsID);
NewTarUnitBinnedFR = cell(1,size(TarUnitBinnedFR,2));
for i=1:size(TarUnitBinnedFR,2) % neuron
    for j = 1:size(TarUnitBinnedFR{1,i},1) % trial
        for k = 1:sliding*TimeGain:size(TarUnitBinnedFR{1,i},2)-(windowsize*TimeGain-1) % time bin
            NewTarUnitBinnedFR{1,i}(j,k) = sum(TarUnitBinnedFR{1,i}(j,k:k+windowsize*TimeGain-1))/windowsize;
        end
    end
end
TarFR = cellfun(@(x) x(:,31:end),NewTarUnitBinnedFR,'UniformOutput', 0);
TarUnitsTrialID = AllUnitTrialID(:,TarUnitsID);
if NeuronCodingType == 1
    TarFR(:,SustainedUnitID) = [];
    TarUnitsTrialID(:,SustainedUnitID) = [];
elseif NeuronCodingType == 2
    TarFR(:,[TransientCodingNeuronID;SwitchedSustainedUnitID]) = [];
    TarUnitsTrialID(:,[TransientCodingNeuronID;SwitchedSustainedUnitID]) = [];
elseif NeuronCodingType == 3
    TarFR = TarFR(:,SustainedUnitID);
    TarUnitsTrialID = TarUnitsTrialID(:,SustainedUnitID);
    UnitsNumUsedForDecoding = length(SustainedUnitID)-10;
elseif NeuronCodingType == 4
    TarFR = TarFR(:,[TransientCodingNeuronID;SwitchedSustainedUnitID]);
    TarUnitsTrialID = TarUnitsTrialID(:,[TransientCodingNeuronID;SwitchedSustainedUnitID]);
elseif NeuronCodingType == 5
    TargetID = setdiff(SortedID,vertcat(SustainedUnitID,TransientCodingNeuronID,SwitchedSustainedUnitID));
    TarFR = TarFR(:,TargetID);
    TarUnitsTrialID = TarUnitsTrialID(:,TargetID);
end
poolobj = gcp('nocreate'); % If no pool, do not create new one.
if isempty(poolobj)
    myCluster = parcluster('local'); myCluster.NumWorkers = WorkerNumber; parpool(myCluster,WorkerNumber);
end
tic

%% Decoding analysis for real data
disp('-----Step 2: start decoding analysis for real data-----')
for iNullDis = 1:NullDistributionNum
    f(iNullDis) = parfeval(@CrossValidationTest,1,TarFR,TarUnitsTrialID,TrialNumForEachCondition...
        ,UnitsNumUsedForDecoding,10,0,IsCTDecoding);
    disp(['---' num2str(iNullDis) '-Computing End---'])
end
for iNullDis = 1:NullDistributionNum
    [~,DECODING_RESULTS] = fetchNext(f);  
    save([RegionCoding '_CVRecording-' num2str(iNullDis)], 'DECODING_RESULTS','-v7.3');
end

%% Decoding analysis for shuffled data
disp('-----Step 3: start decoding analysis for shuffled data-----')
for iNullDis = 1:NullDistributionNum
    f_shuffle(iNullDis) = parfeval(@CrossValidationTest,1,TarFR,TarUnitsTrialID,TrialNumForEachCondition...
        ,UnitsNumUsedForDecoding,ResampRunNum,1,IsCTDecoding);
    disp(['---' num2str(iNullDis) '-Computing End---'])
end
for iNullDis = 1:NullDistributionNum
    [~,ShuffleDECODING_RESULTS] = fetchNext(f_shuffle); 
    save([RegionCoding '_CVShuffled-' num2str(iNullDis)], 'ShuffleDECODING_RESULTS','-v7.3');
end
toc
disp(['Run time: ',num2str(toc)]);

%% Summarize real and shuffled decoding results
% real data
switch(NeuronCodingType)
    case 0
        File = dir('*All_CVRecording*.mat');
    case 1
        File = dir('*Remove sustained_CVRecording*.mat');
    case 2
        File = dir('*Remove_transient_CVRecording*.mat');
    case 3
        File = dir('*Only sustained_CVRecording*.mat');
    case 4
        File = dir('*Only transient_CVRecording*.mat');
    case 5
        File = dir('*Only nonmemory_CVRecording*.mat');
end
data = [];
for i=1:size(File,1)
    filename = File(i,1).name;
    tempdata = load(filename);
    data = [data; tempdata.DECODING_RESULTS];
end
if IsCTDecoding == 0
    save([RegionCoding '_NoCTDecoding_Recording_' Group '_n=' num2str(UnitsNumUsedForDecoding)],'data','-v7.3');
else
    save([RegionCoding '_CTDecoding_Recording_' Group '_n=' num2str(UnitsNumUsedForDecoding)],'data','-v7.3');
end
% shufled data
switch(NeuronCodingType)
    case 0
        File = dir('*All_CVShuffled*.mat');
    case 1
        File = dir('*Remove sustained_CVShuffled*.mat');
    case 2
        File = dir('*Remove_transient_CVShuffled*.mat');
    case 3
        File = dir('*Only sustained_CVShuffled*.mat');
    case 4
        File = dir('*Only transient_CVShuffled*.mat');
    case 5
        File = dir('*Only nonmemory_CVShuffled*.mat');
end
data = [];
for i=1:size(File,1)
    filename = File(i,1).name;
    tempdata = load(filename);
    data = [data; tempdata.ShuffleDECODING_RESULTS];
end
if IsCTDecoding == 0
    save([RegionCoding '_NoCTDecoding_Shuffled_' Group '_n=' num2str(UnitsNumUsedForDecoding)],'data','-v7.3');
else
    save([RegionCoding '_CTDecoding_Shuffled_' Group '_n=' num2str(UnitsNumUsedForDecoding)],'data','-v7.3');
end

%% Delete seperate CV and ShuffledCV results
switch(NeuronCodingType)
    case 0
        CVandShuffledCVFile = dir('*All_CV*.mat');
    case 1
        CVandShuffledCVFile = dir('*Remove sustained_CV*.mat');
    case 2
        CVandShuffledCVFile = dir('*Remove transient_CV*.mat');
    case 3
        CVandShuffledCVFile = dir('*Only sustained_CV*.mat');
    case 4
        CVandShuffledCVFile = dir('*Only transient_CV*.mat');
    case 5
        CVandShuffledCVFile = dir('*Only nonmemory_CV*.mat');
end
for iFile = 1:size(CVandShuffledCVFile,1)
    delete(CVandShuffledCVFile(iFile).name);
end

