%% Decoding accuracy of FCSP events in congruent active and incongruent inactive pairs

clear; clc; close all;

%% Assignment
Reg = 'Within aAIC';
PairsType = [1 4];
PairsName = {'CongAct','InconInact'};
BaseLen = 4; 
SampLen = 1;

%% ID of FC neuronal pairs
load(sprintf('Functional Coupling Coding_5pairs_%s_CtrlGroup.mat',Reg));
AllTypesFCPairID = AllTypesAUC;
for iType = 1:size(AllTypesFCPairID,2)
    for iBin = 1:size(AllTypesFCPairID{iType},2)
        AllTypesFCPairID{iType}{iBin} = horzcat(iBin*ones(size(AllTypesFCPairID{iType}{iBin},1),1),AllTypesFCPairID{iType}{iBin});
    end
    AllTypesFCPairID{iType} = vertcat(AllTypesFCPairID{iType}{:});
end

%% Minimal number of FC pairs between specific types
BsSize = min(horzcat(size(AllTypesFCPairID{1,PairsType(1)},1),size(AllTypesFCPairID{1,PairsType(end)},1)));
BsSize = 10*floor(BsSize/10);

%% Decoding analysis of FCSP events of congruent active and incongruent inactive pairs
DecodingAnalysisForFCSP('CtrlGroup',AllTypesFCPairID,BaseLen,SampLen,PairsType,PairsName,BsSize);

function DecodingAnalysisForFCSP(Group,FCPairID,BaseLen,SampLen,PairsType,PairsName,BootstrapSize)

FCPairID = FCPairID(:,PairsType);

%% Load ID of mPFC and aAIC neurons
load(sprintf('AllNeuronsAUCValue_%s.mat',Group));

%% Load trial ID and trial-based spike rasters of mPFC and aAIC neurons
load(sprintf('NeuronOriginandSpikeRasterInformationfor%s.mat',Group));

%% Trial ID and FCSP events for each FC neuronal pair
TrialID = cell(1,size(FCPairID,2));
FCrate = cell(1,size(FCPairID,2));
for iType = 1:size(FCPairID,2)
    TrialID{iType} = cell(2,size(FCPairID{iType},1));
    FCrate{iType} = cell(1,size(FCPairID{iType},1));
    for iFC = 1:size(FCPairID{iType},1)
        tempID_Bin = FCPairID{iType}(iFC,1);
        tempID_LU = FCPairID{iType}(iFC,2);
        tempID_FU = FCPairID{iType}(iFC,3);
        if ismember(tempID_LU,ID_mPFC)
            tempReg_LU = 'mPFC';
            tempIDinReg_LU = find(ID_mPFC==tempID_LU);
            tempTrialMark = mPFCUnitTrialMark{tempIDinReg_LU};
            tempTS_LU = mPFCSpikeTimeinIndividualTrial{tempIDinReg_LU}; % spike raster
        elseif ismember(tempID_LU,ID_aAIC)
            tempReg_LU = 'aAIC';
            tempIDinReg_LU = find(ID_aAIC==tempID_LU);
            tempTrialMark = aAICUnitTrialMark{tempIDinReg_LU};
            tempTS_LU = aAICSpikeTimeinIndividualTrial{tempIDinReg_LU}; 
        end
        if ismember(tempID_FU,ID_mPFC)
            tempReg_FU = 'mPFC';
            tempIDinReg_FU = find(ID_mPFC==tempID_FU);
            tempTS_FU = mPFCSpikeTimeinIndividualTrial{tempIDinReg_FU};
        elseif ismember(tempID_FU,ID_aAIC)
            tempReg_FU = 'aAIC';
            tempIDinReg_FU = find(ID_aAIC==tempID_FU);
            tempTS_FU = aAICSpikeTimeinIndividualTrial{tempIDinReg_FU}; 
        end
        TrialID{iType}{1,iFC} = find(tempTrialMark(:,2)==1); % ID of Go trials
        TrialID{iType}{2,iFC} = find(tempTrialMark(:,2)==2); % ID of NoGo trials
        FCrate{iType}{iFC} = CalculateFcspRateOfIndividualFcPair(BaseLen,SampLen,tempID_Bin,tempTS_LU,tempTS_FU);
    end
    UsedNeuNum = 2:1:BootstrapSize;
    for iNum = 1:numel(UsedNeuNum)
        % parallel pool
        WorkerNumber = 10;
        NullDistributionNum = 10;
        TrialNumForEachCondition = 75;
        ResamplingNum = 10;
        poolobj = gcp('nocreate');
        if isempty(poolobj)
            myCluster=parcluster('local'); myCluster.NumWorkers=WorkerNumber; parpool(myCluster,WorkerNumber);
        end
        tic
        fprintf('-----Start decoding analysis for real data-UsedNeuNum-%d-%s-%s-----\n',UsedNeuNum(iNum),PairsName{iType},Group);
        for iNullDis = 1:NullDistributionNum
            f(iNullDis) = parfeval(@CrossValidationTestForFcsp,1,FCrate{iType},TrialID{iType},TrialNumForEachCondition...
                ,UsedNeuNum(iNum),ResamplingNum,0);
            disp(['---' num2str(iNullDis) '-Computing End---'])
        end
        for iNullDis = 1:NullDistributionNum
            [~,DecodingResults] = fetchNext(f); % collect the results as they become available.
            save([PairsName{iType} '-FCSPsRecording-' num2str(iNullDis)], 'DecodingResults','-v7.3');
        end
        
        %% Summarize decoding result
        File = dir('*FCSPsRecording*.mat');
        data = [];
        for i=1:size(File,1)
            filename = File(i,1).name;
            tempdata = load(filename);
            data = [data; tempdata.DecodingResults];
        end
        save([PairsName{iType} '-FCSPsDecoding-n=' num2str(UsedNeuNum(iNum)) '-' Group],'data','-v7.3');
        
        %% Delete seperate decoding files
        File = dir('*FCSPsRecording*.mat');
        for iFile = 1:size(File,1)
            delete(File(iFile).name);
        end
        toc
    end
end
end







