%% Selectivity for hit vs. CR at population level.

clear; clc; close all;

%% Assignment
Group = 'CtrlGroup';
TargetBrain = 'aAIC'; 
FANumThres = 10;
CRNumThres = 40;
beforeneuronnum = 0; 
BaseLen = 4; ShownBaseLen = 0.5; SampOdorLen = 1; DelayLen = 4; TestOdorLen = 0.5; RespLen = 0.5;
BinNuminWindow = 5;
BilateralmPFCMiceID = [{'M19'} {'M20'}];
FileID = 1;
AllUnitBinFR = []; 
AllUnitTrialID_Hit = [];
AllUnitTrialID_Miss = [];
AllUnitTrialID_FA = []; 
AllUnitTrialNum_FA = [];
AllUnitTrialID_CR = [];
AllUnitTrialNum_CR = [];
ID_mPFC = []; 
ID_aAIC = [];

%% Target directory
CurrPath = uigetdir; 
AllPath = genpath(CurrPath);
SplitPath = strsplit(AllPath,';');
SubPath = SplitPath';
SubPath = SubPath(2:end-1);
 
%% FR and ID of hit, miss, FA, and CR trials for all neurons
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
                AllUnitBinFR = [AllUnitBinFR {tempsingleunitBinnedFR}];
            end
            % ID and number of hit trials
            TrialID_Hit = repmat({find(TrialMark(:,4)==1)},1,size(SingleUnitList,1)); 
            TrialNum_Hit = repmat(nnz(TrialMark(:,4)==1),1,size(SingleUnitList,1));
            % ID and number of miss trials
            TrialID_Miss = repmat({find(TrialMark(:,4)==2)},1,size(SingleUnitList,1)); 
            TrialNum_Miss = repmat(nnz(TrialMark(:,4)==2),1,size(SingleUnitList,1));
            % ID and number of FA trials
            TrialID_FA = repmat({find(TrialMark(:,4)==3)},1,size(SingleUnitList,1));
            TrialNum_FA = repmat(nnz(TrialMark(:,4)==3),1,size(SingleUnitList,1));
            % ID and number of CR trials
            TrialID_CR = repmat({find(TrialMark(:,4)==4)},1,size(SingleUnitList,1));
            TrialNum_CR = repmat(nnz(TrialMark(:,4)==4),1,size(SingleUnitList,1));
            % matrix of all neurons
            AllUnitTrialID_Hit = [AllUnitTrialID_Hit TrialID_Hit]; % hit
            AllUnitTrialID_Miss = [AllUnitTrialID_Miss TrialID_Miss]; % miss
            AllUnitTrialID_FA = [AllUnitTrialID_FA TrialID_FA]; % FA
            AllUnitTrialNum_FA = [AllUnitTrialNum_FA TrialNum_FA];
            AllUnitTrialID_CR = [AllUnitTrialID_CR TrialID_CR]; % CR
            AllUnitTrialNum_CR = [AllUnitTrialNum_CR TrialNum_CR];
            % mPFC neurons ID
            tempID_mPFC = [];
            Position = [];
            for k = 1:length(BilateralmPFCMiceID)
                Position = [Position regexp(Filename{1,FileID},BilateralmPFCMiceID{1,k})]
            end
            if ~isempty(Position)
                tempID_mPFC = find(SingleUnitList(:,1)>=1 & SingleUnitList(:,1)<=16);  % mice with bilateral mPFC recording
            else
                tempID_mPFC = find((SingleUnitList(:,1)>=9 & SingleUnitList(:,1)<=16) | (SingleUnitList(:,1)>=25 & SingleUnitList(:,1)<=32));
            end
            if ~isempty(tempID_mPFC)
                ID_mPFC = [ID_mPFC tempID_mPFC'+beforeneuronnum];
            end
            % aAIC neurons ID
            tempID_aAIC = [];
            if ~isempty(Position)
                tempID_aAIC = [];
            else
                tempID_aAIC = find((SingleUnitList(:,1)>=1 & SingleUnitList(:,1)<=8) | (SingleUnitList(:,1)>=17 & SingleUnitList(:,1)<=24));
            end
            if ~isempty(tempID_aAIC)
                ID_aAIC = [ID_aAIC tempID_aAIC'+beforeneuronnum];
            end
            beforeneuronnum = beforeneuronnum + size(SingleUnitList,1);
        end
        FileID = FileID + 1;
    end
end
cd(CurrPath);
AllUnitTrialID = [AllUnitTrialID_Hit; AllUnitTrialID_Miss; AllUnitTrialID_FA; AllUnitTrialID_CR];

%% FR and ID of trials for target neurons
windowsize = 0.1; sliding = 0.1; % 100-ms window, sliding 100 ms
if strcmp(TargetBrain,'mPFC')
    TargetBrainUnitsID = ID_mPFC;
else
    TargetBrainUnitsID = ID_aAIC;
end
UnitBinFR = AllUnitBinFR(1,TargetBrainUnitsID);
UnitsTrialID = AllUnitTrialID(:,TargetBrainUnitsID);
UnitTrialNum_FA = AllUnitTrialNum_FA(1,TargetBrainUnitsID);
UnitTrialNum_CR = AllUnitTrialNum_CR(1,TargetBrainUnitsID);

%% ID of bin with significant selectivity for hit vs. CR for qualified neurons
[UnitsFR_Hit,~,~,UnitsFR_CR] = ExtractSpecTrialNumForEachNeuron(UnitBinFR,UnitsTrialID,[],0,1);
QualifiedUnitsID = find(UnitTrialNum_FA >= FANumThres & UnitTrialNum_CR >= CRNumThres);
UnitsFR_Hit = UnitsFR_Hit(:,QualifiedUnitsID);
UnitsFR_CR = UnitsFR_CR(:,QualifiedUnitsID);
[UnitsIDwithSelecHitCR,SigBinsID]= GetPutativeSigUnitID(UnitsFR_Hit,UnitsFR_CR,BinNuminWindow,BaseLen,SampOdorLen,DelayLen,TestOdorLen,ShownBaseLen,RespLen,1);

%% ID of sustained (including switched sustained neurons) and transient neurons
SustainedUnitsID = [];
TransientUnitsID = [];
for i = UnitsIDwithSelecHitCR'
    tempSigBinsID = SigBinsID{i};
    tempSigBinsID = tempSigBinsID(tempSigBinsID >= (BaseLen-ShownBaseLen)*10/BinNuminWindow+2 & tempSigBinsID <= (BaseLen + SampOdorLen + DelayLen + TestOdorLen + RespLen)*10/BinNuminWindow);
    TotalBinsNum = (SampOdorLen+DelayLen+TestOdorLen+RespLen)*10/BinNuminWindow;
    if numel(tempSigBinsID) == TotalBinsNum
        SustainedUnitsID = [SustainedUnitsID i];
    else
        TransientUnitsID = [TransientUnitsID i];
    end
end

%% Significant selectivity of transient neurons
NoTransientUnitsID = [];
ZStatis = cell(1,length(QualifiedUnitsID));
P = cell(1,length(QualifiedUnitsID));
for i = TransientUnitsID
    [tempNoTransientCodingNeuronID,tempSigBinID,tempZStatistics,tempP] = TransientCodingTest_forHitCR(i,UnitsFR_Hit{i},UnitsFR_CR{i},SigBinsID{i},BinNuminWindow,BaseLen,ShownBaseLen,SampOdorLen,DelayLen,TestOdorLen,RespLen);
    NoTransientUnitsID = [NoTransientUnitsID; tempNoTransientCodingNeuronID];
    SigBinsID{i} = tempSigBinID;
    ZStatis{i} = tempZStatistics;
    P{i} = tempP;
    disp(strcat('///Finish testing',num2str(find(TransientUnitsID==i)),'th transient coding unit of total',num2str(numel(TransientUnitsID)),' neurons///'));
end
TransientUnitsID = setdiff(TransientUnitsID,NoTransientUnitsID);
CodingUnitsID = horzcat(QualifiedUnitsID(SustainedUnitsID),QualifiedUnitsID(TransientUnitsID));
CodingUnitsID_NeuronsInRange = horzcat(SustainedUnitsID,TransientUnitsID);
NoTransientUnitsID_NeuronsInRange = NoTransientUnitsID;
NoTransientUnitsID = QualifiedUnitsID(NoTransientUnitsID);

%% Save results
save([TargetBrain 'Data_SelectivityforHitandCR-FANumber' num2str(FANumThres) '-CRNumber' num2str(CRNumThres) '-' Group],'CodingUnitsID','CodingUnitsID_NeuronsInRange',...
    'NoTransientUnitsID','NoTransientUnitsID_NeuronsInRange','P','ZStatis','SigBinsID');

 
 
 
 
 
 
 
 
 

