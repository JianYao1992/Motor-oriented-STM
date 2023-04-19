%% Judge whether FR of delay is significantly higher or lower than baseline activity in hit and CR trials, respectively

clear; clc; close all;

%% Assignment
Group = 'CtrlGroup';
FileID = 1;
WindowSize = 0.1;
Sliding = 0.1;
TimeGain = 10; 
TarPeriod = 31:100; 

%% Target pathway
CurrPath = uigetdir;
AllPath = genpath(CurrPath);
SplitPath = strsplit(AllPath,';'); 
SubPath = SplitPath';
SubPath = SubPath(2:end-1);

%% FR and trial marker of all neurons
AllUnitsBinnedFR = [];
AllUnitsTriMark =[]; 
for iPath = 1:size(SubPath,1) 
    Path = SubPath{iPath,1};
    cd(Path);
    JAVAFiles = dir('*.mat');
    for j = 1:size(JAVAFiles,1)
        Filename{1,FileID} = JAVAFiles(j,1).name;
        load(Filename{1,FileID});
        if ~isempty(SingleUnitList)
            for itru = 1:size(SingleUnitList,1) % neuron
                tempSuBinnedFR = [];
                for itrt = 1:size(TrialMark,1) % trial
                    tempSuBinnedFR = [tempSuBinnedFR; Results{1,itrt}(itru,:)];
                end
                AllUnitsBinnedFR = [AllUnitsBinnedFR {tempSuBinnedFR}];
                AllUnitsTriMark = [AllUnitsTriMark {TrialMark}];
            end
        end
        FileID = FileID + 1;
    end
end
cd(CurrPath);
clearvars -except AllUnitsBinnedFR AllUnitsTrialMark Group

%% Smoothed FR
NewUnitBinnedFR = cell(1,size(AllUnitsBinnedFR,2));
for iUnit = 1:size(AllUnitsBinnedFR,2) % neuron
    for iTrial = 1:size(AllUnitsBinnedFR{1,iUnit},1) % trial
        for k = 1:Sliding*TimeGain:size(AllUnitsBinnedFR{1,iUnit},2)-(WindowSize*TimeGain-1) % time bin
            NewUnitBinnedFR{1,iUnit}(iTrial,k) = sum(AllUnitsBinnedFR{1,iUnit}(iTrial,k:k+WindowSize*TimeGain-1))/WindowSize;
        end
    end
end

%% Averaged FR of each second in hit and CR trials
NewUnitBinnedFR = cellfun(@(x) x(:,TarPeriod),NewUnitBinnedFR,'UniformOutput',0);
FRSigInHitTrials = cell(1,length(NewUnitBinnedFR));
FRSigInCrTrials = cell(1,length(NewUnitBinnedFR)); 
FRinHitTrials = cell(1,length(NewUnitBinnedFR));
FRinCrTrials = cell(1,length(NewUnitBinnedFR));
for iUnit = 1:length(NewUnitBinnedFR)
    % FR in hit trials
    tempFRinHitTrials = NewUnitBinnedFR{1,iUnit}(AllUnitsTriMark{1,iUnit}(:,4)==1,:);
    if ~isempty(tempFRinHitTrials)
        FRSigInHitTrials{iUnit} = [];
        FRinHitTrials{iUnit} = [];
        for iSecond = 3:6
            FRinHitTrials{iUnit}(1,end+1) = mean(mean(tempFRinHitTrials(:,1+TimeGain*(iSecond-1):TimeGain*iSecond),2));
            FRSigInHitTrials{iUnit}(1,end+1) = 4*signrank(mean(tempFRinHitTrials(:,1+TimeGain*(iSecond-1):TimeGain*iSecond),2),mean(tempFRinHitTrials(:,1:10),2))<=0.05;
            if FRSigInHitTrials{iUnit}(1,end) == 1 && mean(mean(tempFRinHitTrials(:,1+TimeGain*(iSecond-1):TimeGain*iSecond),2)) < mean(mean(tempFRinHitTrials(:,1:10),2))
                FRSigInHitTrials{iUnit}(1,end) = -1;
            end
        end
    end
    % FR in CR trials
    tempFRinCrTrials = NewUnitBinnedFR{1,iUnit}(AllUnitsTriMark{1,iUnit}(:,4)==4,:);
    if ~isempty(tempFRinCrTrials)
        FRSigInCrTrials{iUnit} = [];
        for iSecond = 3:6
            FRinCrTrials{iUnit}(1,end+1) = mean(mean(tempFRinCrTrials(:,1+TimeGain*(iSecond-1):TimeGain*iSecond),2));
            FRSigInCrTrials{iUnit}(1,end+1) = 4*signrank(mean(tempFRinCrTrials(:,1+TimeGain*(iSecond-1):TimeGain*iSecond),2),mean(tempFRinCrTrials(:,1:10),2))<=0.05;
            if FRSigInCrTrials{iUnit}(1,end) == 1 && mean(mean(tempFRinCrTrials(:,1+TimeGain*(iSecond-1):TimeGain*iSecond),2)) < mean(mean(tempFRinCrTrials(:,1:10),2))
                FRSigInCrTrials{iUnit}(1,end) = -1;
            end
        end
    end
end
FRinCorrectTrials = vertcat(FRinHitTrials,FRinCrTrials);
FRinCorrectTrials = FRinCorrectTrials';
FRSiginCorrectTrials = vertcat(FRSigInHitTrials,FRSigInCrTrials);
FRSiginCorrectTrials = FRSiginCorrectTrials';

%% Save result
save(['FRSiginCorrectTrialsRelativetoBaseline' Group],'FRSiginCorrectTrials','FRinCorrectTrials','-v7.3');