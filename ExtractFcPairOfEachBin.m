%% Extract FC pairs over the time course of the trial

clear; clc; close all;

%% Assignment
Group = 'ChR2Group';
Target = 'mPFC-aAIC'; % 'Within mPFC', 'Within aAIC', 'mPFC-aAIC', 'aAIC-mPFC'
SecondNum = 7;

%% ID of mPFC and aAIC neurons
load(sprintf('AllNeuronsAUCValue_%s',Group));

%% Neuronal pairs in the target
load(['TestFC_XCORR_stats_-1_0_1msbin_' Group '.mat']);
IsTargetPair = [];
for iPair = 1:length(stats)
    switch Target
        case 'Within mPFC'
            if ~isempty(find(ID_mPFC==stats{iPair}.su1_clusterid)) && ~isempty(find(ID_mPFC==stats{iPair}.su2_clusterid))
                tempIsTargetPair = 1;
            else
                tempIsTargetPair = 0;
            end
        case 'Within aAIC'
            if ~isempty(find(ID_aAIC==stats{iPair}.su1_clusterid)) && ~isempty(find(ID_aAIC==stats{iPair}.su2_clusterid))
                tempIsTargetPair = 1;
            else
                tempIsTargetPair = 0;
            end
        case 'mPFC-aAIC'
            if ~isempty(find(ID_aAIC==stats{iPair}.su1_clusterid)) && ~isempty(find(ID_mPFC==stats{iPair}.su2_clusterid))
                tempIsTargetPair = 1;
            else
                tempIsTargetPair = 0;
            end
        case 'aAIC-mPFC'
            if ~isempty(find(ID_aAIC==stats{iPair}.su1_clusterid)) && ~isempty(find(ID_mPFC==stats{iPair}.su2_clusterid))
                tempIsTargetPair = 1;
            else
                tempIsTargetPair = 0;
            end
    end
    IsTargetPair = [IsTargetPair tempIsTargetPair];
end
TargetPairID = find(IsTargetPair==1);
clearvars stats

%% 'Stat' of FC pair for each bin
TestFcFiles = dir('*TestFC*1ms*.mat');
FcPairOfEachSecond = cell(1,SecondNum);
for iFile = 1:size(TestFcFiles,1) % Different target second
    fprintf('%dth file of total %d files\n',iFile,size(TestFcFiles,1));
    TargetSecond = str2num(TestFcFiles(iFile).name(end-21))+1; % 1:baseline; 2:sample; 3-6: delay; 7:response period.
    load(TestFcFiles(iFile).name);
    Targetstats{iFile} = stats(:,TargetPairID);
    FcPairOfEachSecond{TargetSecond}.s1 = cell(0);
    FcPairOfEachSecond{TargetSecond}.s2 = cell(0);
    for iPair = 1:length(Targetstats{iFile})
        switch Target
            case 'mPFC-aAIC'
                if isfield(Targetstats{iFile}{iPair},'AIs1') && Targetstats{iFile}{iPair}.AIs1 < 0 % significant FC for mPFC-aAIC
                    Targetstats{iFile}{iPair}.PairID = iPair;
                    FcPairOfEachSecond{TargetSecond}.s1{end+1} = Targetstats{iFile}{iPair};
                end
                if isfield(Targetstats{iFile}{iPair},'AIs2') && Targetstats{iFile}{iPair}.AIs2 < 0
                    Targetstats{iFile}{iPair}.PairID = iPair;
                    FcPairOfEachSecond{TargetSecond}.s2{end+1} = Targetstats{iFile}{iPair};
                end
            case 'aAIC-mPFC'
                if isfield(Targetstats{iFile}{iPair},'AIs1') && Targetstats{iFile}{iPair}.AIs1 > 0 % significant FC for aAIC-mPFC
                    Targetstats{iFile}{iPair}.PairID = -1*iPair;
                    FcPairOfEachSecond{TargetSecond}.s1{end+1} = Targetstats{iFile}{iPair};
                end
                if isfield(Targetstats{iFile}{iPair},'AIs2') && Targetstats{iFile}{iPair}.AIs2 > 0
                    Targetstats{iFile}{iPair}.PairID = -1*iPair;
                    FcPairOfEachSecond{TargetSecond}.s2{end+1} = Targetstats{iFile}{iPair};
                end
            otherwise
                if isfield(Targetstats{iFile}{iPair},'AIs1') && Targetstats{iFile}{iPair}.AIs1 ~= 0
                    Targetstats{iFile}{iPair}.PairID = iPair;
                    FcPairOfEachSecond{TargetSecond}.s1{end+1} = Targetstats{iFile}{iPair};
                end
                if isfield(Targetstats{iFile}{iPair},'AIs2') && Targetstats{iFile}{iPair}.AIs2 ~= 0
                    Targetstats{iFile}{iPair}.PairID = iPair;
                    FcPairOfEachSecond{TargetSecond}.s2{end+1} = Targetstats{iFile}{iPair};
                end
        end
    end
    clearvars stats
end

%% Save result
save(['FCofNeurons_PropertyIndividual_1msbin' Target '_' Group '.mat'],'FcPairOfEachSecond','TargetPairID','Targetstats','-v7.3');
