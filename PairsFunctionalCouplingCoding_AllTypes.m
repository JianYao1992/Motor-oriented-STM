%% Memroy-coding ability of FC neuronal pairs for all types

clear; clc; close all;

%% Assignment
Group = 'ChR2Group'; 
Target = 'mPFC-aAIC'; % 'Within mPFC','Within aAIC', 'mPFC-aAIC', and 'aAIC-mPFC'
BinSize = 1; 
BinRange = 1:BinSize:4;
ResamplingTimes = 100;
if strcmp(Target,'Within mPFC')
    targetregion1 = 1; targetregion2 = 1;
elseif strcmp(Target,'Within aAIC')
    targetregion1 = 2; targetregion2 = 2;
elseif strcmp(Target,'mPFC-aAIC') 
    targetregion1 = 1; targetregion2 = 2;
elseif strcmp(Target,'aAIC-mPFC') 
    targetregion1 = 2; targetregion2 = 1;
end

%% Parallel pool
WorkerNum = 10; 
poolobj = gcp('nocreate');
if isempty(poolobj)
    myCluster = parcluster('local'); myCluster.NumWorkers = WorkerNum; parpool(myCluster,WorkerNum);
end

%% Spike rasters in each trial of all neurons
load(sprintf('AllNeuronsAUCValue_%s.mat',Group));
load(sprintf('NeuronOriginandSpikeRasterInformationfor%s.mat',Group));
% neuron ID
ID = [ID_mPFC'; ID_aAIC']; 
ID = num2cell(ID);
% trial marker
UnitTrialMark = horzcat(ID,[mPFCUnitTrialMark; aAICUnitTrialMark]);
SortedUnitTrialMark  = sortrows(UnitTrialMark,1);
% spike raster timestamp within single trial
SpikeTimeIndividualTrial = horzcat(ID,[mPFCSpikeTimeinIndividualTrial; aAICSpikeTimeinIndividualTrial]);
SortedSpikeTimeIndividualTrial = sortrows(SpikeTimeIndividualTrial,1);
clear UnitTrialMark SpikeTimeIndividualTrial

%% Identify various FC neuronal pairs in each delay bin
FCPairID = cell(1,length(BinRange)); 
UniqueFCPairID = cell(1,length(BinRange)); 
for bin = 1:length(BinRange)
    conn_chain_target_Go = [];
    conn_chain_target_Nogo = [];
    load(['1104_conn_chain_delay4s_' num2str(BinRange(bin)) '_' num2str(BinRange(bin)+BinSize) '_GoNogo_' Group '.mat']);
    % Go trials
    for iPair = 1:size(conn_chain_go,1)
        % change the display of FC pairs into positive direction
        if conn_chain_go(iPair,3) < 0
            temp = conn_chain_go(iPair,1); conn_chain_go(iPair,1) = conn_chain_go(iPair,2); conn_chain_go(iPair,2) = temp;
            conn_chain_go(iPair,3) = -1*conn_chain_go(iPair,3);
            temp = pref_chain_go(iPair,1:7); pref_chain_go(iPair,1:7) = pref_chain_go(iPair,8:14); pref_chain_go(iPair,8:14) = temp;
            temp = reg_chain_go(iPair,1); reg_chain_go(iPair,1) = reg_chain_go(iPair,2); reg_chain_go(iPair,2) = temp;
            reg_chain_go(iPair,3) = -1*reg_chain_go(iPair,3);
        end
    end
    % NoGo trials
    for iPair = 1:size(conn_chain_nogo,1)
        % change the display of pairs showing functional coupling
        if conn_chain_nogo(iPair,3) < 0
            temp = conn_chain_nogo(iPair,1); conn_chain_nogo(iPair,1) = conn_chain_nogo(iPair,2); conn_chain_nogo(iPair,2) = temp;
            conn_chain_nogo(iPair,3) = -1*conn_chain_nogo(iPair,3);
            temp = pref_chain_nogo(iPair,1:7); pref_chain_nogo(iPair,1:7) = pref_chain_nogo(iPair,8:14); pref_chain_nogo(iPair,8:14) = temp;
            temp = reg_chain_nogo(iPair,1); reg_chain_nogo(iPair,1) = reg_chain_nogo(iPair,2); reg_chain_nogo(iPair,2) = temp;
            reg_chain_nogo(iPair,3) = -1*reg_chain_nogo(iPair,3);
        end
    end
    conn_chain_target_Go = conn_chain_go(reg_chain_go(:,1) == targetregion1 & reg_chain_go(:,2) == targetregion2,:);
    conn_chain_target_Nogo = conn_chain_nogo(reg_chain_nogo(:,1) == targetregion1 & reg_chain_nogo(:,2) == targetregion2,:);
    % FC pairs in Go and Nogo trials
    FCPairID{bin} = vertcat({conn_chain_target_Go},{conn_chain_target_Nogo});
    % unique FC pairs
    UniqueFCPairID{bin} = unique(vertcat(conn_chain_target_Go(:,1:2),conn_chain_target_Nogo(:,1:2)),'rows');
end

%% FCSPs rate in Go and NoGo trials for FC pairs
FCEventRate = CalculateFunctionalCouplingEventRate(FCPairID,SortedUnitTrialMark,SortedSpikeTimeIndividualTrial);

%% Calculate auROC of FCSP events of each neuron pair for Go and NoGo trials
UniqueFCPairAUC = CalculatePopulationFcAUC_WithNonFcspAUC(BinRange,UniqueFCPairID,SortedUnitTrialMark,SortedSpikeTimeIndividualTrial,ResamplingTimes);
save(sprintf('FC Rate and Coding_%s_%s',Target,Group),'FCEventRate','UniqueFCPairAUC','-v7.3');


