%% This code adds ID of distance from the optical fiber and ID of putative GABAergic fast spiking (FS) neurons or putative pyramidal neurons.

clear; clc; close all;

%% Enter path
CurrPath = uigetdir;

%% Mat file
MatFile = dir('*.mat');

%% Distance ID and neuron type ID
% //////distance ID////// 0: 0 mm; 1: 1mm; 2: 2 mm; 3: 3 mm
% //////neuron type ID////// 0: unclassified; 1: pyramidal neurons; 2: FS neurons
for iFile = 1:size(MatFile,1)
    fprintf('Analyzing %dth file of total %d files\n',iFile,size(MatFile,1));
    if size(MatFile,1) == 1
        load(MatFile.name);
    else
        load(MatFile(iFile,1).name);
    end
    %% Distance ID
    DistID = zeros(size(SingleUnitList,1),1);
    DistID(SingleUnitList(:,1)>=9 & SingleUnitList(:,1)<=16) = 1;
    DistID(SingleUnitList(:,1)>=17 & SingleUnitList(:,1)<=24) = 2;
    DistID(SingleUnitList(:,1)>=25 & SingleUnitList(:,1)<=32) = 3;

    %% Neuron type ////// FS neurons peak to trough < 350 us; pyramidal neurons peak to trough > 450 us
    % *Spatiotemporal constraints on optogenetic inactivation in cortical circuits-Karel Svoboda*
    UnitTypeID = zeros(size(SingleUnitList,1),1);
    PeakTroughInterval = [];
    for iUnit = 1:size(NewWaveForm)
       [tempPeakTroughInterval,~,~,~,~] = CalculatePeakTroughDuration(NewWaveForm(iUnit,:));
       PeakTroughInterval(end+1,1) = tempPeakTroughInterval;
    end
    UnitTypeID(PeakTroughInterval>450) = 1;
    UnitTypeID(PeakTroughInterval<350) = 2;
    save(MatFile(iFile,1).name,'LaserResults','LaserRGResults','Newdata','NewWaveForm','PeakTroughInterval','SingleUnitList','SpikeTime','UnitTypeID','DistID','-v7.3');
end










