
%% Cross-correlogram analysis for neuronal pairs of simultaneously recorded neurons

clear; clc; close all;

%% Assignment
Group = 'CtrlGroup'; 
LearningDay = 6; 
beforeneuronnum = 0;
BeforeOdorDuration = 4;
BilateralmPFCMiceID = [{'M19'} {'M20'}]; 

%% Add path of FieldTrip
addpath('D:\YJ\fieldtrip-20200326');
CurrPath = 'H:\20190811DRT_RQP_LaserAndNoLaser\SelectivityAndLaserFiring\Ctrl';
AllPath = genpath(CurrPath);
SplitPath = strsplit(AllPath,';'); 
SubPath = SplitPath';
SubPath = SubPath(2:end-1);

%% Spike rasters of trials for all neurons
MiceID = []; 
mPFC_ID = [];
aAIC_ID = [];
trial = []; 
time = [];
trialinfo = [];
trialtime = [];
LearningDayID = [];
FromBrainArea = []; % 1: mPFC 2: aAIC
for iPath = 1:size(SubPath,1)
    Path = SubPath{iPath,1}; 
    cd(Path);
    JAVAFiles = dir('*.mat');
    for j = 1:size(JAVAFiles,1)
        Filename{1,FileID} = JAVAFiles(j,1).name;
        load(Filename{1,FileID});
        if ~isempty(SingleUnitList)
            RGResults = cellfun(@(x) (x-BeforeOdorDuration)', RGResults, 'UniformOutput',false);
            MiceID = [MiceID; repmat(str2num(Filename{1,FileID}(8:9)),size(SingleUnitList,1),1)];
            LearningDayID = [LearningDayID; repmat(str2num(Filename{1,FileID}(end-4)),size(SingleUnitList,1),1)];
            trialinfo = [trialinfo; repmat({TrialMark(:,2:end)},size(SingleUnitList,1),1)];
            temptrialtime = repmat([-1*BeforeOdorDuration 6.5],size(TrialMark,1),1);
            trialtime = [trialtime; repmat({temptrialtime},size(SingleUnitList,1),1)];
            Position = [];
            for k = 1:length(BilateralmPFCMiceID)
                Position = [Position regexp(Filename{1,FileID},BilateralmPFCMiceID{1,k})];
            end
            % ID of aAIC neurons
            if ~isempty(Position)
                tempaAICID = [];
            else
                tempaAICID = find((SingleUnitList(:,1)>=1 & SingleUnitList(:,1)<=8) | (SingleUnitList(:,1)>=17 & SingleUnitList(:,1)<=24));
            end
            if ~isempty(tempaAICID)
                aAIC_ID = [aAIC_ID tempaAICID'+beforeneuronnum];
                tempaAICRaster = RGResults(tempaAICID,:);
                for iUnit = 1:size(tempaAICRaster,1)
                    temptrial = []; temptime = [];
                    for iTrial = 1:size(tempaAICRaster,2)
                        temptrial = [temptrial iTrial*ones(1,length(tempaAICRaster{iUnit,iTrial}))];
                        temptime = [temptime tempaAICRaster{iUnit,iTrial}];
                    end
                    trial{1,tempaAICID(iUnit)+beforeneuronnum} = temptrial;
                    time{1,tempaAICID(iUnit)+beforeneuronnum} = temptime;
                end
                FromBrainArea = [FromBrainArea; repmat(2,length(tempaAICID),1)];
            end
            % ID of mPFC neurons
            if ~isempty(Position)
                tempmPFCID = find(SingleUnitList(:,1)>=1 & SingleUnitList(:,1)<=16);  % mice with bilateral mPFC recording
            else
                tempmPFCID = find((SingleUnitList(:,1)>=9 & SingleUnitList(:,1)<=16) | (SingleUnitList(:,1)>=25 & SingleUnitList(:,1)<=32));
            end
            if ~isempty(tempmPFCID)
                mPFC_ID = [mPFC_ID tempmPFCID'+beforeneuronnum];
                tempmPFCRaster = RGResults(tempmPFCID,:);
                for iUnit = 1:size(tempmPFCRaster,1)
                    temptrial = []; temptime = [];
                    for iTrial = 1:size(tempmPFCRaster,2)
                        temptrial = [temptrial iTrial*ones(1,length(tempmPFCRaster{iUnit,iTrial}))];
                        temptime = [temptime tempmPFCRaster{iUnit,iTrial}];
                    end
                    trial{1,tempmPFCID(iUnit)+beforeneuronnum} = temptrial;
                    time{1,tempmPFCID(iUnit)+beforeneuronnum} = temptime;
                end
                FromBrainArea = [FromBrainArea; repmat(1,length(tempmPFCID),1)];
            end
            beforeneuronnum = beforeneuronnum + size(SingleUnitList,1);
            disp(num2str(beforeneuronnum ));
        end
        FileID = FileID + 1; 
    end
end
cd(CurrPath);
OriginInfo = [MiceID LearningDayID FromBrainArea];

%% Memory mPFC neurons
load(['mPFCSelectivityData_' Group]); 
SustainedID_mPFC = mPFC_ID(SustainedUnitID); TransientCodingNeuronID = vertcat(TransientCodingNeuronID,SwitchedSustainedUnitID);
TransientID_mPFC = mPFC_ID(TransientCodingNeuronID);
NocodingID_mPFC = setdiff(mPFC_ID, horzcat(SustainedID_mPFC,TransientID_mPFC));

%% Memory aAIC neurons
load(['aAICSelectivityData_' Group]); 
SustainedID_aAIC = aAIC_ID(SustainedUnitID); TransientCodingNeuronID = vertcat(TransientCodingNeuronID,SwitchedSustainedUnitID);
TransientID_aAIC = aAIC_ID(TransientCodingNeuronID);
NocodingID_aAIC = setdiff(aAIC_ID, horzcat(SustainedID_aAIC,TransientID_aAIC));

%% Summarize all sustained and transient neurons
SustainedID_Both = [SustainedID_mPFC SustainedID_aAIC]; 
TransientID_Both = [TransientID_mPFC TransientID_aAIC]; 

%% Basic information of all neurons
load(['UnitsSummary_' Group]);
UnitsSummary = vertcat(UnitsSummary.mPFC, UnitsSummary.aAIC);
UnitsInfo = [];
for i = 1:size(UnitsSummary,1)
    UnitsInfo(UnitsSummary{i,1}(1),:) = UnitsSummary(i,:);
end

%% Distribution of lags between spike trains of neurons in pair
tic;
for iBin = -1:1:5 
	bin_range = [iBin iBin+1]; % target period /// [-1 0]: baseline; [0 1]: sample; [1 2]: 1st second of delay; [2 3]: 2nd second of delay; [3 4]: 3rd second of delay; [4 5]: 4th second of delay; [5 6]: response ///
	FileID = 1;
    Sums = [];
	for i = 1:length(SubPath)
		iMouse = str2num(SubPath{i}(end-1:end));
	    if iMouse >= 10
	    	MiceIDLabel = 'M';
	    else
	    	MiceIDLabel = 'M0';
	    end
	    tempPath = SubPath{i};
	    for iDay = 1:LearningDay % learning day
	        for tempfile = 1:length(Filename)
	            if ~isempty(strfind(Filename{tempfile},strcat(MiceIDLabel,num2str(iMouse)))) && ~isempty(strfind(Filename{tempfile},strcat('Day0',num2str(iDay))))
	                folder = strcat(tempPath,'\',Filename{tempfile});
	                break;
	            end
            end
            NeuronID_Simul = find(OriginInfo(:,1)==iMouse & OriginInfo(:,2)==iDay); % containing both mPFC and aAIC
            SustainedID_Simul = [];
            TransientID_Simul = [];
            NonmemoryID_Simul = [];
            % classify neuron into sustained, transient or non-memory one
            for iUnit = 1:length(NeuronID_Simul)
                if ismember(NeuronID_Simul(iUnit),SustainedID_Both)
                    SustainedID_Simul = [SustainedID_Simul NeuronID_Simul(iUnit)];
                elseif ismember(NeuronID_Simul(iUnit),TransientID_Both)
                    TransientID_Simul = [TransientID_Simul NeuronID_Simul(iUnit)];
                else
                    NonmemoryID_Simul = [NonmemoryID_Simul NeuronID_Simul(iUnit)];
                end
            end
	        if ~isempty(NeuronID_Simul)
	            spikeTrials = [];                
	            for iUnit = 1:length(NeuronID_Simul)
	                spikeTrials.label{iUnit,1} = num2str(NeuronID_Simul(iUnit));
	            end
	            spikeTrials.time = time(:,NeuronID_Simul);
	            spikeTrials.trial = trial(:,NeuronID_Simul);
	            spikeTrials.trialtime = trialtime{NeuronID_Simul(1,1),1};
	            spikeTrials.trialinfo = trialinfo{NeuronID_Simul(1,1),1};
	            [xc_s1,xcshuf_s1,xc_s2,xcshuf_x2] = plotxcorr(spikeTrials,bin_range); 
                for iUnit = 1:length(NeuronID_Simul)
                    xc_s1.label{iUnit,2} = UnitsInfo{NeuronID_Simul(iUnit),1}; % waveform stats
                    xc_s1.label{iUnit,3} = UnitsInfo{NeuronID_Simul(iUnit),2}; % waveform voltage value
                    xc_s1.label{iUnit,4} = UnitsInfo{NeuronID_Simul(iUnit),3}; % sample prefer
                    xc_s1.label{iUnit,5} = UnitsInfo{NeuronID_Simul(iUnit),4}; % FR in hit trials
                    xc_s1.label{iUnit,6} = UnitsInfo{NeuronID_Simul(iUnit),5}; % FR in CR trials
                    xc_s1.label{iUnit,7} = UnitsInfo{NeuronID_Simul(iUnit),6}; % brain region
                end
           		tempsums = {FileID,folder,SustainedID_Simul,TransientID_Simul,NonmemoryID_Simul,xc_s1,xcshuf_s1,xc_s2,xcshuf_x2};
                Sums = [Sums; tempsums];
	            fprintf('%d of %d\n',FileID,length(Filename))
	        end
	        FileID = FileID + 1;
	    end
    end
    save(sprintf(['XCORR_%d_%d_2msbin_Projection_' Group '.mat'],bin_range(1),bin_range(2)),'Sums','-v7.3'); %prefix
end
toc;

return
function [Xc_S1,Xshuff_S1,Xc_S2,Xshuff_S2] = plotxcorr(spikeTrials,bin_range)
% https://www.nature.com/articles/nn799
% A role for inhibition in shaping the temporal flow of information in prefrontal cortex
% Christos Constantinidis, Graham V. Williams & Patricia S. Goldman-Rakic 
% Nature Neuroscience volume 5, pages175-180(2002)
% 
% Neuron, Volume 76
% Functional microcircuit recruited during retrieval of object association memory in monkey perirhinal cortex
% Toshiyuki Hirabayashi, Daigo Takeuchi, Keita Tamura, and Yasushi Miyashita

cfg             = [];
cfg.maxlag      = 0.1; % maximum 100 ms
cfg.binsize     = 0.002; % bins of 2 ms
cfg.outputunit  = 'raw'; 
cfg.latency     = bin_range; % time bin based on sample onset
cfg.vartriallen = 'no'; % allow variable trial lengths
cfg.debias      = 'no';

%% Trials in context 1
cfg.trials      = find(spikeTrials.trialinfo(:,1)==1 & (spikeTrials.trialinfo(:,3)==1|spikeTrials.trialinfo(:,3)==4));
if numel(cfg.trials) < 2
    Xc_S1 = [];
    Xshuff_S1 = [];
else
    cfg.method      = 'xcorr'; % normal cross-correlogram
    Xc_S1 = ft_spike_xcorr(cfg,spikeTrials);
    cfg.method      = 'shiftpredictor'; % shift predictor
    Xshuff_S1 = ft_spike_xcorr(cfg,spikeTrials);
end

%% Trials in context 2
cfg.trials      = find(spikeTrials.trialinfo(:,1)==2 & (spikeTrials.trialinfo(:,3)==1|spikeTrials.trialinfo(:,3)==4));
if numel(cfg.trials) < 2 
    Xc_S2 = [];
    Xshuff_S2 = [];
else
    cfg.method      = 'xcorr'; % normal cross-correlogram
    Xc_S2 = ft_spike_xcorr(cfg,spikeTrials);
    cfg.method      = 'shiftpredictor'; % shift predictor
    Xshuff_S2 = ft_spike_xcorr(cfg,spikeTrials);
end
end



