
%% Index

% lick                  0;
% Automaticfeeding      94;
% trial start           58;
% laser                 65;
% odor A                9;
% odor B                10;
% response odor         66;
% reward period         8;
% hit                   7;
% miss                  6;
% false                 4;
% correct rejection     5;
% session               61;
% train                 62;
% ITI start             59;

clear; clc; close all;

%% Read data
JAVAFiles = dir('*.ser');
FileNum = size(JAVAFiles,1);

%% Assignment
IsByBlock = 1; % 1: Block; 0: Window
LickThres = 31; MaxTrialsNum = 160; MaxWindowNum = 5; HitRateThres = 0.4;
Beforetrial = 2; DRTSampleDura = 1; DelayDura = 10; DRTResponseDura = 0.5; RewardDura = 1; Aftertrial = 3.0;
TrialsNuminBlock = 20;

%% Behavioral measurement
for num = 1:FileNum
    GoOdorMarker1 = [9 83]; GoOdorMarker2 = [1 3 5];% Go odor marker
    NogoOdorMarker1 = [10 84]; NogoOdorMarker2 = [2 4 6]; % NoGo odor marker
    ResponseOdorMarker1 = 66; ResponseOdorMarker2 = [3 4 5]; % response odor marker
    if FileNum > 1 Filename{1,num} = JAVAFiles(num,1).name; Data = double(ser2mat(Filename{1,num}));
    else Filename = JAVAFiles.name; Data = double(ser2mat(Filename));
    end
    CD = cd;
    % extract information about the Go/NoGo odor
    GoOdor = find(ismember(Data(:,3),GoOdorMarker1) & ismember(Data(:,4),GoOdorMarker2)); GoOdorStartStamp = Data(GoOdor,1);
    NogoOdor = find(ismember(Data(:,3),NogoOdorMarker1) & ismember(Data(:,4),NogoOdorMarker2)); NogoOdorStartStamp = Data(NogoOdor,1);
    SampleOdor = Data(sortrows(cat(1,GoOdor,NogoOdor)),3);
    SampleOdor(ismember(SampleOdor,GoOdorMarker1)) = 1;
    SampleOdor(ismember(SampleOdor,NogoOdorMarker1)) = 2;
    SampleOdorStartStamp = sortrows(cat(1,GoOdorStartStamp,NogoOdorStartStamp));
    % extract information about the response odor
    ResponseOdor = find(ismember(Data(:,3),ResponseOdorMarker1) & ismember(Data(:,4),ResponseOdorMarker2)); ResponseOdorStartStamp = Data(ResponseOdor,1);
    ResponseOdor = Data(ResponseOdor,3);
    ResponseOdor(ResponseOdor == ResponseOdorMarker1) = 3;
    % convert hit, miss, FA, CR into 1, 2, 3, and 4
    DRTOutcome = Data(Data(:,3) == 7 | Data(:,3) == 6 | Data(:,3) == 4 | Data(:,3) == 5,3);
    DRTOutcome(DRTOutcome(:,1) == 7) = 1;
    DRTOutcome(DRTOutcome(:,1) == 6) = 2;
    DRTOutcome(DRTOutcome(:,1) == 4) = 3;
    DRTOutcome(DRTOutcome(:,1) == 5) = 4;
    DRTOutcome1 = DRTOutcome;
    % extract licking time
    LickTime = Data(Data(:,3) == 0,1);
    DiffLick = diff(LickTime);
    RealLickID = [0 (find(DiffLick > LickThres))']+1;
    LickTime = LickTime(RealLickID,1);
    % start time
    OriginTime = min([min(LickTime) min(SampleOdorStartStamp)]);
    % new time
    LickTime = (LickTime - OriginTime)/1000;
    SampleOdorStartStamp = (SampleOdorStartStamp - OriginTime)/1000;
    ResponseOdorStartStamp = (ResponseOdorStartStamp - OriginTime)/1000;
    % Calculate hit, miss, FA, and CR rates
    if IsByBlock == 0
        % plot performance in each window
        TrialsNuminWindow = fix(MaxTrialsNum/MaxWindowNum);
        for i=1:floor(size(DRTOutcome1,1)/TrialsNuminWindow)
            HitNuminWindow(i) = length(find(DRTOutcome1((TrialsNuminWindow*(i-1)+1:TrialsNuminWindow*i),1) == 1));
            MissNuminWindow(i) = length(find(DRTOutcome1((TrialsNuminWindow*(i-1)+1:TrialsNuminWindow*i),1) == 2));
            FalseNuminWindow(i) = length(find(DRTOutcome1((TrialsNuminWindow*(i-1)+1:TrialsNuminWindow*i),1) == 3));
            CRNuminWindow(i) = length(find(DRTOutcome1((TrialsNuminWindow*(i-1)+1:TrialsNuminWindow*i),1) == 4));
            HitRateinWindow(i) = HitNuminWindow(i)/(HitNuminWindow(i)+MissNuminWindow(i));
            CRRateinWindow(i) =CRNuminWindow(i)/(FalseNuminWindow(i)+CRNuminWindow(i));
            PerfinWindow(i) = (HitNuminWindow(i)+CRNuminWindow(i))/TrialsNuminWindow;
        end
        if (floor(size(DRTOutcome1,1)/TrialsNuminWindow) ~= size(DRTOutcome1,1)/TrialsNuminWindow) & (size(DRTOutcome1,1)-TrialsNuminWindow*floor(size(DRTOutcome1,1)/TrialsNuminWindow) >= TrialsNuminWindow/2)
            HitNuminWindow(floor(size(DRTOutcome1,1)/TrialsNuminWindow)+1) = length(find(DRTOutcome1(TrialsNuminWindow*floor(size(DRTOutcome1,1)/TrialsNuminWindow)+1:end,1) == 1));
            MissNuminWindow(floor(size(DRTOutcome1,1)/TrialsNuminWindow)+1) = length(find(DRTOutcome1(TrialsNuminWindow*floor(size(DRTOutcome1,1)/TrialsNuminWindow)+1:end,1) == 2));
            FalseNuminWindow(floor(size(DRTOutcome1,1)/TrialsNuminWindow)+1) = length(find(DRTOutcome1(TrialsNuminWindow*floor(size(DRTOutcome1,1)/TrialsNuminWindow)+1:end,1) == 3));
            CRNuminWindow(floor(size(DRTOutcome1,1)/TrialsNuminWindow)+1) = length(find(DRTOutcome1(TrialsNuminWindow*floor(size(DRTOutcome1,1)/TrialsNuminWindow)+1:end,1) == 4));
            HitRateinWindow(floor(size(DRTOutcome1,1)/TrialsNuminWindow)+1) =  HitNuminWindow(floor(size(DRTOutcome1,1)/TrialsNuminWindow)+1)/( HitNuminWindow(floor(size(DRTOutcome1,1)/TrialsNuminWindow)+1)+ MissNuminWindow(floor(size(DRTOutcome1,1)/TrialsNuminWindow)+1));
            CRRateinWindow(floor(size(DRTOutcome1,1)/TrialsNuminWindow)+1) =  CRNuminWindow(floor(size(DRTOutcome1,1)/TrialsNuminWindow)+1)/( FalseNuminWindow(floor(size(DRTOutcome1,1)/TrialsNuminWindow)+1)+ CRNuminWindow(floor(size(DRTOutcome1,1)/TrialsNuminWindow)+1));
            PerfinWindow(floor(size(DRTOutcome1,1)/TrialsNuminWindow)+1) =  (HitNuminWindow(floor(size(DRTOutcome1,1)/TrialsNuminWindow)+1)+CRNuminWindow(floor(size(DRTOutcome1,1)/TrialsNuminWindow)+1))/length(DRTOutcome1(TrialsNuminWindow*floor(size(DRTOutcome1,1)/TrialsNuminWindow)+1:end,1));
        end
        BelowThresWindowID = find(HitRateinWindow <= HitRateThres); 
        if ~isempty(BelowThresWindowID) 
            HitRateinWindow(:,BelowThresWindowID) = [];
            CRRateinWindow(:,BelowThresWindowID) = [];
            PerfinWindow(:,BelowThresWindowID) = [];
            for i=1:size(BelowThresWindowID,2)
                if BelowThresWindowID(1,i) < length(HitNuminWindow)
                    DRTOutcome1(TrialsNuminWindow*(BelowThresWindowID(1,i)-1)+1:TrialsNuminWindow*BelowThresWindowID(1,i),1)=0;
                    ResponseOdorStartStamp(TrialsNuminWindow*(BelowThresWindowID(1,i)-1)+1:TrialsNuminWindow*BelowThresWindowID(1,i),1) = 0;
                else
                    DRTOutcome1(TrialsNuminWindow*(BelowThresWindowID(1,i)-1)+1:end,1)=0;
                    ResponseOdorStartStamp(TrialsNuminWindow*(BelowThresWindowID(1,i)-1)+1:end,1) = 0;
                end
            end
            DRTOutcome1(DRTOutcome1 == 0) = [];
            ResponseOdorStartStamp(ResponseOdorStartStamp == 0) = [];
        end
        HitNuminSession = length(find(DRTOutcome1 == 1));
        MissNuminSession = length(find(DRTOutcome1 == 2));
        FalseNuminSession = length(find(DRTOutcome1 == 3));
        CRNuminSession = length(find(DRTOutcome1 == 4));
        HitRateinSession = HitNuminSession / (HitNuminSession + MissNuminSession);
        CRRateinSession =  CRNuminSession / (FalseNuminSession + CRNuminSession);
        PerfinSession = (HitNuminSession + CRNuminSession) / size(DRTOutcome1,1);
        AbortedTrialsRatio = (length(SampleOdor) - length(ResponseOdor)) / length(SampleOdor);
    else
        for i=1:floor(size(DRTOutcome1,1)/TrialsNuminBlock)
            HitNuminBlock(i) = length(find(DRTOutcome1((TrialsNuminBlock*(i-1)+1:TrialsNuminBlock*i),1) == 1));
            MissNuminBlock(i) = length(find(DRTOutcome1((TrialsNuminBlock*(i-1)+1:TrialsNuminBlock*i),1) == 2));
            FalseNuminBlock(i) = length(find(DRTOutcome1((TrialsNuminBlock*(i-1)+1:TrialsNuminBlock*i),1) == 3));
            CRNuminBlock(i) = length(find(DRTOutcome1((TrialsNuminBlock*(i-1)+1:TrialsNuminBlock*i),1) == 4));
            HitRateinBlock(i) = HitNuminBlock(i) / (HitNuminBlock(i) + MissNuminBlock(i));
            CRRateinBlock(i) = CRNuminBlock(i) / (FalseNuminBlock(i) + CRNuminBlock(i));
            PerfinBlock(i) = (HitNuminBlock(i) + CRNuminBlock(i)) / TrialsNuminBlock;
        end
        if (floor(size(DRTOutcome1,1)/TrialsNuminBlock) ~= size(DRTOutcome1,1)/TrialsNuminBlock) & (size(DRTOutcome1,1)-TrialsNuminBlock*floor(size(DRTOutcome1,1)/TrialsNuminBlock) >= TrialsNuminBlock/2)
            HitNuminBlock(floor(size(DRTOutcome1,1)/TrialsNuminBlock)+1) = length(find(DRTOutcome1(TrialsNuminBlock*floor(size(DRTOutcome1,1)/TrialsNuminBlock)+1:end,1) == 1));
            MissNuminBlock(floor(size(DRTOutcome1,1)/TrialsNuminBlock)+1) = length(find(DRTOutcome1(TrialsNuminBlock*floor(size(DRTOutcome1,1)/TrialsNuminBlock)+1:end,1) == 2));
            FalseNuminBlock(floor(size(DRTOutcome1,1)/TrialsNuminBlock)+1) = length(find(DRTOutcome1(TrialsNuminBlock*floor(size(DRTOutcome1,1)/TrialsNuminBlock)+1:end,1) == 3));
            CRNuminBlock(floor(size(DRTOutcome1,1)/TrialsNuminBlock)+1) = length(find(DRTOutcome1(TrialsNuminBlock*floor(size(DRTOutcome1,1)/TrialsNuminBlock)+1:end,1) == 4));
            HitRateinBlock(floor(size(DRTOutcome1,1)/TrialsNuminBlock)+1) =  HitNuminBlock(floor(size(DRTOutcome1,1)/TrialsNuminBlock)+1)/( HitNuminBlock(floor(size(DRTOutcome1,1)/TrialsNuminBlock)+1)+ MissNuminBlock(floor(size(DRTOutcome1,1)/TrialsNuminBlock)+1));
            CRRateinBlock(floor(size(DRTOutcome1,1)/TrialsNuminBlock)+1) =  CRNuminBlock(floor(size(DRTOutcome1,1)/TrialsNuminBlock)+1)/( FalseNuminBlock(floor(size(DRTOutcome1,1)/TrialsNuminBlock)+1)+ CRNuminBlock(floor(size(DRTOutcome1,1)/TrialsNuminBlock)+1));
            PerfinBlock(floor(size(DRTOutcome1,1)/TrialsNuminBlock)+1) =  ( HitNuminBlock(floor(size(DRTOutcome1,1)/TrialsNuminBlock)+1)+CRNuminBlock(floor(size(DRTOutcome1,1)/TrialsNuminBlock)+1))/length(DRTOutcome1(TrialsNuminBlock*floor(size(DRTOutcome1,1)/TrialsNuminBlock)+1:end,1));
        end
        if rem(size(HitRateinBlock,2),2) == 1
            HitRateinBlock(:,end) = [];
            CRRateinBlock(:,end) = [];
            PerfinBlock(:,end) = [];
            DRTOutcome1(TrialsNuminBlock*(size(HitRateinBlock,2)-1)+1:end,:) = [];
            ResponseOdorStartStamp(TrialsNuminBlock*(size(HitRateinBlock,2)-1)+1:end,:) = [];
        end
        BelowThresWindowID = find(HitRateinBlock <= 0.5); 
        for i = 1:length(BelowThresWindowID)
            if rem(BelowThresWindowID(i),2) == 1
                HitRateinBlock(1,BelowThresWindowID(i):BelowThresWindowID(i)+1) = -1;
                CRRateinBlock(1,BelowThresWindowID(i):BelowThresWindowID(i)+1) = -1;
                PerfinBlock(1,BelowThresWindowID(i):BelowThresWindowID(i)+1) = -1;
            else
                HitRateinBlock(1,BelowThresWindowID(i)-1:BelowThresWindowID(i)) = -1;
                CRRateinBlock(1,BelowThresWindowID(i)-1:BelowThresWindowID(i)) = -1;
                PerfinBlock(1,BelowThresWindowID(i)-1:BelowThresWindowID(i)) = -1;
            end
        end
        MissedBlockID = find(HitRateinBlock == -1);
        if ~isempty(MissedBlockID)
            for i=1:length(MissedBlockID)
                if MissedBlockID(1,i)~=length(HitRateinBlock) DRTOutcome1(TrialsNuminBlock*(MissedBlockID(1,i)-1)+1:TrialsNuminBlock*MissedBlockID(1,i),1)=0;
                else DRTOutcome1(TrialsNuminBlock*(MissedBlockID(1,i)-1)+1:end,1)=0;
                end
            end
            DRTOutcome1(DRTOutcome1 == 0) = [];
            ResponseOdorStartStamp(DRTOutcome1 == 0) = [];
            PerfinBlock(:,PerfinBlock == -1) = [];
            HitRateinBlock(:,HitRateinBlock == -1) = [];
            CRRateinBlock(:,CRRateinBlock == -1) = [];
        end
        % aborted trials 
        TrialStr = []; TrialStr(:,1) = SampleOdorStartStamp; TrialStr(:,2:3) = zeros(length(SampleOdorStartStamp),2); FinishedTrialID = 0;
        for iSample = 1:length(SampleOdorStartStamp)
            temp = find(ResponseOdorStartStamp < SampleOdorStartStamp(iSample) + DRTSampleDura + DelayDura + 0.5 & ResponseOdorStartStamp > SampleOdorStartStamp(iSample) + DRTSampleDura + DelayDura - 0.5);
            if ~isempty(temp) % defining borders of block further
                TrialStr(iSample,2) = 1; FinishedTrialID = FinishedTrialID + 1;
                if rem(FinishedTrialID,TrialsNuminBlock) ~= 0
                    TrialStr(iSample,3) = fix(FinishedTrialID / TrialsNuminBlock) + 1;
                else
                    TrialStr(iSample,3) = fix(FinishedTrialID / TrialsNuminBlock);
                end
            end
        end
        BorderTrialID = [];
        for iWindow = 1:max(TrialStr(:,3))
            if iWindow > 1
                TrialStr(max(find(TrialStr(:,3) == iWindow - 1)) + 1:max(find(TrialStr(:,3) == iWindow)),4) = iWindow;
            else
                TrialStr(1:max(find(TrialStr(:,3) == iWindow)),4) = iWindow;
            end
        end
        AbortedTrialsRatio_laseroff = nnz(TrialStr(:,2) == 0 & rem(TrialStr(:,4),2) == 1) / nnz(rem(TrialStr(:,4),2) == 1);
        AbortedTrialsRatio_laseron = nnz(TrialStr(:,2) == 0 & rem(TrialStr(:,4),2) == 0) / nnz(rem(TrialStr(:,4),2) == 0);
    end
    % put lick into categorization by trials. 
    LickPSTH = []; Bin = 0.1; LickinTrial = []; TrialID = 1;
    for i = 1:size(ResponseOdorStartStamp,1) % trial
        LickNuminBin = [];
        for t = ResponseOdorStartStamp(i,1) - DelayDura - DRTSampleDura - Beforetrial:Bin:(ResponseOdorStartStamp(i,1) + DRTResponseDura  + RewardDura + Aftertrial - Bin)
            tempLickNumofBin = size(find((LickTime > t) & (LickTime <= t + Bin )),1);
            LickNuminBin = [LickNuminBin tempLickNumofBin];
        end
        LickPSTH = [LickPSTH; LickNuminBin];
        LickinTrial{TrialID} =  LickTime((LickTime > ResponseOdorStartStamp(i,1) - DelayDura - DRTSampleDura - Beforetrial) & (LickTime <= ResponseOdorStartStamp(i,1) + DRTResponseDura  + RewardDura + Aftertrial))-(ResponseOdorStartStamp(i,1) - DelayDura - DRTSampleDura - Beforetrial);
        TrialID = TrialID + 1;
    end
    LickPSTH_Hit = mean(LickPSTH(DRTOutcome1 == 1,:),1); % Lick in Hit trials
    LickPSTH_False = mean(LickPSTH(DRTOutcome1 == 3,:),1); % Lick in False trials
%     %% Reaction time
%     Reactiontime = []; a = 1;
%     for i = 1:size(Lickintrial_Hit,2) % trial
%         for j = 1:size(Lickintrial_Hit{1,i},1)
%             if Lickintrial_Hit{1,i}(j,1) > beforetrial + cueodorlen + delay & length(find(Lickintrial_Hit{1,i}<=Lickintrial_Hit{1,i}(j,1)+0.5 & Lickintrial_Hit{1,i}>= Lickintrial_Hit{1,i}(j,1))) >=3
%                 tempReactiontimepoint = Lickintrial_Hit{1,i}(j,1) - (beforetrial + cueodorlen + delay);
%                 Reactiontime = [Reactiontime tempReactiontimepoint];
%                 break;
%             end
%         end
%         if rem(i,Block/2) == 0 | i == size(Lickintrial_Hit,2)
%             BlockReactiontime{1,a} = Reactiontime;
%             a = a + 1;
%             Reactiontime = [];
%         end
%     end
%     BlockReactiontime = cellfun(@mean, BlockReactiontime, 'UniformOutput', false);
%     BlockReactiontime = (vertcat(BlockReactiontime{:}))';

    %% Save  variable into the mat file.
    if IsByBlock == 0 save([num2str(Filename{1,num}(end-8:end-4))],'HitRateinSession','CRRateinSession','PerfinSession','AbortedTrialsRatio','DRTOutcome1','LickinTrial','LickPSTH','LickPSTH_Hit','LickPSTH_False');
    else
        if FileNum>1 save(Filename{1,num}(end-8:end-4),'HitRateinBlock','CRRateinBlock','PerfinBlock','DRTOutcome1','LickinTrial','AbortedTrialsRatio_laseroff','AbortedTrialsRatio_laseron');
        else save(Filename(end-8:end-4),'HitRateinBlock','CRRateinBlock','PerfinBlock','DRTOutcome1','LickinTrial','AbortedTrialsRatio_laseroff','AbortedTrialsRatio_laseron');
        end
    end
    clearvars -except FileNum IsByBlock LickThres MaxTrialsNum MaxWindowNum HitRateThres Beforetrial DRTSampleDura DelayDura DRTResponseDura RewardDura Aftertrial TrialsNuminBlock 
    clc;
    JAVAFiles = dir('*.ser');
end

%% Summarize all-session results of one mouse
% load mat file
File = dir('*.mat'); FileNum = size(File,1);
for num = 1:size(File,1)
    if FileNum>1 Filename{1,num} = File(num,1).name; Dataofmouse{1,num} = load(Filename{1,num});
    else Filename = File.name; Dataofmouse = load(Filename);
    end
end
CD = cd;

%% Extract results
if IsByBlock == 0
    AbortedTrialsRatio = [];
    HitRate = [];
    CRRate = [];
    Perf = [];
end
for i = 1:size(Dataofmouse,2) % session
    if IsByBlock == 1
        if FileNum > 1
            HitRate{i} = Dataofmouse{1,i}.HitRateinBlock;
            CRRate{i} = Dataofmouse{1,i}.CRRateinBlock;
            Perf{i} = Dataofmouse{1,i}.PerfinBlock;
            DRTTrialMarker{i} = Dataofmouse{1,i}.DRTOutcome1;
            LickinTrial{i} = Dataofmouse{1,i}.LickinTrial;
            AbortedTrialsRatio_laseroff(1,i) = Dataofmouse{1,i}.AbortedTrialsRatio_laseroff;
            AbortedTrialsRatio_laseron(1,i) = Dataofmouse{1,i}.AbortedTrialsRatio_laseron;
        else
            HitRate = Dataofmouse.HitRateinBlock;
            CRRate = Dataofmouse.CRRateinBlock;
            Perf = Dataofmouse.PerfinBlock;
            DRTTrialMarker = Dataofmouse.DRTOutcome1;
            LickinTrial = Dataofmouse.LickinTrial;
            AbortedTrialsRatio_laseroff = Dataofmouse.AbortedTrialsRatio_laseroff;
            AbortedTrialsRatio_laseron = Dataofmouse.AbortedTrialsRatio_laseron;
        end
    else
        DRTTrialMarker{i} = Dataofmouse{1,i}.DRTOutcome1;
        Perf = [Perf Dataofmouse{1,i}.PerfinSession];
        HitRate = [HitRate Dataofmouse{1,i}.HitRateinSession];
        CRRate = [CRRate Dataofmouse{1,i}.CRRateinSession];
        AbortedTrialsRatio = [AbortedTrialsRatio Dataofmouse{1,i}.AbortedTrialsRatio];
        LickinTrial{i} = Dataofmouse{1,i}.LickinTrial;
        LickPSTH{i} = Dataofmouse{1,i}.LickPSTH;
        LickPSTH_Hit{i} = Dataofmouse{1,i}.LickPSTH_Hit;
        LickPSTH_False{i} = Dataofmouse{1,i}.LickPSTH_False;
    end
end
if IsByBlock == 0
    save([CD(end-4:end) '_All_day_Performance'],'Perf','HitRate','CRRate', 'LickinTrial','LickPSTH','LickPSTH_Hit','LickPSTH_False','DRTTrialMarker','AbortedTrialsRatio');
else
    save([CD(end-4:end) '_All_day_Performance'],'Perf','HitRate','CRRate','DRTTrialMarker','LickinTrial','AbortedTrialsRatio_laseroff','AbortedTrialsRatio_laseron');
end
movefile([CD(end-4:end) '_All_day_Performance.mat'],[CD(1:end-6)]);
































