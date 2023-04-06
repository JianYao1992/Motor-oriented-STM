
%% Index representation

% lick                  0;
% Automaticfeeding      94;
% trial start           58;
% test tone             90;
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

%% Assignment
TrialsNuminWindow = 20;
LickThres = 31; HitRateThres = 0.5;
Beforetrial = 2; DRTSampleDura = 1; DelayDura = 10; DRTResponseDura = 0.5; RewardDura = 1; Aftertrial = 3.0;

%% Behavioral measurements
% load result
JAVAFiles = dir('*.ser');
FileNum = size(JAVAFiles,1);
for num = 1:FileNum
    if FileNum > 1 Filename{1,num} = JAVAFiles(num,1).name; Data = double(ser2mat(Filename{1,num}));
    else Filename = JAVAFiles.name; Data = double(ser2mat(Filename));
    end
    CD = cd;
    % change hit, miss, FA, and CR into 1, 2, 3, and 4
    GNGOutcome = Data(Data(:,3) == 7 | Data(:,3) == 6 | Data(:,3) == 4 | Data(:,3) == 5,3);
    GNGOutcome(GNGOutcome(:,1) == 7) = 1;
    GNGOutcome(GNGOutcome(:,1) == 6) = 2;
    GNGOutcome(GNGOutcome(:,1) == 4) = 3;
    GNGOutcome(GNGOutcome(:,1) == 5) = 4;
    GNGOutcome1 = GNGOutcome;
    % analysis in each window
    for i=1:floor(size(GNGOutcome1,1)/TrialsNuminWindow)
        HitNuminWindow(i) = length(find(GNGOutcome1((TrialsNuminWindow*(i-1)+1:TrialsNuminWindow*i),1) == 1));
        MissNuminWindow(i) = length(find(GNGOutcome1((TrialsNuminWindow*(i-1)+1:TrialsNuminWindow*i),1) == 2));
        FalseNuminWindow(i) = length(find(GNGOutcome1((TrialsNuminWindow*(i-1)+1:TrialsNuminWindow*i),1) == 3));
        CRNuminWindow(i) = length(find(GNGOutcome1((TrialsNuminWindow*(i-1)+1:TrialsNuminWindow*i),1) == 4));
        WindowHitRate(i) = HitNuminWindow(i)/(HitNuminWindow(i)+MissNuminWindow(i));
        WindowCRRate(i) = CRNuminWindow(i)/(FalseNuminWindow(i)+CRNuminWindow(i));
    end
    if (floor(size(GNGOutcome1,1)/TrialsNuminWindow) ~= size(GNGOutcome1,1)/TrialsNuminWindow) & (size(GNGOutcome1,1)-TrialsNuminWindow*floor(size(GNGOutcome1,1)/TrialsNuminWindow) >= TrialsNuminWindow/2)
        HitNuminWindow(floor(size(GNGOutcome1,1)/TrialsNuminWindow)+1) = length(find(GNGOutcome1(TrialsNuminWindow*floor(size(GNGOutcome1,1)/TrialsNuminWindow)+1:end,1) == 1));
        MissNuminWindow(floor(size(GNGOutcome1,1)/TrialsNuminWindow)+1) = length(find(GNGOutcome1(TrialsNuminWindow*floor(size(GNGOutcome1,1)/TrialsNuminWindow)+1:end,1) == 2));
        FalseNuminWindow(floor(size(GNGOutcome1,1)/TrialsNuminWindow)+1) = length(find(GNGOutcome1(TrialsNuminWindow*floor(size(GNGOutcome1,1)/TrialsNuminWindow)+1:end,1) == 3));
        CRNuminWindow(floor(size(GNGOutcome1,1)/TrialsNuminWindow)+1) = length(find(GNGOutcome1(TrialsNuminWindow*floor(size(GNGOutcome1,1)/TrialsNuminWindow)+1:end,1) == 4));
        WindowHitRate(floor(size(GNGOutcome1,1)/TrialsNuminWindow)+1) =  HitNuminWindow(floor(size(GNGOutcome1,1)/TrialsNuminWindow)+1)/( HitNuminWindow(floor(size(GNGOutcome1,1)/TrialsNuminWindow)+1)+ MissNuminWindow(floor(size(GNGOutcome1,1)/TrialsNuminWindow)+1));
        WindowCRRate(floor(size(GNGOutcome1,1)/TrialsNuminWindow)+1) =  CRNuminWindow(floor(size(GNGOutcome1,1)/TrialsNuminWindow)+1)/( FalseNuminWindow(floor(size(GNGOutcome1,1)/TrialsNuminWindow)+1)+ CRNuminWindow(floor(size(GNGOutcome1,1)/TrialsNuminWindow)+1));
    end
    BelowThresWindowID = find(WindowHitRate <= HitRateThres);
    if ~isempty(BelowThresWindowID)
        WindowHitRate(:,BelowThresWindowID) = [];
        WindowCRRate(:,BelowThresWindowID) = [];
        for i=1:size(BelowThresWindowID,2)
            if BelowThresWindowID(1,i) < length(HitNuminWindow)
                GNGOutcome1(TrialsNuminWindow*(BelowThresWindowID(1,i)-1)+1:TrialsNuminWindow*BelowThresWindowID(1,i),1)=0;
            else
                GNGOutcome1(TrialsNuminWindow*(BelowThresWindowID(1,i)-1)+1:end,1)=0;
            end
        end
        GNGOutcome1(GNGOutcome1 == 0) = [];
    end
    % save  variable into the mat file.
    if FileNum > 1
        save([num2str(Filename{1,num}(end-8:end-4))],'WindowHitRate','WindowCRRate','WindowPerf','DRTOutcome1');
    else
        save([num2str(Filename(end-8:end-4))],'WindowHitRate','WindowCRRate','WindowPerf','DRTOutcome1');
    end
    clearvars -except FileNum LickThres HitRateThres Beforetrial DRTSampleDura DelayDura DRTResponseDura RewardDura Aftertrial TrialsNuminWindow
    clc;
    JAVAFiles = dir('*.ser');
end

%% Summarize all-session results
% load all mat file
File = dir('*.mat'); FileNum = size(File,1);
for num = 1:size(File,1)
    if FileNum>1 Filename{1,num} = File(num,1).name; Dataofmouse{1,num} = load(Filename{1,num});
    else Filename = File.name; Dataofmouse = load(Filename);
    end
end
CD = cd;
% extract results
for i = 1:size(Dataofmouse,2) % Day
    if iscell(Dataofmouse)
        GNGTrialMarker{i} = Dataofmouse{1,i}.GNGOutcome1;
        WindowHitRate{1,i} = Dataofmouse{1,i}.HitRateinWindow;
        WindowCRRate{1,i} = Dataofmouse{1,i}.CRRateinWindow;
    else
        GNGTrialMarker{i} = Dataofmouse.GNGOutcome1;
        WindowHitRate{1,i} = Dataofmouse.WindowHitRate;
        WindowCRRate{1,i} = Dataofmouse.WindowCRRate;
    end
end
save([CD(end-4:end) '_All_day_Performance'],'WindowHitRate','WindowCRRate','GNGTrialMarker');
movefile([CD(end-4:end) '_All_day_Performance.mat'],[CD(1:end-6)]);
































