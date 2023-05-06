%% Collect information of each tetrode, to produce one file for one session of one mouse

clear; clc; close all;

%% Target pathway
tic
CurrPath = uigetdir;
AllPath = genpath(CurrPath);
SplitPath = strsplit(AllPath,';');
SubPath = SplitPath';
SubPath = SubPath(2:end-1);

%% Remove first-class pathway
FirstGradeFolder = dir(CurrPath);
ToRemoveFileID = [];
for i = 3:size(FirstGradeFolder,1)
    % pathway name to remove
    ToDeleteFolderName = strcat(CurrPath,'\',FirstGradeFolder(i).name);
    % pathway ID to remove
    for j = 1:length(SubPath)
        if strcmp(ToDeleteFolderName,SubPath(j)) == 1
            ToRemoveFileID = [ToRemoveFileID; j];
        end
    end
end
% remove pathway unrelated to data
SubPath(ToRemoveFileID)= [];

%% Merge cross-channel data, and detect events timestamp
for iPath = 1:size(SubPath,1)
    % go into directory
    Path = SubPath{iPath,1};
    cd(Path);
    % extracellular file
    MatlabFile = dir('*.mat');
    % pl2 file to detect event timestamp
    Pl2File = dir('*.pl2');
    % merge cross-channel data
    Data = [];
    for i = 1:size(MatlabFile,1)
        % load extracellular file
        if size(MatlabFile,1) == 1
            load(MatlabFile.name);
        else
            load(MatlabFile(i,1).name);
        end
        if exist('adc045') == 1
            for x = 1:4:9
                if exist(strcat('adc00',num2str(x)))
                    eval(['Data=[Data;adc00',num2str(x),'];']);
                end
            end
            for x = 13:4:97
                if exist(strcat('adc0',num2str(x)))
                    eval(['Data=[Data;adc0',num2str(x),'];']);
                end
            end
            for x = 101:4:125
                if exist(strcat('adc',num2str(x)))
                    eval(['Data=[Data;adc',num2str(x),'];']);
                end
            end
        else
            for x = 65:4:97
                if exist(strcat('adc0',num2str(x)))
                    eval(['Data=[Data;adc0',num2str(x),'];']);
                end
            end
            for x = 101:4:125
                if exist(strcat('adc',num2str(x)))
                    eval(['Data=[Data;adc',num2str(x),'];']);
                end
            end
        end
    end
    switch exist('adc045')
        case 0
            pl2name = Pl2File.name;
            data = Data;
            Cue1 = PL2EventTs(pl2name,'EVT02'); % S1
            Cue2 = PL2EventTs(pl2name,'EVT03'); % S2 or R
            Lick = PL2EventTs(pl2name,'EVT01'); % lick
            Laser = PL2EventTs(pl2name,'EVT04'); % laser
            Laser = Laser.Ts;
            CorrectLaser = [];
            LaserTrialNum = 80;
            % exclude special case, such as suddenly unknow laser signals
            if ~isempty(Laser) && length(Laser) > 10 
                for j = 1:length(Laser) - 1
                    if isempty(CorrectLaser) == 1
                        if Laser(j+1)-Laser(j) < 14.0423+0.005 && Laser(j+1)-Laser(j) > 14.04231-0.005 % resting-state laser
                            CorrectLaser = [CorrectLaser; Laser(j)];
                        end
                    elseif ~isempty(CorrectLaser)
                        if Laser(j)-CorrectLaser(end) < 14.0423+0.005 && Laser(j)-CorrectLaser(end) > 14.0423-0.005
                            CorrectLaser = [CorrectLaser; Laser(j)];
                        elseif Laser(j+1)-Laser(j) < 14.0423+0.005 && Laser(j+1)-Laser(j) > 14.0423-0.005
                            CorrectLaser = [CorrectLaser; Laser(j)];
                        end
                    end
                end
                if Laser(end)-CorrectLaser(end) < 14.0423+0.005 && Laser(end)-CorrectLaser(end) > 14.0423-0.005
                    CorrectLaser = [CorrectLaser; Laser(end)];
                end
            end
            Laser = CorrectLaser;
            if length(Laser) <= LaserTrialNum || isempty(Laser)
                save(['CIO_' pl2name(1:end-4)],'data','Cue1','Cue2','Lick','Laser','-v7.3');
                AlignRecordingData(strcat('CIO_',pl2name(1:end-4)),LaserTrialNum);
                movefile(['short_' pl2name(1:end-4) '.mat'],Path(1:end-9));
            end
            disp(strcat('///Finish organizing data from_',pl2name(1:end-4),'///'));
            clearvars -except MatlabFile Pl2File iPath SubPath
        case 1
            pl2name = Pl2File.name;
            data = Data(Data(:,1) <= 48,:);
%             pl2 = PL2GetFileIndex(pl2name);
            Cue1 = PL2EventTs(pl2name,'EVT07'); % S1, 3, or 5
            Cue2 = PL2EventTs(pl2name,'EVT08'); % S2, 4, or 6
            Lick = PL2EventTs(pl2name,'EVT06'); % lick
            Laser = PL2EventTs(pl2name,'EVT09'); % laser
            Laser = Laser.Ts;
            CorrectLaser = [];
            for j = 1:length(Laser)-1
                if isempty(CorrectLaser) == 1
                    if (Laser(j+1)-Laser(j) < 0.0701+0.005 && Laser(j+1)-Laser(j)>0.0701-0.005) || (Laser(j+1)-Laser(j) < 10.086+0.005 && Laser(j+1)-Laser(j) > 10.086-0.005)
                        CorrectLaser = [CorrectLaser; Laser(j)];
                    end
                elseif ~isempty(CorrectLaser)
                    if (Laser(j)-CorrectLaser(end) < 0.0701+0.005 && Laser(j)-CorrectLaser(end) > 0.0701-0.005) || (Laser(j)-CorrectLaser(end) < 10.086+0.005 && Laser(j)-CorrectLaser(end) > 10.086-0.005)
                        CorrectLaser = [CorrectLaser; Laser(j)];
                    elseif Laser(j+1)-Laser(j) < 0.0701+0.005 && Laser(j+1)-Laser(j) > 0.0701-0.005
                        CorrectLaser = [CorrectLaser; Laser(j)];
                    end
                end
            end
            if (Laser(end)-CorrectLaser(end) < 0.0701+0.005 && Laser(end)-CorrectLaser(end) > 0.0701-0.005) || (Laser(end)-CorrectLaser(end) < 10.086+0.005 && Laser(end)-CorrectLaser(end) > 10.086-0.005)
                CorrectLaser = [CorrectLaser; Laser(end)];
            end
            if length(CorrectLaser) == 800
                Laser = CorrectLaser(length(CorrectLaser)/2+1:end,1);
                if ~isempty(data)
                    save(['CIO_' Pl2File.name(4:6) Pl2File.name(8:end-4)],'data','Cue1','Cue2','Lick','Laser','-v7.3');
                    AlignRecordingData(strcat('CIO_',Pl2File.name(4:6),Pl2File.name(8:end-4)),400);
                    movefile(['short_' Pl2File.name(4:6) Pl2File.name(8:end-4) '.mat'],Path(1:end-9));
                end
            end
            disp(strcat('///Finish organizing data from_',Pl2File.name(4:6),Pl2File.name(8:end-4),'///'));
            data = Data(Data(:,1) > 48,:);
            Cue1 = PL2EventTs(pl2name,'EVT07'); % S1, 3, or 5
            Cue2 = PL2EventTs(pl2name,'EVT08'); % S2, 4, or 6
            Lick = PL2EventTs(pl2name,'EVT06'); % lick
            if length(CorrectLaser)==800
                Laser = CorrectLaser(1:length(CorrectLaser)/2,1);
                if ~isempty(data)
                    save(['CIO_' Pl2File.name(1:3) Pl2File.name(8:end-4)],'data','Cue1','Cue2','Lick','Laser','-v7.3');
                    AlignRecordingData(strcat('CIO_',Pl2File.name(1:3),Pl2File.name(8:end-4)),400);
                    movefile(['short_' Pl2File.name(1:3) Pl2File.name(8:end-4) '.mat'],Path(1:end-9));
                end
            end
            disp(strcat('///Finish organizing data from_',Pl2File.name(1:3),Pl2File.name(8:end-4),'///'));
            clearvars -except MatlabFile Pl2File iPath SubPath
    end
end
toc
disp(['Run time: ',num2str(toc)]);
% %% Intan Recording
%     Lick = Event{1,1}; % Lick
%     Cue1 = Event{1,2}; % Sample or Test Odor 1 in DPA
%     Cue2 = Event{1,3}; % Sample or Test Odor 2 in DPA or response odor in ODR
%     Laser = Event{1,6}; % Laser



