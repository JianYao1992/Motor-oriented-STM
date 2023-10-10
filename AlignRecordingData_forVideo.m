function AlignRecordingData_forVideo(FileName, LaserTrialNum)

%% Load recording data
load(FileName);

%% Parameter definition
FirstOdorLen = 1;
Delay = 4;
RespOdorLen = 0.5;
Waterlen = 1;
ITI = 10;
TimeGain = 10; % number of bins of 1 sec
WaveformStartPosition = 4; % start column ID for waveform information of recording data

%% Units aligned with tetrodes
TetrodeList = [];
UnitsList = [];
WaveForm = [];
for itr = data(1,1):data(end,1) % 128 channels
    TetrodeList = [TetrodeList; itr*ones(max(data(data(:,1)==itr,2)),1)];
    UnitsList = [UnitsList 1:max(data(data(:,1)==itr,2))];
    for itr1 = 1:max(data(data(:,1)==itr,2))
        WaveForm = [WaveForm; mean(data(data(:,1)==itr & data(:,2)==itr1,WaveformStartPosition:size(data,2)))]; % waveform
    end
end
UnitsList = [TetrodeList UnitsList'];
data(:,1) = data(:,1)-32;
TetrodeList = TetrodeList-32;
UnitsList(:,1) = UnitsList(:,1)-32;

%% Extract timestamps of A and B
A = Cue1.Ts;
B = Cue2.Ts;
if ~isempty(A)
    A(:,2) = 1;
    B(:,2) = 2;
end
rawA = A;
rawB = B;
for itr = 1:length(rawA)
    if isempty(find(rawB(:,1) == rawA(itr,1)));
    else
        B(B(:,1)==rawA(itr,1),2) = 0;
        A(itr,2) = 0;
    end
    if itr < length(rawA) && rawA(itr+1,1) - rawA(itr,1) < 0.5
        A(itr,2) = 0;
    end
end
for itr = 1:length(rawB)
    if itr < length(rawB) && rawB(itr+1,1) - rawB(itr,1) < 0.5
        B(itr,2) = 0;
    end
end
if ~isempty(A)
    A(A(:,2)==0,:) = [];
    B(B(:,2)==0,:) = [];
end

%% All timestamps of sample odor onset
RawTS = sortrows([Cue1.Ts;Cue2.Ts]);
PutativeSampleOnsetTS = [];
for i = 1:size(RawTS,1)-1
    temp = (RawTS(i+1)-RawTS(i)>FirstOdorLen+Delay-0.5 & RawTS(i+1)-RawTS(i)<FirstOdorLen+Delay+0.5) | (RawTS(i+1)-RawTS(i)>FirstOdorLen+Delay+RespOdorLen+Waterlen+ITI-0.5 & RawTS(i+1)-RawTS(i)<FirstOdorLen+Delay+RespOdorLen+Waterlen+ITI+0.5);
    if temp == 1
        PutativeSampleOnsetTS = [PutativeSampleOnsetTS RawTS(i)];
    end
end

%% Units FR
Newdata = data;
FiringRate = zeros(length(UnitsList),TimeGain*ceil(max(Newdata(:,3))));
for itr = 1:length(Newdata)
    if Newdata(itr,2) ~= 0
        FiringRate(Newdata(itr,2)+length(find(TetrodeList<Newdata(itr,1))),ceil(TimeGain*Newdata(itr,3))) = FiringRate(Newdata(itr,2)+length(find(TetrodeList<Newdata(itr,1))),ceil(TimeGain*Newdata(itr,3))) + 1;
    end
end
% units number
for itr = 1:length(Newdata)
    if Newdata(itr,2) ~= 0
        Newdata(itr,2) = Newdata(itr,2) + length(find(TetrodeList<Newdata(itr,1)));
    end
end

%% Filter lick signals
if ~isempty(Lick.Ts)
    NewLick = Lick.Ts(1,:);
    for itr = 2:length(Lick.Ts)
        if Lick.Ts(itr,1)-Lick.Ts(itr-1,1) > 0.031
            NewLick = [NewLick; Lick.Ts(itr,1)];
        end
        if Lick.Ts(itr,1)-NewLick(length(NewLick)) > 0.1
            NewLick = [NewLick; Lick.Ts(itr,1)];
        end
    end
end

%% Checking period
Odor = sortrows([A;B]);
if ~isempty(Odor)
    % abortion ID
    abort = find(diff(Odor(:,1)) > ITI + FirstOdorLen + Delay + RespOdorLen + Waterlen - 0.5 & diff(Odor(:,1)) < ITI + FirstOdorLen + Delay + RespOdorLen + Waterlen + 0.5);
    % remove unreasonable abortion
    if isempty(abort) ~= 1
        for i = 1:size(abort,1)
            if abort(i,1) ~= 1
                a = Odor(abort(i,1),1) - Odor(abort(i,1)-1,1);
                if a < FirstOdorLen + Delay + 0.5 && a > FirstOdorLen + Delay - 0.5
                    abort(i,1) = 0;
                end
            end
        end
        abort(abort(:,1)==0,:) = [];
        Odor(abort,:) = [];
    end
    % filter further, to delete noise
    i = 1;
    while i < size(Odor,1)
        if (i == 1) || (isempty(find(Odor(1:i-1,2)~=0)) && i > 1)
            if Odor(i+1,1) - Odor(i,1) > FirstOdorLen + Delay - 0.5 && Odor(i+1,1) - Odor(i,1) < FirstOdorLen + Delay + 0.5
                i = i + 2;
            else
                Odor(i,2) = 0;
                i = i + 1;
            end
        else
            if (Odor(i,1) - Odor(max(find(Odor(1:i-1,2)~=0)),1) > ITI + RespOdorLen + Waterlen - 0.5) && Odor(i+1,1) - Odor(i,1) > FirstOdorLen + Delay - 0.5 && Odor(i+1,1) - Odor(i,1) < FirstOdorLen + Delay + 0.5
                i = i + 2;
            else
                Odor(i,2) = 0;
                i = i + 1;
            end
        end
    end
    remain = find(Odor(1:end-1,2)~=0);
    if Odor(remain(end,1),1)-Odor(remain(end-1,1),1) > FirstOdorLen+Delay-0.3 && Odor(remain(end,1),1)-Odor(remain(end-1,1),1) < FirstOdorLen+Delay+0.3 % special condition
        Odor(end,2) = 0;
    end
    Odor(Odor(:,2)==0,:) = [];
end
% trial structure
Odor = Odor';
OdorStart = [];
for i = 1:size(Odor,2)
    if rem(i,2) == 0
        OdorStart = [OdorStart Odor(:,i)];
    end
end
if ~isempty(OdorStart)
    CheckStart = OdorStart(1,:) + RespOdorLen;
    CheckEnd = CheckStart + Waterlen;
else
    CheckStart = [];
    CheckEnd = [];
end
% timestamp of laser application
if ~isempty(Laser) && length(Laser) <= LaserTrialNum
    LaserStart = Laser;
end
if ~isempty(Odor)
    SampleOnsetTs = Odor(1,1:2:end);
    TestOnsetTS = Odor(1,2:2:end);
    SampleTestOnsetTs = horzcat(SampleOnsetTs(:),TestOnsetTS(:));
end

%% ID of unaborted trials in all trials
AllTrialsNumber = numel(PutativeSampleOnsetTS);
[~,CompletedTrialID,~] = intersect(PutativeSampleOnsetTS,SampleOnsetTs,'stable');
% %% FR of baseline period
% PseudoLength = 10;
% PseudoTime = 290;
% PseudoTrialRG = [];
% for itr = 1:size(FiringRate,1)
%     PseudoTrialRG = [PseudoTrialRG; [{[]}]];
% end
% for itr = 1:PseudoLength:PseudoTime-PseudoLength+1
%     for itr1 = 1:size(FiringRate,1)
%         PseudoTrialRG{itr1,1} = [PseudoTrialRG{itr1,1}; {Newdata(Newdata(:,2) == itr1 & Newdata(:,3) > itr & Newdata(:,3) < itr+PseudoLength,3)-itr}];
%     end
% end

%% Maximal ID of trials with activity
AfterLaser = 9;
MaxTrialNum = find(CheckEnd*TimeGain > size(FiringRate,2));
if ~isempty(MaxTrialNum)
    MaxTrialNum = min(MaxTrialNum)-1;
else
    MaxTrialNum = size(CheckEnd,2);
end
if exist('LaserStart')
    LaserStart = LaserStart(LaserStart<(size(FiringRate,2)/TimeGain-AfterLaser-0.5));
end

%% Align FR with events
NewWaveForm = [];
SingleUnitList = [];
TrialMark = [];
SpikeTime = [];
Results = [];
RGResults = [];
if exist('LaserStart')
    LaserResults = [];
    LaserRGResults = [];
end
BeforeFirstOdor = 9;
LickTime = [];
LickRate = [];
SingleUnitCount = 0;
UnitID = 1;
DiscardUnitID = [];
for itr = 1:size(FiringRate,1)
    if (~isempty(CheckEnd) && mean(FiringRate(itr,1:round(TimeGain*CheckEnd(1,MaxTrialNum))),2)*TimeGain >= 2) || (isempty(CheckEnd) && mean(FiringRate(itr,:),2)*TimeGain >= 2)
        NewWaveForm = [NewWaveForm; WaveForm(itr,:)];
        SpikeTrace = FiringRate(itr,:); % unit FR
        SingleUnitList = [SingleUnitList; UnitsList(itr,:)];
        SingleUnitCount = SingleUnitCount + 1;
        tempdata = Newdata(Newdata(:,2)==itr,:); % unit spike data
        SP = [];
        RG = [];
        SpikeTime{UnitID,1} = tempdata(:,3); % for further cross correlogram analysis, here not considering S1 and S2, respectively
        UnitID = UnitID + 1;
        for itr1 = 1:MaxTrialNum % all trials
            % sample and response odor in the trial
            SampleOdor = Odor(2,abs(Odor(1,:) + FirstOdorLen + Delay + RespOdorLen - CheckStart(1,itr1)) < 1);
            ResponseOdor = Odor(2,abs(Odor(1,:) + RespOdorLen - CheckStart(1,itr1)) < 1);
            SP = [SP {SpikeTrace(round(CheckEnd(1,itr1)*TimeGain)-ceil((BeforeFirstOdor+ FirstOdorLen + Delay + RespOdorLen + Waterlen)*TimeGain):round(CheckEnd(1,itr1)*TimeGain))}];
            RG = [RG {tempdata(tempdata(:,3) > CheckEnd(1,itr1) - (BeforeFirstOdor+ FirstOdorLen + Delay + RespOdorLen + Waterlen) & tempdata(:,3) < CheckEnd(1,itr1),3)-(CheckEnd(1,itr1) - (BeforeFirstOdor+ FirstOdorLen + Delay + RespOdorLen + Waterlen))}];
            % population variable is registered only one time
            licknuminbin = [];
            if SingleUnitCount == 1
                for t = (CheckEnd(1,itr1) - (BeforeFirstOdor + FirstOdorLen + Delay + RespOdorLen + Waterlen)):1/TimeGain:CheckEnd(1,itr1)-1/TimeGain
                    templicknumofbin = length(find((NewLick>t) & (NewLick<=t+1/TimeGain)));
                    licknuminbin = [licknuminbin templicknumofbin];
                end
                LickRate = [LickRate; licknuminbin];
                LickTime = [LickTime {NewLick(NewLick > CheckEnd(1,itr1) - (BeforeFirstOdor+ FirstOdorLen + Delay + RespOdorLen + Waterlen) & NewLick < CheckEnd(1,itr1),1) - (CheckEnd(1,itr1) - (BeforeFirstOdor+ FirstOdorLen + Delay + RespOdorLen + Waterlen))}];
                if SampleOdor == 1
                    if isempty(find(NewLick > CheckStart(1,itr1) & NewLick < CheckEnd(1,itr1)))
                        TrialMark = [TrialMark; [SampleTestOnsetTs(itr1,:) 1 ResponseOdor 2]];
                    else
                        TrialMark = [TrialMark; [SampleTestOnsetTs(itr1,:) 1 ResponseOdor 1]];
                    end
                else
                    if isempty(find(NewLick > CheckStart(1,itr1) & NewLick < CheckEnd(1,itr1)))
                        TrialMark = [TrialMark; [SampleTestOnsetTs(itr1,:) 2 ResponseOdor 4]];
                    else
                        TrialMark = [TrialMark; [SampleTestOnsetTs(itr1,:) 2 ResponseOdor 3]];
                    end
                end
            end
        end
        Results = CellCombine(Results, SP);
        RGResults = [RGResults; RG];

        %%%%%%%%%%%%%%%%% Laser on-off %%%%%%%%%%%%%%%%%%%%
        if exist('LaserStart') % only for files where laser was applied
            Laser_SP = [];
            Laser_RG = [];
            for itr1 = 1:size(LaserStart,1) % all laser trials
                Laser_SP = [Laser_SP {SpikeTrace(round(LaserStart(itr1,1)*TimeGain-BeforeFirstOdor*TimeGain):(round(LaserStart(itr1,1)*TimeGain)+(Delay+AfterLaser)*TimeGain))}];
                Laser_RG = [Laser_RG {tempdata(tempdata(:,3) > LaserStart(itr1,1)-BeforeFirstOdor & tempdata(:,3) < LaserStart(itr1,1)+Delay+AfterLaser,3)- (LaserStart(itr1,1)-BeforeFirstOdor)}];
            end
            LaserResults = CellCombine(LaserResults, Laser_SP);
            LaserRGResults = [LaserRGResults; Laser_RG];
        end
    else
        DiscardUnitID = [DiscardUnitID itr];
    end
end
if ~isempty(DiscardUnitID)
    for iUnit = 1:length(DiscardUnitID)
        Newdata(Newdata(:,2)==DiscardUnitID(iUnit),:) = [];
    end
    a = unique(Newdata(:,2));
    for i = 1:length(a)
        tempUnitID = a(i);
        Newdata(Newdata(:,2)==tempUnitID,2) = i;
    end
end
if exist('LaserStart')
    save(['short_' FileName(5:end)],'AllTrialsNumber','CompletedTrialID','TrialMark','SingleUnitList','NewWaveForm','LaserResults','LaserRGResults','Results','RGResults','LickRate','LickTime','SpikeTime','Newdata','-v7.3');
else
    save(['short_' FileName(5:end)],'AllTrialsNumber','CompletedTrialID','TrialMark','SingleUnitList','NewWaveForm','Results','RGResults','LickRate','LickTime','SpikeTime','Newdata','-v7.3');
end




