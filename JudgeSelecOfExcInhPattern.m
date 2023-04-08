function [SelecDirection,NeuType,GroupID] = JudgeSelecOfExcInhPattern(SecondSelectivity,IsSigModulation,DelayFR)

%% Maximal difference in delay-period FR between S1 and S2 trials
DelayFRdiff = DelayFR(1,:)-DelayFR(2,:);
SigSelecBin = find(SecondSelectivity~=0);
[~,MaxDiffFRBinID] = max(abs(DelayFRdiff(SigSelecBin)));
I = SigSelecBin(MaxDiffFRBinID); % ID of bin with the largest difference in FR
DiffFRpeak = DelayFRdiff(I); % largetst value of difference in FR

%% Group type of neuron
GroupID = 7;
if DiffFRpeak > 0 % Go-preferred
    SelecDirection = 1;
    % default setting
    NeuType = '-GoSelec-UnModulated';
    % number of bins in Go trials showing FR modulation relative to baseline period
    ExiBinsNum_Go = length(find(IsSigModulation(2,2:5)==1));
    InhBinsNum_Go = length(find(IsSigModulation(2,2:5)==-1));
    % ID of bin with the largest and smallest value of FR in Go trials
    [~,MaxDelayBin_Go] = max(DelayFR(1,:));
    [~,MinDelayBin_Go] = min(DelayFR(1,:));
    % all possible conditions
    if length(SigSelecBin) == 1 && IsSigModulation(2,I+1) == 1 || (IsSigModulation(2,I+1) ~= -1 && IsSigModulation(2,MaxDelayBin_Go+1) == 1 && InhBinsNum_Go < 2) 
        GroupID = 1;
        NeuType = '-GoSelec-GoExi';
    elseif length(SigSelecBin) == 1 && (IsSigModulation(2,I+1) == -1 || InhBinsNum_Go >= 3 || (InhBinsNum_Go == 2 && ExiBinsNum_Go < 2 && IsSigModulation(2,I+1) ~= 1)) || (IsSigModulation(2,I+1) ~= 1 && IsSigModulation(2,MinDelayBin_Go+1) == -1 && ExiBinsNum_Go < 2) 
        GroupID = 2;
        NeuType = '-GoSelec-GoInh';
    elseif length(SigSelecBin) == 1 && IsSigModulation(3,I+1) == -1
        GroupID = 3;
        NeuType = '-GoSelec-NogoInh';
    elseif IsSigModulation(2,end) == 1 || IsSigModulation(2,I+1) == 1 || (IsSigModulation(2,MaxDelayBin_Go+1) == 1 && InhBinsNum_Go < 2) 
        if IsSigModulation(2,I+1) == -1
            GroupID = 2;
            NeuType = '-GoSelec-GoInh';
        else
            GroupID = 1;
            NeuType = '-GoSelec-GoExi';
        end
    elseif IsSigModulation(2,end) == -1 || IsSigModulation(2,I+1) == -1 || (IsSigModulation(2,MinDelayBin_Go+1) == -1 && ExiBinsNum_Go < 2)
        if IsSigModulation(2,I+1) == 1
            GroupID = 1;
            NeuType = '-GoSelec-GoExi';
        else
            GroupID = 2;
            NeuType = '-GoSelec-GoInh';
        end
    elseif IsSigModulation(3,end) == -1 || IsSigModulation(3,I+1) == -1 
        GroupID = 3;
        NeuType = '-GoSelec-NogoInh';        
    end
elseif DiffFRpeak < 0 % NoGo-preferred
    SelecDirection = 2;
    % default setting
    NeuType = '-NogoSelec-UnModulated';
    % number of bins in NoGo trials showing FR modulation relative to baseline period
    ExiBinNum_Nogo = length(find(IsSigModulation(3,2:5)==1));
    InhBinNum_Nogo = length(find(IsSigModulation(3,2:5)==-1));
    % ID of bin with the largest and smallest value of FR in NoGo trials
    [~,MaxDelayBin_Nogo] = max(DelayFR(2,:));
    [~,MinDelayBin_Nogo] = min(DelayFR(2,:));
    % all possible conditions
    if length(SigSelecBin) == 1 && IsSigModulation(3,I+1) == 1 || (IsSigModulation(3,MaxDelayBin_Nogo+1) == 1 && InhBinNum_Nogo < 2 && IsSigModulation(3,I+1) ~= -1)
        GroupID = 4;
        NeuType = '-NogoSelec-NogoExi';
    elseif length(SigSelecBin) == 1 && (IsSigModulation(3,I+1) == -1 || InhBinNum_Nogo >= 3 || (InhBinNum_Nogo == 2 && ExiBinNum_Nogo < 2 && IsSigModulation(3,I+1) ~= 1))||(IsSigModulation(3,MinDelayBin_Nogo+1) == -1 && ExiBinNum_Nogo < 2 && IsSigModulation(3,I+1) ~= 1)
        GroupID = 5;
        NeuType = '-NogoSelec-NogoInh';
    elseif length(SigSelecBin) == 1 && IsSigModulation(2,I+1) == -1
        GroupID = 6;
        NeuType = '-NogoSelec-GoInh';
    elseif IsSigModulation(3,end) == 1 || IsSigModulation(3,I+1) ==1 || (IsSigModulation(3,MaxDelayBin_Nogo+1) == 1 && InhBinNum_Nogo < 2)
        if IsSigModulation(3,I+1) == -1
            GroupID = 5;
            NeuType = '-NogoSelec-NogoInh';
        else
            GroupID = 4;
            NeuType = '-NogoSelec-NogoExi';
        end
    elseif IsSigModulation(3,end) == -1 || IsSigModulation(3,I+1) == -1 || (IsSigModulation(3,MinDelayBin_Nogo+1) == -1 && ExiBinNum_Nogo < 2)
        if IsSigModulation(3,I+1) == 1
            GroupID = 4;
            NeuType = '-NogoSelec-NogoExi';
        else
            GroupID = 5;
            NeuType = '-NogoSelec-NogoInh';
        end
    elseif IsSigModulation(2,end) == -1 || IsSigModulation(2,I+1) == -1
        GroupID = 6;
        NeuType = '-NogoSelec-GoInh';
    end
end