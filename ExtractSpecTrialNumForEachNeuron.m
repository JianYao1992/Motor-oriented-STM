function [RandPickedFR_Hit,RandPickedFR_Miss,RandPickedFR_FA,RandPickedFR_CR] = ExtractSpecTrialNumForEachNeuron(tempFR,TrialID,TrialNumForEachCondition,IsShuffleDecoding,IsConsiderCorrectError)

% set FRs in hit, miss, FA and CR trials
RandPickedFR_Hit = cell(size(tempFR));
RandPickedFR_Miss = cell(size(tempFR));
RandPickedFR_FA = cell(size(tempFR));
RandPickedFR_CR = cell(size(tempFR));
% FRs in hit, miss, false and CR trials
for iUnit = 1:size(tempFR,2)
    if IsConsiderCorrectError == 0
        if size(TrialID,1) == 2
            TrialIndex_Go = TrialID{1,iUnit}; % Go trials
            TrialIndex_Nogo = TrialID{2,iUnit}; % NoGo trials
        else
            TrialIndex_Go = vertcat(TrialID{1,iUnit},TrialID{2,iUnit}); % Go trials
            TrialIndex_Nogo = vertcat(TrialID{3,iUnit},TrialID{4,iUnit}); % NoGo trials
        end
        if IsShuffleDecoding == 1
            TrialIndex = [TrialIndex_Go; TrialIndex_Nogo];
            temp = randperm(length(TrialIndex));
            TrialIndex_Go = TrialIndex(temp(1:length(TrialIndex_Go)));
            TrialIndex_Nogo = TrialIndex(temp(length(TrialIndex_Go)+1:end));
        end
        temp = randperm(length(TrialIndex_Go)); 
        if ~isempty(TrialNumForEachCondition)
            NewTrialIndex_Go = TrialIndex_Go(temp(1:TrialNumForEachCondition));
        else
            NewTrialIndex_Go = TrialIndex_Go;
        end
        RandPickedFR_Hit{1,iUnit} = tempFR{iUnit}(NewTrialIndex_Go,:);
        temp = randperm(length(TrialIndex_Nogo));
        if ~isempty(TrialNumForEachCondition)
            NewTrialIndex_Nogo = TrialIndex_Nogo(temp(1:TrialNumForEachCondition));
        else
            NewTrialIndex_Nogo = TrialIndex_Nogo;
        end
        RandPickedFR_CR{1,iUnit} = tempFR{iUnit}(NewTrialIndex_Nogo,:);
    else
        TrialIndex_Hit = TrialID{1,iUnit}; % ID of hit trials
        TrialIndex_Miss = TrialID{2,iUnit}; % ID of miss trials
        TrialIndex_FA = TrialID{3,iUnit}; % ID of FA trials
        TrialIndex_CR = TrialID{4,iUnit}; % ID of CR trials
        if IsShuffleDecoding == 1
            TrialIndex = [TrialIndex_FA; TrialIndex_CR];
            temp = randperm(length(TrialIndex));
            TrialIndex_FA = TrialIndex(temp(1:length(TrialIndex_FA)));
            TrialIndex_CR = TrialIndex(temp(length(TrialIndex_FA)+1:end));
        end
        % FRs in hit trials
        temp = randperm(length(TrialIndex_Hit));
        if ~isempty(TrialNumForEachCondition)
            NewTrialIndex_Hit = TrialIndex_Hit(temp(1:TrialNumForEachCondition));
        else
            NewTrialIndex_Hit = TrialIndex_Hit;
        end
        RandPickedFR_Hit{1,iUnit} = tempFR{iUnit}(NewTrialIndex_Hit,:);
        % FRs in miss trials
        temp = randperm(length(TrialIndex_Miss));
        if ~isempty(TrialNumForEachCondition)
            NewTrialIndex_Miss = TrialIndex_Miss(temp(1:TrialNumForEachCondition));
        else
            NewTrialIndex_Miss = TrialIndex_Miss;
        end
        RandPickedFR_Miss{1,iUnit} = tempFR{iUnit}(NewTrialIndex_Miss,:);
        % FRs in FA trials
        temp = randperm(length(TrialIndex_FA));
        if ~isempty(TrialNumForEachCondition)
            NewTrialIndex_FA = TrialIndex_FA(temp(1:TrialNumForEachCondition));
        else
            NewTrialIndex_FA = TrialIndex_FA;
        end
        RandPickedFR_FA{1,iUnit} = tempFR{iUnit}(NewTrialIndex_FA,:);
        % FRs in CR trials
        temp = randperm(length(TrialIndex_CR));
        if ~isempty(TrialNumForEachCondition)
            NewTrialIndex_CR = TrialIndex_CR(temp(1:TrialNumForEachCondition));
        else
            NewTrialIndex_CR = TrialIndex_CR;
        end
        RandPickedFR_CR{1,iUnit} = tempFR{iUnit}(NewTrialIndex_CR,:);
    end
end
