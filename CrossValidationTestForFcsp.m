function Decoding_Results = CrossValidationTestForFcsp(AllPairsFcspRate,AllPairsTrialID,num_trial_ForEachCondition,num_pair_forDecoding...
    ,num_resample_runs,IsShuffleDecoding)

if IsShuffleDecoding == 0
    Decoding_Results = zeros(num_resample_runs,1);
else
    ShuffleDecoding_Results = zeros(num_resample_runs,1);
end

for i = 1:num_resample_runs % neurons and trials change for each resampling
    temp = randperm(length(AllPairsFcspRate));
    targetPairID = temp(1:num_pair_forDecoding); % randomly pick up the neuron ID
    ResampAllPairsFCSPrate = AllPairsFcspRate(:,targetPairID);
    ReasmpAllPairsTrialID = AllPairsTrialID(:,targetPairID);
    
    %% Extract specified trial number for each neuron
    [tempRule1Rate,~,~,tempRule2Rate] = ExtractSpecTrialNumForEachNeuron(ResampAllPairsFCSPrate,ReasmpAllPairsTrialID,num_trial_ForEachCondition,0,0);
    RandPickedTestTriID = randperm(num_trial_ForEachCondition);
    if IsShuffleDecoding == 0
        CVResults = zeros(2*length(RandPickedTestTriID),1);
    else
        [tempRule1Rate_shuffle,~,~,tempRule2Rate_shuffle] = ExtractSpecTrialNumForEachNeuron(ResampAllPairsFCSPrate,ReasmpAllPairsTrialID,num_trial_ForEachCondition,1,0);
        CVResults_shuffle = zeros(2*length(RandPickedTestTriID),1);
    end
    for CVnum = 1:length(RandPickedTestTriID)
        TestTrialID = RandPickedTestTriID(CVnum);
        templateTrialID = setdiff(RandPickedTestTriID,TestTrialID);
        % construct test trial
        TestTrialID = repmat({TestTrialID},size(tempRule1Rate));
        if IsShuffleDecoding == 0
            Rate_test1 = cellfun(@(x,y) x(y,:),tempRule1Rate,TestTrialID,'UniformOutput', 0);
            Rate_test1 = vertcat(Rate_test1{:});
            Rate_test2 = cellfun(@(x,y) x(y,:),tempRule2Rate,TestTrialID,'UniformOutput', 0);
            Rate_test2 = vertcat(Rate_test2{:});
        else
            ShuffleRate_test1 = cellfun(@(x,y) x(y,:),tempRule1Rate_shuffle,TestTrialID,'UniformOutput', 0);
            ShuffleRate_test1 = vertcat(ShuffleRate_test1{:});
            ShuffleRate_test2 = cellfun(@(x,y) x(y,:),tempRule2Rate_shuffle,TestTrialID,'UniformOutput', 0);
            ShuffleRate_test2 = vertcat(ShuffleRate_test2{:});
        end
        % construct template trial
        templateTrialID = repmat({templateTrialID},size(tempRule1Rate));
        if IsShuffleDecoding == 0
            Rate_template1 = cellfun(@(x,y) x(y,:),tempRule1Rate,templateTrialID,'UniformOutput', 0);
            Rate_template2 = cellfun(@(x,y) x(y,:),tempRule2Rate,templateTrialID,'UniformOutput', 0);
        else
            ShuffleRate_template1 = cellfun(@(x,y) x(y,:),tempRule1Rate_shuffle,templateTrialID,'UniformOutput', 0);
            ShuffleRate_template2 = cellfun(@(x,y) x(y,:),tempRule2Rate_shuffle,templateTrialID,'UniformOutput', 0);
        end
        if IsShuffleDecoding == 0
            Rate_template = cellfun(@(x,y) [x;y],Rate_template1,Rate_template2,'UniformOutput', 0);
            Mean = cellfun(@mean,Rate_template,'UniformOutput', 0);
        else
            ShuffleRate_template = cellfun(@(x,y) [x;y],ShuffleRate_template1,ShuffleRate_template2,'UniformOutput', 0);
            Mean_shuffle = cellfun(@mean,ShuffleRate_template,'UniformOutput', 0);
        end
        if IsShuffleDecoding == 1
            Rate_template1 = ShuffleRate_template1;
        end
        for iPair = 1:length(Rate_template1) % go through each FC neuronal pair
            if IsShuffleDecoding == 0
                Rate_template1{iPair} = Rate_template1{iPair}-repmat(Mean{iPair},size(Rate_template1{1},1),1);
                Rate_template2{iPair} = Rate_template2{iPair}-repmat(Mean{iPair},size(Rate_template2{1},1),1);
                Rate_test1(iPair,:) = Rate_test1(iPair,:)-Mean{iPair};
                Rate_test2(iPair,:) = Rate_test2(iPair,:)-Mean{iPair};
            else
                ShuffleRate_template1{iPair} = ShuffleRate_template1{iPair}-repmat(Mean_shuffle{iPair},size(ShuffleRate_template1{1},1),1);
                ShuffleRate_template2{iPair} = ShuffleRate_template2{iPair}-repmat(Mean_shuffle{iPair},size(ShuffleRate_template2{1},1),1);
                ShuffleRate_test1(iPair,:) = ShuffleRate_test1(iPair,:)-Mean_shuffle{iPair};
                ShuffleRate_test2(iPair,:) = ShuffleRate_test2(iPair,:)-Mean_shuffle{iPair};
            end
        end
        if IsShuffleDecoding == 0
            Rate_template1 = cellfun(@mean,Rate_template1,'UniformOutput', 0);
            Rate_template2 = cellfun(@mean,Rate_template2,'UniformOutput', 0);
            Rate_template1 = vertcat(Rate_template1{:});
            Rate_template2 = vertcat(Rate_template2{:});
        else
            ShuffleRate_template1 = cellfun(@mean,ShuffleRate_template1,'UniformOutput', 0);
            ShuffleRate_template2 = cellfun(@mean,ShuffleRate_template2,'UniformOutput', 0);
            ShuffleRate_template1 = vertcat(ShuffleRate_template1{:});
            ShuffleRate_template2 = vertcat(ShuffleRate_template2{:});
        end
        
        %% Correlation coefficient between test trial and template
        if IsShuffleDecoding == 1
            Rate_template1 = ShuffleRate_template1;
        end
        test1template1corr = [];
        test1template2corr = [];
        test2template1corr = [];
        test2template2corr = [];
        test1template1corr_shuffle = [];
        test1template2corr_shuffle = [];
        test2template1corr_shuffle =[];
        test2template2corr_shuffle = [];
        if IsShuffleDecoding == 0
            test1template1corr = min(min(corrcoef(Rate_test1,Rate_template1)));
            test1template2corr = min(min(corrcoef(Rate_test1,Rate_template2)));
            test2template1corr = min(min(corrcoef(Rate_test2,Rate_template1)));
            test2template2corr = min(min(corrcoef(Rate_test2,Rate_template2)));
        else
            test1template1corr_shuffle = min(min(corrcoef(ShuffleRate_test1,ShuffleRate_template1)));
            test1template2corr_shuffle = min(min(corrcoef(ShuffleRate_test1,ShuffleRate_template2)));
            test2template1corr_shuffle = min(min(corrcoef(ShuffleRate_test2,ShuffleRate_template1)));
            test2template2corr_shuffle = min(min(corrcoef(ShuffleRate_test2,ShuffleRate_template2)));
        end
        if IsShuffleDecoding == 0
            CVResults(2*(CVnum-1)+1:2*CVnum,:) = [ceil((test1template1corr-test1template2corr)/10000); ceil((test2template2corr-test2template1corr)/10000)];
        else
            CVResults_shuffle(2*(CVnum-1)+1:2*CVnum,:) = [ceil((test1template1corr_shuffle-test1template2corr_shuffle)/10000); ceil((test2template2corr_shuffle-test2template1corr_shuffle)/10000)];
        end
    end
    if IsShuffleDecoding == 0
        Decoding_Results(i,:) = mean(CVResults);
    else
        ShuffleDecoding_Results(i,:) = mean(CVResults_shuffle);
        Decoding_Results = ShuffleDecoding_Results;
    end
end