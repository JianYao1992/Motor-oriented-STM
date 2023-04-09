function DECODING_RESULTS = CrossValidationTest(TrialBinnedFR,AllNeuSampleTrialID,num_trial_ForEachCondition,num_neuron_forDecoding,num_resample_runs,IsShuffleDecoding,IsCTDecoding)

BinNum = size(TrialBinnedFR{1},2);
if IsCTDecoding == 0
    if IsShuffleDecoding == 0
        DECODING_RESULTS = zeros(num_resample_runs,BinNum);
    else
        ShuffleDECODING_RESULTS = zeros(num_resample_runs,BinNum);
    end
else
    if IsShuffleDecoding == 0
        DECODING_RESULTS = zeros(num_resample_runs,BinNum,BinNum);
    else
        ShuffleDECODING_RESULTS = zeros(num_resample_runs,BinNum,BinNum);
    end
end

%% Decoding analysis for randomly selected neurons and trials
for i = 1:num_resample_runs % neurons and trials change at each resampling run
    temp = randperm(length(TrialBinnedFR));
    targetNeuronID = temp(1:num_neuron_forDecoding);
    temp_AllUnitBinnedFR = TrialBinnedFR(targetNeuronID);
    temp_AllUnitTrialID = AllNeuSampleTrialID(:,targetNeuronID);
    
    %% Extract specified trials for each neuron
    [tempRule1FR,~,~,tempRule2FR] = ExtractSpecTrialNumForEachNeuron(temp_AllUnitBinnedFR,temp_AllUnitTrialID,num_trial_ForEachCondition,0,0);
    [tempRule1FR_shuffle,~,~,tempRule2FR_shuffle] = ExtractSpecTrialNumForEachNeuron(temp_AllUnitBinnedFR,temp_AllUnitTrialID,num_trial_ForEachCondition,1,0);
    randompickedtesttrialID = randperm(num_trial_ForEachCondition);
    
    %% Leave one out
    if IsCTDecoding == 1
        if IsShuffleDecoding == 0
            CVResults = zeros(2*length(randompickedtesttrialID),size(tempRule1FR{1,1},2),size(tempRule1FR{1,1},2));
        else
            CVResults_shuffle = zeros(2*length(randompickedtesttrialID),size(tempRule1FR{1,1},2),size(tempRule1FR{1,1},2));
        end
    else
        if IsShuffleDecoding == 0
            CVResults = zeros(2*length(randompickedtesttrialID),size(tempRule1FR{1,1},2));
        else
            CVResults_shuffle = zeros(2*length(randompickedtesttrialID),size(tempRule1FR{1,1},2));
        end
    end
    for CVnum = 1:length(randompickedtesttrialID)
        testtrialID = randompickedtesttrialID(CVnum);
        templatetrialID = setdiff(randompickedtesttrialID,testtrialID);
        % construct test trial
        testtrialID = repmat({testtrialID},size(tempRule1FR));
        if IsShuffleDecoding == 0
            test1BinnedFR = cellfun(@(x,y) x(y,:),tempRule1FR,testtrialID,'UniformOutput', 0);
            test1BinnedFR = vertcat(test1BinnedFR{:});
            test2BinnedFR = cellfun(@(x,y) x(y,:),tempRule2FR,testtrialID,'UniformOutput', 0);
            test2BinnedFR = vertcat(test2BinnedFR{:});
        else
            test1BinnedFR_shuffle = cellfun(@(x,y) x(y,:),tempRule1FR_shuffle,testtrialID,'UniformOutput', 0);
            test1BinnedFR_shuffle = vertcat(test1BinnedFR_shuffle{:});
            test2BinnedFR_shuffle = cellfun(@(x,y) x(y,:),tempRule2FR_shuffle,testtrialID,'UniformOutput', 0);
            test2BinnedFR_shuffle = vertcat(test2BinnedFR_shuffle{:});
        end
        % construct template trial
        templatetrialID = repmat({templatetrialID},size(tempRule1FR));
        if IsShuffleDecoding == 0
            template1BinnedFR = cellfun(@(x,y) x(y,:),tempRule1FR,templatetrialID,'UniformOutput', 0);
            template2BinnedFR = cellfun(@(x,y) x(y,:),tempRule2FR,templatetrialID,'UniformOutput', 0);
        else
            template1BinnedFR_shuffle = cellfun(@(x,y) x(y,:),tempRule1FR_shuffle,templatetrialID,'UniformOutput', 0);
            template2BinnedFR_shuffle = cellfun(@(x,y) x(y,:),tempRule2FR_shuffle,templatetrialID,'UniformOutput', 0);
        end
        % matrix of FR of template and test trials in decoding
        if IsShuffleDecoding == 0
            templateBinnedFR = cellfun(@(x,y) [x;y],template1BinnedFR,template2BinnedFR,'UniformOutput', 0);
            Mean = cellfun(@mean,templateBinnedFR,'UniformOutput', 0);
            Std = cellfun(@std,templateBinnedFR,'UniformOutput', 0);
            Std = cellfun(@(x) x+~x,Std,'UniformOutput', 0);
        else
            templateBinnedFR_shuffle = cellfun(@(x,y) [x;y],template1BinnedFR_shuffle,template2BinnedFR_shuffle,'UniformOutput', 0);
            Mean_shuffle = cellfun(@mean,templateBinnedFR_shuffle,'UniformOutput', 0);
            Std_shuffle = cellfun(@std,templateBinnedFR_shuffle,'UniformOutput', 0);
            Std_shuffle = cellfun(@(x) x+~x,Std_shuffle,'UniformOutput', 0);
        end
        if IsShuffleDecoding == 1
            template1BinnedFR = template1BinnedFR_shuffle;
        end
        for itru = 1:length(template1BinnedFR) % neuron
            if IsShuffleDecoding == 0
                template1BinnedFR{itru} = (template1BinnedFR{itru}-repmat(Mean{itru},size(template1BinnedFR{1},1),1))./repmat(Std{itru},size(template1BinnedFR{1},1),1);
                template2BinnedFR{itru} = (template2BinnedFR{itru}-repmat(Mean{itru},size(template2BinnedFR{1},1),1))./repmat(Std{itru},size(template2BinnedFR{1},1),1);
                test1BinnedFR(itru,:) = (test1BinnedFR(itru,:)-Mean{itru})./ Std{itru};
                test2BinnedFR(itru,:) = (test2BinnedFR(itru,:)-Mean{itru})./Std{itru};
            else
                template1BinnedFR_shuffle{itru} = (template1BinnedFR_shuffle{itru}-repmat(Mean_shuffle{itru},size(template1BinnedFR_shuffle{1},1),1))./repmat(Std_shuffle{itru},size(template1BinnedFR_shuffle{1},1),1);
                template2BinnedFR_shuffle{itru} = (template2BinnedFR_shuffle{itru}-repmat(Mean_shuffle{itru},size(template2BinnedFR_shuffle{1},1),1))./repmat(Std_shuffle{itru},size(template2BinnedFR_shuffle{1},1),1);
                test1BinnedFR_shuffle(itru,:) = (test1BinnedFR_shuffle(itru,:)-Mean_shuffle{itru})./Std_shuffle{itru};
                test2BinnedFR_shuffle(itru,:) = (test2BinnedFR_shuffle(itru,:)-Mean_shuffle{itru})./Std_shuffle{itru};
            end
        end
        if IsShuffleDecoding == 0
            template1BinnedFR = cellfun(@mean,template1BinnedFR,'UniformOutput', 0);
            template2BinnedFR = cellfun(@mean,template2BinnedFR,'UniformOutput', 0);
            template1BinnedFR = vertcat(template1BinnedFR{:});
            template2BinnedFR = vertcat(template2BinnedFR{:});
        else
            template1BinnedFR_shuffle = cellfun(@mean,template1BinnedFR_shuffle,'UniformOutput', 0);
            template2BinnedFR_shuffle = cellfun(@mean,template2BinnedFR_shuffle,'UniformOutput', 0);
            template1BinnedFR_shuffle = vertcat(template1BinnedFR_shuffle{:});
            template2BinnedFR_shuffle = vertcat(template2BinnedFR_shuffle{:});
        end
        % correlation coefficient between test trial and template
        if IsShuffleDecoding == 1
            template1BinnedFR = template1BinnedFR_shuffle;
        end
        if IsCTDecoding == 1 
            for itrtrain = 1:size(template1BinnedFR,2)
                for itrtest = 1:size(template1BinnedFR,2)
                    if IsShuffleDecoding == 0
                        test1template1corr(itrtrain,itrtest) = min(min(corrcoef(test1BinnedFR(:,itrtest),template1BinnedFR(:,itrtrain))));
                        test1template2corr(itrtrain,itrtest) = min(min(corrcoef(test1BinnedFR(:,itrtest),template2BinnedFR(:,itrtrain))));
                        test2template1corr(itrtrain,itrtest) = min(min(corrcoef(test2BinnedFR(:,itrtest),template1BinnedFR(:,itrtrain))));
                        test2template2corr(itrtrain,itrtest) = min(min(corrcoef(test2BinnedFR(:,itrtest),template2BinnedFR(:,itrtrain))));
                    else
                        test1template1corr_shuffle(itrtrain,itrtest) = min(min(corrcoef(test1BinnedFR_shuffle(:,itrtest),template1BinnedFR_shuffle(:,itrtrain))));
                        test1template2corr_shuffle(itrtrain,itrtest) = min(min(corrcoef(test1BinnedFR_shuffle(:,itrtest),template2BinnedFR_shuffle(:,itrtrain))));
                        test2template1corr_shuffle(itrtrain,itrtest) = min(min(corrcoef(test2BinnedFR_shuffle(:,itrtest),template1BinnedFR_shuffle(:,itrtrain))));
                        test2template2corr_shuffle(itrtrain,itrtest) = min(min(corrcoef(test2BinnedFR_shuffle(:,itrtest),template2BinnedFR_shuffle(:,itrtrain))));
                    end
                end
            end
            if IsShuffleDecoding == 0
                CVResults(2*(CVnum-1)+1,:,:) = ceil((test1template1corr-test1template2corr)/10000);
                CVResults(2*CVnum,:,:) = ceil((test2template2corr-test2template1corr)/10000);
            else
                CVResults_shuffle(2*(CVnum-1)+1,:,:) = ceil((test1template1corr_shuffle-test1template2corr_shuffle)/10000);
                CVResults_shuffle(2*CVnum,:,:) = ceil((test2template2corr_shuffle-test2template1corr_shuffle)/10000);
            end
        else
            for itrbin = 1:size(template1BinnedFR,2)
                if IsShuffleDecoding == 0
                    test1template1corr(itrbin) = min(min(corrcoef(test1BinnedFR(:,itrbin),template1BinnedFR(:,itrbin))));
                    test1template2corr(itrbin) = min(min(corrcoef(test1BinnedFR(:,itrbin),template2BinnedFR(:,itrbin))));
                    test2template1corr(itrbin) = min(min(corrcoef(test2BinnedFR(:,itrbin),template1BinnedFR(:,itrbin))));
                    test2template2corr(itrbin) = min(min(corrcoef(test2BinnedFR(:,itrbin),template2BinnedFR(:,itrbin))));
                else
                    test1template1corr_shuffle(itrbin) = min(min(corrcoef(test1BinnedFR_shuffle(:,itrbin),template1BinnedFR_shuffle(:,itrbin))));
                    test1template2corr_shuffle(itrbin) = min(min(corrcoef(test1BinnedFR_shuffle(:,itrbin),template2BinnedFR_shuffle(:,itrbin))));
                    test2template1corr_shuffle(itrbin) = min(min(corrcoef(test2BinnedFR_shuffle(:,itrbin),template1BinnedFR_shuffle(:,itrbin))));
                    test2template2corr_shuffle(itrbin) = min(min(corrcoef(test2BinnedFR_shuffle(:,itrbin),template2BinnedFR_shuffle(:,itrbin))));
                end
            end
            if IsShuffleDecoding == 0
                CVResults(2*(CVnum-1)+1:2*CVnum,:) = [ceil((test1template1corr-test1template2corr)/10000);ceil((test2template2corr-test2template1corr)/10000)];
            else
                CVResults_shuffle(2*(CVnum-1)+1:2*CVnum,:) = [ceil((test1template1corr_shuffle-test1template2corr_shuffle)/10000);ceil((test2template2corr_shuffle-test2template1corr_shuffle)/10000)];
            end
        end
    end
    if IsCTDecoding == 1
        if IsShuffleDecoding == 0
            CVResults = mean(CVResults);
            CVResults = reshape(CVResults,size(CVResults,2),size(CVResults,3));
            DECODING_RESULTS(i,:,:) = CVResults;
        else
            CVResults_shuffle = mean(CVResults_shuffle);
            CVResults_shuffle = reshape(CVResults_shuffle,size(CVResults_shuffle,2),size(CVResults_shuffle,3));
            ShuffleDECODING_RESULTS(i,:,:) = CVResults_shuffle;
            DECODING_RESULTS = ShuffleDECODING_RESULTS;
        end
    else
        if IsShuffleDecoding == 0
            DECODING_RESULTS(i,:) = mean(CVResults);
        else
            ShuffleDECODING_RESULTS(i,:) = mean(CVResults_shuffle);
            DECODING_RESULTS = ShuffleDECODING_RESULTS;
        end
    end
end