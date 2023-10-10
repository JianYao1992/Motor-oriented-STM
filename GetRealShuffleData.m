function [Realdata,Shuffledata] = GetRealShuffleData(UnitBinnedFR,UnittrialID,NumTrialForEachCondition)

%% Real data
[temp1FR,~,~,temp2FR] = ExtractSpecTrialNumForEachNeuron(UnitBinnedFR,UnittrialID,NumTrialForEachCondition,0,0);
Mean_temp1FR = cellfun(@mean,temp1FR,'UniformOutput', 0);
Mean_temp2FR = cellfun(@mean,temp2FR,'UniformOutput', 0);
Realdata = cellfun(@(x,y) abs(x-y)./(x+y),Mean_temp1FR,Mean_temp2FR,'UniformOutput', 0);
Realdata = vertcat(Realdata{:});
a = isnan(Realdata);
for i = 1:size(a,1)
    for j = 1:size(a,2)
        if a(i,j) == 1
            Realdata(i,j) = 0;
        end
    end
end

%% Shuffled data
Shuffledata = [];
for itr = 1:1000
    [temp1FR_shuffle,~,~,temp2FR_shuffle]=ExtractSpecTrialNumForEachNeuron(UnitBinnedFR,UnittrialID,NumTrialForEachCondition,1,0);
    Mean_temp1FR_shuffle = cellfun(@mean,temp1FR_shuffle,'UniformOutput', 0);
    Mean_temp2FR_shuffle = cellfun(@mean,temp2FR_shuffle,'UniformOutput', 0);
    tempShuffledata = cellfun(@(x,y) abs(x-y)./(x+y),Mean_temp1FR_shuffle,Mean_temp2FR_shuffle,'UniformOutput', 0);
    mean_tempshuffledata = vertcat(tempShuffledata{:});
    a = isnan(mean_tempshuffledata);
    for i = 1:size(a,1)
        for j = 1:size(a,2)
            if a(i,j) == 1
                mean_tempshuffledata(i,j) = 0;
            end
        end
    end
    Shuffledata = [Shuffledata; mean(mean_tempshuffledata,1)];
    disp(strcat('///Finish calculating shuffled selectivity_',num2str(itr),'times///'));
end




