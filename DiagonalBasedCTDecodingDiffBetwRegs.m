%% Difference in CTD accuracy between mPFC and aAIC based on the separation from the diagonal

clear; clc; close all;

%% Assignment
TimeGain = 10;

%% Load CTD results of all mPFC and aAIC neurons
load('mPFC Neurons_All_CTDecoding_Recording_CtrlGroup_n=300.mat');
CTDecoding_mPFC = data;
load('aAIC Neurons_All_CTDecoding_Recording_CtrlGroup_n=300.mat');
CTDecoding_aAIC = data;
X = 1:size(data,2); 
TimeLampSampDelayBinID = find(X>10 & X<66);
SampDelayDecodingResult_mPFC = CTDecoding_mPFC(:,TimeLampSampDelayBinID,TimeLampSampDelayBinID);
SampDelayDecodingResult_aAIC = CTDecoding_aAIC(:,TimeLampSampDelayBinID,TimeLampSampDelayBinID);

%% Difference in decoding accuracy between mPFC and aAIC neurons
SampDelayBinNum = numel(TimeLampSampDelayBinID);
ReSampleTimes = size(SampDelayDecodingResult_mPFC,1);
Decoding_mPFC = zeros(ReSampleTimes,SampDelayBinNum);
Decoding_aAIC = zeros(ReSampleTimes,SampDelayBinNum);
DecodingDiff = zeros(ReSampleTimes,SampDelayBinNum);
for iResample = 1:ReSampleTimes
    tempDecoding_mPFC = SampDelayDecodingResult_mPFC(iResample,:,:);
    tempDecoding_mPFC = reshape(tempDecoding_mPFC,SampDelayBinNum,SampDelayBinNum);
    tempDecoding_aAIC = SampDelayDecodingResult_aAIC(iResample,:,:);
    tempDecoding_aAIC = reshape(tempDecoding_aAIC,SampDelayBinNum,SampDelayBinNum);
    tempDiff = tempDecoding_aAIC - tempDecoding_mPFC;
    % difference for statistical significance test //////bootstrap//////
    tempVariousIntervalDecoding_mPFC = cell(1,SampDelayBinNum);
    tempVariousIntervalDecoding_aAIC = cell(1,SampDelayBinNum);
    tempDecodingDiff = cell(1,SampDelayBinNum);
    for iTrainBin = 1:SampDelayBinNum % each train bin
        for iTestBin = 1:SampDelayBinNum % each test bin
            TimeSeparation = abs(iTestBin-iTrainBin);
            tempVariousIntervalDecoding_mPFC{TimeSeparation+1} = [tempVariousIntervalDecoding_mPFC{TimeSeparation+1}; tempDecoding_mPFC(iTrainBin,iTestBin)];
            tempVariousIntervalDecoding_aAIC{TimeSeparation+1} = [tempVariousIntervalDecoding_aAIC{TimeSeparation+1} tempDecoding_aAIC(iTrainBin,iTestBin)];
            tempDecodingDiff{TimeSeparation+1} = [tempDecodingDiff{TimeSeparation+1}; tempDiff(iTrainBin,iTestBin)];
        end
    end
    AverDecoding_mPFC = cellfun(@mean,tempVariousIntervalDecoding_mPFC);
    AverDecoding_aAIC = cellfun(@mean,tempVariousIntervalDecoding_aAIC);
    AverDecodingDiff = cellfun(@mean,tempDecodingDiff);
    Decoding_mPFC(iResample,:) = 100*AverDecoding_mPFC;
    Decoding_aAIC(iResample,:) = 100*AverDecoding_aAIC;
    DecodingDiff(iResample,:) = 100*AverDecodingDiff;
end

%% Plot decoding accuracy aligning to diagonal
figure('OuterPosition',[219 303 420 534]);
seperationtime = 1:size(Decoding_mPFC,2);
Time = [seperationtime/TimeGain, fliplr(seperationtime/TimeGain)];
% mPFC
CI_mPFC = prctile(Decoding_mPFC,[2.5 97.5],1);
average_mPFC = mean(Decoding_mPFC,1);
Highervalue_mPFC = smooth(CI_mPFC(2,:),5)';
Lowervalue_mPFC = smooth(CI_mPFC(1,:),5)';
value_mPFC = [Highervalue_mPFC, fliplr(Lowervalue_mPFC)];
a = fill(Time,value_mPFC,[1 1 0],'edgecolor','none');
alpha(a,0.2); hold on
plot(seperationtime/TimeGain,smooth(average_mPFC,5), 'color',[1 1 0],'linewidth',5); hold on
% aAIC
CI_aAIC = prctile(Decoding_aAIC,[2.5 97.5],1);
average_aAIC = mean(Decoding_aAIC,1);
Highervalue_aAIC = smooth(CI_aAIC(2,:),5)';
Lowervalue_aAIC = smooth(CI_aAIC(1,:),5)';
value_aAIC = [Highervalue_aAIC, fliplr(Lowervalue_aAIC)];
a = fill(Time,value_aAIC,[1 0 0],'edgecolor','none');
alpha(a,0.2); hold on
plot(seperationtime/TimeGain,smooth(average_aAIC,5), 'color',[1 0 0],'linewidth',5); hold on
% bin ID with significant difference
for iBin = 1:size(DecodingDiff,2)
    if nnz(DecodingDiff(:,iBin)<=0)/size(DecodingDiff,1)<=0.05 || nnz(DecodingDiff(:,iBin)>=0)/size(DecodingDiff,1)<=0.05
        patch([iBin-0.5 iBin-0.5 iBin+0.5 iBin+0.5]./TimeGain,[98 100 100 98],[0 0 0],'edgecolor','none');
        hold on
    end
end
set(gca,'XTick',0.1:1:7.1,'XTickLabel',num2cell(0:1:6),'xlim',[0.1 5.5]);
set(gca,'YTick',50:10:100,'YTickLabel',num2cell(50:10:100),'FontName','Arial','FontSize',16,'ylim',[50 100]);
box off
xlabel('Time separation from diagonal (s)');
ylabel('Decoding accuracy (%)');
set(gcf,'Renderer', 'Painter'); saveas(gcf,'Decoding accuracy aligning to diagonal','fig'); close;


