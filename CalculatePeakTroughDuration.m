function [PeakTroughDuration,PeakTroughInfo,FWHM,FWHMLeftPoint,FWHMRightPoint] = CalculatePeakTroughDuration(WaveForm)

WBFs = 40000; %40000 for plexon
WFSampleNum = length(WaveForm)/4;

[~,i] = min(WaveForm);
LargestChannel = ceil(i/WFSampleNum);
LargestChannelWF = WaveForm(WFSampleNum*(LargestChannel-1)+1:WFSampleNum*LargestChannel);
LargestChannelWF1 = interp1(1:length(LargestChannelWF),LargestChannelWF,1:0.5:length(LargestChannelWF)+0.5);%12.5us per sample point

% subplot(1,2,1)
% plot(LargestChannelWF)
% subplot(1,2,2)
% plot(LargestChannelWF1)
% plot(WaveForm)
% hold on
% plot(LargestChannelWF)
%close all
% [MinAmp,MinID]=min(LargestChannelWF);
% [MaxAmp,MaxID]=max(LargestChannelWF(MinID:end));
% MaxID=MaxID+MinID-1;
% PeakTroughInfo=[MaxAmp WFSampleNum*(LargestChannel-1)+MaxID MinAmp WFSampleNum*(LargestChannel-1)+MinID];
[MinAmp,MinID] = min(LargestChannelWF1);
[MaxAmp,MaxID] = max(LargestChannelWF1(MinID:end));
MaxID = MaxID + MinID - 1;
PeakTroughInfo=[MaxAmp WFSampleNum*(LargestChannel-1)+(MaxID+1)/2 MinAmp WFSampleNum*(LargestChannel-1)+(MinID+1)/2];
PeakTroughDuration=(MaxID-MinID)/WBFs/2*10^6;

[~,LeftPoint,RightPoint] = FullWidthHalfMaximumAnalysis(LargestChannelWF1);
FWHM = (RightPoint-LeftPoint)*1/WBFs/2*10^6;
FWHMLeftPoint=WFSampleNum*(LargestChannel-1)+(LeftPoint+1)/2;
FWHMRightPoint=WFSampleNum*(LargestChannel-1)+(RightPoint+1)/2;

% plot(LargestChannelWF1)
% hold on
% plot([LeftPoint LeftPoint],[min(LargestChannelWF1) max(LargestChannelWF1)],'--r')
% plot([RightPoint RightPoint],[min(LargestChannelWF1) max(LargestChannelWF1)],'--r')