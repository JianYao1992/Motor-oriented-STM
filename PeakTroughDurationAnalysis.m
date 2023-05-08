function [PeakTroughDuration,LargestChannelWF1,PeakTroughInfo,FWHM,FWHMLeftPoint,FWHMRightPoint] = PeakTroughDurationAnalysis(WaveForm)

WBFs = 40000;
WFSampleNum = length(WaveForm)/4;
[~,i] = min(WaveForm);
LargestChannel = ceil(i/WFSampleNum);
LargestChannelWF = WaveForm(WFSampleNum*(LargestChannel-1)+1:WFSampleNum*LargestChannel);
LargestChannelWF1 = interp1(1:length(LargestChannelWF),LargestChannelWF,1:0.5:length(LargestChannelWF)+0.5); % 12.5 usec per sample point
[MinAmp,MinID] = min(LargestChannelWF1);
[MaxAmp,MaxID] = max(LargestChannelWF1(MinID:end));
MaxID = MaxID + MinID - 1;
PeakTroughInfo = [MaxAmp WFSampleNum*(LargestChannel-1)+(MaxID+1)/2 MinAmp WFSampleNum*(LargestChannel-1)+(MinID+1)/2];
PeakTroughDuration = (MaxID-MinID)/WBFs/2*10^6;
[~,LeftPoint,RightPoint] = FullWidthHalfMaximumAnalysis(LargestChannelWF1);
FWHM = (RightPoint-LeftPoint)*1/WBFs/2*10^6;
FWHMLeftPoint = WFSampleNum*(LargestChannel-1)+(LeftPoint+1)/2;
FWHMRightPoint = WFSampleNum*(LargestChannel-1)+(RightPoint+1)/2;
