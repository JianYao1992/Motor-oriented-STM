function SNR = SignalNoiseRatioAna(LargestChannelWF)

AllSpikeAmplitude = zeros(size(LargestChannelWF,1),1);
Baseline = zeros(size(LargestChannelWF,1),1);
for i = 1:size(LargestChannelWF,1)
    tempSpike = LargestChannelWF(i,:);
    [Amplitude,I] = min(tempSpike);
    AllSpikeAmplitude(i) = Amplitude;
    Baseline(i) = mean(tempSpike(1:8));
end
AllSpikeAmplitude = abs(AllSpikeAmplitude);
MeanAmp = mean(AllSpikeAmplitude);
STD = std(Baseline-mean(Baseline));
SNR = MeanAmp/STD;