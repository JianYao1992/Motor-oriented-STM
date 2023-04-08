%% Properties of each unit, including.

clear; clc; close all

%% Assignment
Group = 'CtrlGroup';
beforeneuronnum = 0; 
SamplingRate = 40000;
X = [-19.5:-0.5 0.5:19.5];
WaveForm = [];
MeanWaveform = [];
FASNR = []; 
AutoCorr = []; 
BilateralmPFCMiceID = [{'M19'} {'M20'}]; 
Filename = [];
Filename_CIO = [];
mPFC_ID = [];
aAIC_ID = [];

%% Target directory
CurrPath = uigetdir;
AllPath = genpath(CurrPath);
SplitPath = strsplit(AllPath,';');
SubPath = SplitPath';
SubPath = SubPath(2:end-1);

%% Waveform, FA rate, and SNR of all neurons
for iPath = 1:size(SubPath,1)
    fprintf('IsProcessing_%d_th Path\n',iPath);
    Path = SubPath{iPath,1};
    cd(Path);
    JAVAFiles = dir('*short*.mat');
    JAVAFiles_CIO = dir('*CIO*.mat');
    for j = 1:size(JAVAFiles,1)
        fprintf('IsProcessing_%d_th file\n',j);
        Filename{1,end+1} = JAVAFiles(j,1).name;
        load(Filename{1,end});
        if ~isempty(SingleUnitList)
            % ID of mPFC neurons
            tempmPFCID = [];
            Position = [];
            for k = 1:length(BilateralmPFCMiceID)
                Position = [Position regexp(Filename{1,end},BilateralmPFCMiceID{1,k})]
            end
            if ~isempty(Position)
                tempmPFCID = find(SingleUnitList(:,1)>=1 & SingleUnitList(:,1)<=16);  % mice with recording from bilateral mPFC neurons
            else
                tempmPFCID = find((SingleUnitList(:,1)>=9 & SingleUnitList(:,1)<=16) | (SingleUnitList(:,1)>=25 & SingleUnitList(:,1)<=32));
            end
            if ~isempty(tempmPFCID)
                mPFC_ID = [mPFC_ID tempmPFCID'+beforeneuronnum];
            end
            % ID of aAIC neurons
            tempAIID = [];
            if ~isempty(Position)
                tempAIID = [];
            else
                tempAIID = find((SingleUnitList(:,1)>=1 & SingleUnitList(:,1)<=8) | (SingleUnitList(:,1)>=17 & SingleUnitList(:,1)<=24));
            end
            if ~isempty(tempAIID)
                aAIC_ID = [aAIC_ID tempAIID'+beforeneuronnum];
            end
            beforeneuronnum = beforeneuronnum + size(SingleUnitList,1);
           
           %% Waveform, FA rate, and SNR of individual neuron
            Filename_CIO{1,end+1} = JAVAFiles_CIO(j,1).name;
            load(Filename_CIO{1,end});
            data(:,1) = data(:,1) - 32;
            for iUnit = 1:size(SingleUnitList)
                fprintf('processing %d th neuron\n',iUnit);
                % FA rate in autocorr plot
                WFSampleNum = size(data,2) - 3;
                SpikeTime = data(data(:,1)==SingleUnitList(iUnit,1) & data(:,2)==SingleUnitList(iUnit,2),3);
                tempNewData = data(data(:,1)==SingleUnitList(iUnit,1) & data(:,2)==SingleUnitList(iUnit,2),:);
                AutoCorr{1,end+1} = PlotAutoCorr(tempNewData);
                FalseAlarmRate=sum(AutoCorr{1,end}(21:22))/length(SpikeTime);
                FASNR(1,end+1) = FalseAlarmRate;
                % waveform
                WaveForm{1,end+1} = data(data(:,1)==SingleUnitList(iUnit,1) & data(:,2)==SingleUnitList(iUnit,2),4:end);
                % SNR
                [~,I] = min(NewWaveForm(iUnit,:));
                LargestChannel = ceil(I/(WFSampleNum/4));
                LargestChannelWF = WaveForm{1,end}(:,WFSampleNum/4*(LargestChannel-1)+1:WFSampleNum/4*LargestChannel);
                SNR = SignalNoiseRatioAna(LargestChannelWF);
                FASNR(2,end) = SNR;
            end
            MeanWaveform = [MeanWaveform; NewWaveForm];
        end
    end
end
Region = ones(1,length(WaveForm));
Region(aAIC_ID) = 2;

%% Autocorrelogram plots for qualified neurons
SingCaseNeuronID = find(FASNR(1,:) < 0.0009 & FASNR(2,:) > 10);
cd(SplitPath{1,1});
for iUnit = 1:length(SingCaseNeuronID)
    tempUnitID = SingCaseNeuronID(iUnit);
    if Region(tempUnitID) == 1
        C = [0 0 0];
        Target = 'mPFC';
    else
        C = [0 0 1];
        Target = 'aAIC';
    end
    figure;
    bar(X,AutoCorr{tempUnitID}(1:40),'FaceColor',C,'edgeColor',C); hold on
    rectangle('Position',[-2,0,4,max(AutoCorr{tempUnitID})],'FaceColor',[255 153 199]/255,'EdgeColor','none'); hold on;
    text(-4,max(AutoCorr{tempUnitID}(1:40))*0.97,['SpikeNum=' num2str(size(WaveForm{tempUnitID},1))],'fontsize',12); hold on
    text(-4,max(AutoCorr{tempUnitID}(1:40))*0.9,['FA=' num2str(FASNR(1,tempUnitID))],'fontsize',12); hold on
    text(-4,max(AutoCorr{tempUnitID}(1:40))*0.83,['SNR=' num2str(FASNR(2,tempUnitID))],'fontsize',12); hold on
    box off
    xlabel('Time(ms)','fontsize',12)
    xlim([-20 20])
    ylabel('Counts','fontsize',12)
    ylim([0 max(AutoCorr{tempUnitID}(1:40))*1.01])
    title(['-' num2str(tempUnitID)])
    set(gcf,'Renderer','Painter'); saveas(gcf,[Target 'Unit' num2str(tempUnitID) 'Distribution'],'fig');close all;
end

%% Scatter plot of FA rate and SNR for all neurons
figure;
for iUnit = 1:length(WaveForm)
    if Region(iUnit)==1
        C = [0 0 0];
    else
        C = [0 0 1];
    end
    plot(FASNR(1,iUnit),FASNR(2,iUnit),'marker','o','markerfacecolor',C,'markeredgecolor','none','markersize',3); hold on
end
set(gca,'XTick',0:0.0005:0.0015,'XTickLabel',{'0','0.0005','0.0010','0.0015'},'xlim',[-0.00003 0.0015],'YTick',0:10:30,...
    'YTickLabel',{'0','10','20','30'},'ylim',[0 30]);
box off;
set(gcf,'Renderer','Painter'); saveas(gcf,['UnitsFAandSNR_' Group],'fig');close all;

%% Save results
save([Group 'UnitsWaveformProperties'],'Region','WaveForm','MeanWaveform','FASNR','AutoCorr','-v7.3');

%% Waveform of example neuron
ExamUnitID = 126; % set ID of example neuron
ExamWaveform = WaveForm{ExamUnitID};
ExamAveWaveform = mean(ExamWaveform);
[~,I] = min(ExamAveWaveform);
WFSampleNum = numel(ExamAveWaveform);
LargestChannel = ceil(I/(WFSampleNum/4));
LargestChannelWF = ExamWaveform(:,WFSampleNum/4*(LargestChannel-1)+1:WFSampleNum/4*LargestChannel);
SNR = SignalNoiseRatioAna(LargestChannelWF);
BinSize = 1/SamplingRate*1000;
Bins = BinSize:BinSize:BinSize*WFSampleNum/4;
for i = 2001:2300
    plot(Bins,ExamWaveform(i,WFSampleNum/4*(LargestChannel-1)+1:WFSampleNum/4*LargestChannel),'k')
    hold on
end
text(Bins(1),max(ExamAveWaveform(WFSampleNum/4*(LargestChannel-1)+1:WFSampleNum/4*LargestChannel))*0.97...
    ,['SNR=' num2str(SNR)],'color',[1 0 0])
text(Bins(1),-150,['SpikeNum=' num2str(size(ExamWaveform,1))],'color',[1 0 0])
set(gcf,'renderer','painter'); saveas(gcf,['Unit-' num2str(ExamUnitID) '-Waveform-SNR'],'fig'); close all


