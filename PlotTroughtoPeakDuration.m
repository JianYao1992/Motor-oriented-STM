%% Calculate trough-to-peak duration to define broad/narrow spiking neurons.
function PlotTroughtoPeakDuration(Waveform,BrainRegion,BrainRegionUnitsID,threshold,BinWidth,Group)

for iUnit = 1:size(Waveform,1)
    [PeakTroughDuration,~,~,~,~]=CalculatePeakTroughDuration(Waveform(iUnit,:));
    AllUnitsPeakTroughDuration(1,iUnit) = PeakTroughDuration;
end
BrainRegionUnitsPeakTroughDuration = AllUnitsPeakTroughDuration(BrainRegionUnitsID); 
a=min(BrainRegionUnitsPeakTroughDuration); b=max(BrainRegionUnitsPeakTroughDuration);
a=(floor(a/BinWidth)-1)*BinWidth;  b=(ceil(b/BinWidth)+1)*BinWidth; Bin = a:BinWidth:b;
[counts,centers] = hist(BrainRegionUnitsPeakTroughDuration,Bin);
bar(centers,counts/length(BrainRegionUnitsPeakTroughDuration),1,'FaceColor',[0 0 0],'EdgeColor',[.3 .3 .3],'LineWidth',1.5);
hold on
plot([350 350],[0 max(counts/length(BrainRegionUnitsPeakTroughDuration))],'r--','linewidth',4);
set(gca,'XTick',0:400:1200,'XTickLabel',{'0','0.4','0.8','1.2'},'FontName','Arial','FontSize',16,'xlim',[0 1200]);
set(gca,'YTick',0:0.05:0.2,'YTickLabel',{'0','5','10','15','20'},'FontName','Arial','FontSize',16,'ylim',[0 0.2]);
xlabel('Peak-to-valley time (ms)');  ylabel('Proportion (%)'); 
box off;
set(gcf,'Renderer','Painter'); saveas(gcf,[BrainRegion 'PeaktoValleyDurationDistribution'],'fig'); close;
PCID = find(BrainRegionUnitsPeakTroughDuration>threshold); PCDuration = BrainRegionUnitsPeakTroughDuration(PCID);
FSIID = find(BrainRegionUnitsPeakTroughDuration<=threshold); FSIDuration = BrainRegionUnitsPeakTroughDuration(FSIID);
% load control group 
load([BrainRegion 'FRinS1S2_' Group]); % ctrl group: mPFC
FRinAllTrialsforCtrl = cell(1,length(TargetBrainUnitsFRinS1));
for iUnit = 1:length(TargetBrainUnitsFRinS1) % combine S1 and S2 trials
    FRinAllTrialsforCtrl{1,iUnit} = [TargetBrainUnitsFRinS1{iUnit};TargetBrainUnitsFRinS2{iUnit}];
end
for iUnit = 1:length(PCID)
    PCFiring(iUnit) = mean(mean(FRinAllTrialsforCtrl{PCID(iUnit)}(:,1:20)));
end
for iUnit = 1:length(FSIID)
    FSIFiring(iUnit) = mean(mean(FRinAllTrialsforCtrl{FSIID(iUnit)}(:,1:20)));
end
% % plot each neuron point ( x:peak-to-valley duration, y:firing rate )
% figure;
% plot(FSIDuration, FSIFiring, 'marker','o','markerfacecolor',[1 0 0],'markeredgecolor','none','linestyle','none');
% hold on
% plot(PCDuration,PCFiring, 'marker','o','markerfacecolor',[0 0 0],'markeredgecolor','none','linestyle','none');
% hold on
% set(gca,'XTick',0:400:1200,'XTickLabel',{'0','0.4','0.8','1.2'},'FontName','Arial','FontSize',16,'xlim',[0 1200]);
% set(gca,'YTick',0:10:40,'YTickLabel',{'0','10','20','30','40'},'FontName','Arial','FontSize',16,'ylim',[0 35]);
% xlabel('Peak-to-valley time (ms)');  ylabel('Firing rate (Hz)'); 
% box off;
% set(gcf,'Renderer','Painter'); saveas(gcf,[BrainRegion 'PeaktoValleyDurationandbaselineFiringRate'],'fig'); close;
save(['Putative' BrainRegion 'PCandFSIID'],'PCID','FSIID','PCDuration','FSIDuration','PCFiring','FSIFiring','BrainRegionUnitsPeakTroughDuration','-v7.3');

