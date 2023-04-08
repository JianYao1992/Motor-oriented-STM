
function PlotHeatMapOfSigSelecForMemoryNeurons(SelectivityDataFileName)

load(SelectivityDataFileName);
ConstructSelectivity = [];

%% Sustained neurons /// sorting by selecvitivy value ///
if ~isempty(SustainedUnitID)
    RowID = []; Selectivity = [];
    for iUnit = 1:length(SustainedUnitID)
        tempRowID = find(SortedID==SustainedUnitID(iUnit));
        RowID = [RowID; tempRowID];
        Selectivity = [Selectivity; TargetBrainSigSelectivity(tempRowID,:)];
    end
    ConstructSelectivity = sortrows([-1*mean(Selectivity(:,21:60),2) SustainedUnitID RowID Selectivity]);
end

%% Neurons with 4-second selectivity of swithced pattern /// sorting by selecvitivy value ///
if ~isempty(SwitchedSustainedUnitID)
    RowID = []; Selectivity = [];
    for iUnit = 1:length(SwitchedSustainedUnitID)
        tempRowID = find(SortedID==SwitchedSustainedUnitID(iUnit));
        RowID = [RowID; tempRowID];
        Selectivity = [Selectivity; TargetBrainSigSelectivity(tempRowID,:)];
    end
    ConstructSelectivity = [ConstructSelectivity; sortrows([-1*mean(Selectivity(:,21:60),2) SwitchedSustainedUnitID RowID Selectivity])];
end

%% Neurons with 3-second selectivity /// sort by selecvitivy onset,then value ///
if ~isempty(UnitIDWith3sCoding)
    RowID_codingfrom2ndsec = []; RowID_codingfrom1stsec = []; 
    Selectivity_codingfrom2ndsec = []; Selectivity_codingfrom1stsec =[];
    UnitID_codingfrom2ndsec = []; UnitID_codingfrom1stsec = [];
    for iUnit = 1:length(UnitIDWith3sCoding)
        tempRowID = find(SortedID==UnitIDWith3sCoding(iUnit));
        if TargetBrainSigSelectivity(tempRowID,21) ~= 0
            UnitID_codingfrom1stsec =  [UnitID_codingfrom1stsec; UnitIDWith3sCoding(iUnit)];
            RowID_codingfrom1stsec = [RowID_codingfrom1stsec; tempRowID];
            Selectivity_codingfrom1stsec = [Selectivity_codingfrom1stsec; TargetBrainSigSelectivity(tempRowID,:)];
        else
            UnitID_codingfrom2ndsec =  [UnitID_codingfrom2ndsec; UnitIDWith3sCoding(iUnit)];
            RowID_codingfrom2ndsec = [RowID_codingfrom2ndsec; tempRowID];
            Selectivity_codingfrom2ndsec = [Selectivity_codingfrom2ndsec; TargetBrainSigSelectivity(tempRowID,:)];
        end
    end
    if ~isempty(UnitID_codingfrom2ndsec)
        ConstructSelectivity = [ConstructSelectivity; sortrows([-1*mean(Selectivity_codingfrom2ndsec(:,21:60),2) UnitID_codingfrom2ndsec RowID_codingfrom2ndsec Selectivity_codingfrom2ndsec])];
    end
    if ~isempty(UnitID_codingfrom1stsec)
        ConstructSelectivity = [ConstructSelectivity; sortrows([-1*mean(Selectivity_codingfrom1stsec(:,21:60),2) UnitID_codingfrom1stsec RowID_codingfrom1stsec Selectivity_codingfrom1stsec])];
    end
end

%% Neurons with 2-second selectivity /// sort by selecvitivy onset,then value ///
if ~isempty(UnitIDWith2sCoding)
    RowID_codingfrom3rdsec = []; RowID_codingfrom2ndsec = []; RowID_codingfrom1stsec = []; 
    Selectivity_codingfrom3rdsec = []; Selectivity_codingfrom2ndsec = []; Selectivity_codingfrom1stsec =[];
    UnitID_codingfrom3rdsec = []; UnitID_codingfrom2ndsec = []; UnitID_codingfrom1stsec = [];
    for iUnit = 1:length(UnitIDWith2sCoding)
        tempRowID = find(SortedID==UnitIDWith2sCoding(iUnit));
        if TargetBrainSigSelectivity(tempRowID,21) ~= 0
            UnitID_codingfrom1stsec =  [UnitID_codingfrom1stsec; UnitIDWith2sCoding(iUnit)];
            RowID_codingfrom1stsec = [RowID_codingfrom1stsec; tempRowID];
            Selectivity_codingfrom1stsec = [Selectivity_codingfrom1stsec; TargetBrainSigSelectivity(tempRowID,:)];
        elseif TargetBrainSigSelectivity(tempRowID,31) ~= 0
            UnitID_codingfrom2ndsec =  [UnitID_codingfrom2ndsec; UnitIDWith2sCoding(iUnit)];
            RowID_codingfrom2ndsec = [RowID_codingfrom2ndsec; tempRowID];
            Selectivity_codingfrom2ndsec = [Selectivity_codingfrom2ndsec; TargetBrainSigSelectivity(tempRowID,:)];
        elseif TargetBrainSigSelectivity(tempRowID,41) ~= 0
            UnitID_codingfrom3rdsec =  [UnitID_codingfrom3rdsec; UnitIDWith2sCoding(iUnit)];
            RowID_codingfrom3rdsec = [RowID_codingfrom3rdsec; tempRowID];
            Selectivity_codingfrom3rdsec = [Selectivity_codingfrom3rdsec; TargetBrainSigSelectivity(tempRowID,:)];
        end
    end
    if ~isempty(UnitID_codingfrom3rdsec)
        ConstructSelectivity = [ConstructSelectivity; sortrows([-1*mean(Selectivity_codingfrom3rdsec(:,21:60),2) UnitID_codingfrom3rdsec RowID_codingfrom3rdsec Selectivity_codingfrom3rdsec])];
    end
    if ~isempty(UnitID_codingfrom2ndsec)
        ConstructSelectivity = [ConstructSelectivity; sortrows([-1*mean(Selectivity_codingfrom2ndsec(:,21:60),2) UnitID_codingfrom2ndsec RowID_codingfrom2ndsec Selectivity_codingfrom2ndsec])];
    end
    if ~isempty(UnitID_codingfrom1stsec)
        ConstructSelectivity = [ConstructSelectivity; sortrows([-1*mean(Selectivity_codingfrom1stsec(:,21:60),2) UnitID_codingfrom1stsec RowID_codingfrom1stsec Selectivity_codingfrom1stsec])];
    end
end

%% Neurons with 1-second selectivity /// sort by selecvitivy onset,then value ///
if ~isempty(UnitIDWith1sCoding)
    RowID_codingfrom4thsec = []; RowID_codingfrom3rdsec = []; RowID_codingfrom2ndsec = []; RowID_codingfrom1stsec = [];
    Selectivity_codingfrom4thsec = []; Selectivity_codingfrom3rdsec = []; Selectivity_codingfrom2ndsec = []; Selectivity_codingfrom1stsec =[];
    UnitID_codingfrom4thsec = []; UnitID_codingfrom3rdsec = []; UnitID_codingfrom2ndsec = []; UnitID_codingfrom1stsec = [];
    for iUnit = 1:length(UnitIDWith1sCoding)
        tempRowID = find(SortedID==UnitIDWith1sCoding(iUnit));
        if TargetBrainSigSelectivity(tempRowID,21) ~= 0
            UnitID_codingfrom1stsec =  [UnitID_codingfrom1stsec; UnitIDWith1sCoding(iUnit)];
            RowID_codingfrom1stsec = [RowID_codingfrom1stsec; tempRowID];
            Selectivity_codingfrom1stsec = [Selectivity_codingfrom1stsec; TargetBrainSigSelectivity(tempRowID,:)];
        elseif TargetBrainSigSelectivity(tempRowID,31) ~= 0
            UnitID_codingfrom2ndsec =  [UnitID_codingfrom2ndsec; UnitIDWith1sCoding(iUnit)];
            RowID_codingfrom2ndsec = [RowID_codingfrom2ndsec; tempRowID];
            Selectivity_codingfrom2ndsec = [Selectivity_codingfrom2ndsec; TargetBrainSigSelectivity(tempRowID,:)];
        elseif TargetBrainSigSelectivity(tempRowID,41) ~= 0
            UnitID_codingfrom3rdsec =  [UnitID_codingfrom3rdsec; UnitIDWith1sCoding(iUnit)];
            RowID_codingfrom3rdsec = [RowID_codingfrom3rdsec; tempRowID];
            Selectivity_codingfrom3rdsec = [Selectivity_codingfrom3rdsec; TargetBrainSigSelectivity(tempRowID,:)];
        elseif TargetBrainSigSelectivity(tempRowID,51) ~= 0
            UnitID_codingfrom4thsec =  [UnitID_codingfrom4thsec; UnitIDWith1sCoding(iUnit)];
            RowID_codingfrom4thsec = [RowID_codingfrom4thsec; tempRowID];
            Selectivity_codingfrom4thsec = [Selectivity_codingfrom4thsec; TargetBrainSigSelectivity(tempRowID,:)];
        end
    end
    if ~isempty(UnitID_codingfrom4thsec)
        ConstructSelectivity = [ConstructSelectivity; sortrows([-1*mean(Selectivity_codingfrom4thsec(:,21:60),2) UnitID_codingfrom4thsec RowID_codingfrom4thsec Selectivity_codingfrom4thsec])];
    end
    if ~isempty(UnitID_codingfrom3rdsec)
        ConstructSelectivity = [ConstructSelectivity; sortrows([-1*mean(Selectivity_codingfrom3rdsec(:,21:60),2) UnitID_codingfrom3rdsec RowID_codingfrom3rdsec Selectivity_codingfrom3rdsec])];
    end
    if ~isempty(UnitID_codingfrom2ndsec)
        ConstructSelectivity = [ConstructSelectivity; sortrows([-1*mean(Selectivity_codingfrom2ndsec(:,21:60),2) UnitID_codingfrom2ndsec RowID_codingfrom2ndsec Selectivity_codingfrom2ndsec])];
    end
    if ~isempty(UnitID_codingfrom1stsec)
        ConstructSelectivity = [ConstructSelectivity; sortrows([-1*mean(Selectivity_codingfrom1stsec(:,21:60),2) UnitID_codingfrom1stsec RowID_codingfrom1stsec Selectivity_codingfrom1stsec])];
    end
end

%% Heatmap of significant selectivity
figure('OuterPosition',[219 303 420 534]);
ColorStep = 0:0.1:1;
IncreaseColor = [ones(length(ColorStep),1) ColorStep' ColorStep'];
ColorStep1 = 1:-0.1:0;
DecreaseColor = [ColorStep1' ColorStep1' ones(length(ColorStep),1)];
Map = [IncreaseColor;DecreaseColor(2:end,:)];
Map = flipud(Map);
colormap(Map);
imagesc(ConstructSelectivity,[-1 1]);
colorbar;
hold on
set(gca,'XTick',[13.5 23.5 33.5 43.5 53.5 63.5 73.5],'XTickLabel',{'0','1','2','3','4', '5','6'},'FontName','Arial','FontSize',16,'xlim',[8.5 73.5]);
set(gca,'YTick',[0 100 200 300 400],'YTickLabel',{'0','100','200','300','400'},'FontName','Arial','FontSize',16);
plot([13.5 13.5],[0 size(ConstructSelectivity,1)+1],'k','linestyle','- -','linewidth',1.5);
plot([23.5 23.5],[0 size(ConstructSelectivity,1)+1],'k','linestyle','- -','linewidth',1.5);
plot([63.5 63.5],[0 size(ConstructSelectivity,1)+1],'k','linestyle','- -','linewidth',1.5);
plot([68.5 68.5],[0 size(ConstructSelectivity,1)+1],'k','linestyle','- -','linewidth',1.5);
set(gcf,'Renderer','Painter'); saveas(gcf,[SelectivityDataFileName(1:4) 'SelectivityHeatMap_SigBinPlot'],'fig');
close;