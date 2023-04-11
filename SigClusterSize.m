function [ResultsSigClusSize,ResultsSigClusID] = SigClusterSize(TimelampIsSigBin)

TestBinNum = length(TimelampIsSigBin);
ResultsSigClusSize = [];
ResultsSigClusID = [];
SigClusSize = 0;
SigClusBinID = [];
for iTestBin = 2:TestBinNum
    if iTestBin == 2 && TimelampIsSigBin(iTestBin-1) ~= 0
        SigClusSize = 1;
        SigClusBinID = [SigClusBinID 1];
    end
    if TimelampIsSigBin(iTestBin) ~= 0 && TimelampIsSigBin(iTestBin-1) ~= 0
        SigClusSize = SigClusSize+1;
        SigClusBinID = [SigClusBinID iTestBin];
        if iTestBin == TestBinNum
            ResultsSigClusSize = [ResultsSigClusSize SigClusSize];
            ResultsSigClusID = [ResultsSigClusID {SigClusBinID}];
        end
    elseif TimelampIsSigBin(iTestBin) ~= 0 && TimelampIsSigBin(iTestBin-1) == 0
        SigClusSize = 1;
        SigClusBinID = iTestBin;
    else
        if SigClusSize > 0
            ResultsSigClusSize = [ResultsSigClusSize SigClusSize];
            ResultsSigClusID = [ResultsSigClusID {SigClusBinID}];
        end
        SigClusSize = 0;
        SigClusBinID = [];
    end
end
