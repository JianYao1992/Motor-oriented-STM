function FrDiff = CalculateFcPairsFrDiff(BinID,PairsID_S1,PairsID_S2,FR)

FrDiff = [];
if ~isempty(PairsID_S1) && ~isempty(PairsID_S2)
    FcPairsID_Both = intersect(PairsID_S1,PairsID_S2,'rows');
    FcPairsID_S1 = setdiff(PairsID_S1,FcPairsID_Both,'rows');
    FcPairsID_S2 = setdiff(PairsID_S2,FcPairsID_Both,'rows');
elseif ~isempty(PairsID_S1) && isempty(PairsID_S2)
    FcPairsID_Both = [];
    FcPairsID_S1 = PairsID_S1;
    FcPairsID_S2 = [];
elseif isempty(PairsID_S1) && ~isempty(PairsID_S2)
    FcPairsID_Both = [];
    FcPairsID_S1 = [];
    FcPairsID_S2 = PairsID_S2;
end
FcPairsID = [{FcPairsID_S1},{FcPairsID_S2},{FcPairsID_Both}];
for iType = 1:length(FcPairsID)
    if ~isempty(FcPairsID{iType})
        for iPair = 1:size(FcPairsID{iType},1)
            Unit1 = FcPairsID{iType}(iPair,1); Unit2 = FcPairsID{iType}(iPair,2);
            switch iType
                case 1
                    diffFR_Unit1 = FR{Unit1}(1,BinID+2) - FR{Unit1}(2,BinID+2);
                    diffFR_Unit2 = FR{Unit2}(1,BinID+2) - FR{Unit2}(2,BinID+2);
                    FrDiff(end+1,:) = mean(vertcat(diffFR_Unit1,diffFR_Unit2));
                case 2
                    diffFR_Unit1 = FR{Unit1}(2,BinID+2) - FR{Unit1}(1,BinID+2);
                    diffFR_Unit2 = FR{Unit2}(2,BinID+2) - FR{Unit2}(1,BinID+2);
                    FrDiff(end+1,:) = mean(vertcat(diffFR_Unit1,diffFR_Unit2));
                case 3
                    diffFR_Unit1 = FR{Unit1}(1,BinID+2) - FR{Unit1}(2,BinID+2);
                    diffFR_Unit2 = FR{Unit2}(1,BinID+2) - FR{Unit2}(2,BinID+2);
                    FrDiff(end+1,:) = abs(mean(vertcat(diffFR_Unit1,diffFR_Unit2)));
            end
        end
    end
end
