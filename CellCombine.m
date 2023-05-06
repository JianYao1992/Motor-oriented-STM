%% Combine two cell matrix (just like joint matrix)

function CellMat = CellCombine(CellA, CellB)
CellMat = [];
if length(CellA) > length(CellB)
    for itr1 = 1:length(CellB)
        if length(CellA{itr1}) ~= length(CellB{itr1})
            CellMat = [];
            break
        else
            CellMat{itr1} = [CellA{itr1}; CellB{itr1}];
        end
    end
    for itr1 = length(CellB)+1:length(CellA)
        CellMat{itr1} = CellA{itr1};
    end
elseif length(CellB) > length(CellA)
    for itr1 = 1:length(CellA)
        if length(CellA{itr1}) ~= length(CellB{itr1})
            CellMat = [];
            break
        else
            CellMat{itr1}=[CellA{itr1}; CellB{itr1}];
        end
    end
    for itr1 = length(CellA)+1:length(CellB)
        CellMat{itr1} = CellB{itr1};
    end   
else
    for itr1 = 1:length(CellA)
        if length(CellA{itr1}) ~= length(CellB{itr1})
            CellMat=[];
            break
        else
            CellMat{itr1} = [CellA{itr1}; CellB{itr1}];
        end
    end
end
end