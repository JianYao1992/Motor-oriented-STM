function FCinEachBin = ConstructFcMatrix(conn_chain,pref_chain,mouse_chain,bin,columnnum,CouplinginEachBin)

for iPair = 1:size(conn_chain,1)
    if ~isempty(CouplinginEachBin)
        idx = find(CouplinginEachBin(:,1) == conn_chain(iPair,1) & CouplinginEachBin(:,2) == conn_chain(iPair,2));
    else
        idx = [];
    end
    if pref_chain(iPair,bin+2) > 0 && pref_chain(iPair,bin+9) > 0
        selectype = 2;
    else
        selectype = 1;
    end
    if isempty(idx)
        CouplinginEachBin(end+1,:) = zeros(1,columnnum);
        CouplinginEachBin(end,1:2) = conn_chain(iPair,1:2);
        CouplinginEachBin(end,bin+2) = selectype;
        CouplinginEachBin(end,columnnum) = mouse_chain(iPair);
    else
        CouplinginEachBin(idx,bin+2) = selectype;
    end
end
FCinEachBin = CouplinginEachBin;