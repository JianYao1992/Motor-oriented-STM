function sumI = SortFcPairs(CouplingValue)

if ~isempty(CouplingValue)
    dura = [];
    for i = 1:size(CouplingValue,1)
        dura = [dura; nnz(CouplingValue(i,3:end-1))];
    end
    sort_mat = CouplingValue(:,3:end-1) > 0;
    act_mat = CouplingValue(:,3:end-1) > 1;
    if size(sort_mat,2) == 4
        sortsum = act_mat*(2.^(0:1:3)') + sort_mat*(2.^(4:1:7)') + dura*(2^8);
    else
        sortsum = act_mat*(2.^(0:1:7)') + sort_mat*(2.^(8:1:15)' + dura*(2^16));
    end
    [~,sumI] = sort(sortsum,'descend');
end  