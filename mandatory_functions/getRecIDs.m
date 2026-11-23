function ids = getRecIDs(idtable)
% Devuelve cellstr (n x 1) con el nombre de cada registro a partir de idtable.
    vn  = idtable.Properties.VariableNames;
    pat = {'id','name','rec','sess','folder','animal'};
    col = 1;
    for k = 1:numel(pat)
        hit = find(contains(lower(vn), pat{k}), 1);
        if ~isempty(hit), col = hit; break; end
    end
    x   = idtable.(vn{col});
    ids = cellstr(string(x));
    ids = ids(:);
end