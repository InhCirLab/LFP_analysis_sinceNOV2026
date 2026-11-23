function pow = extractEpochPower(intervals, powts, ...
    rawT, rawL, rawH, dbT, dbL, dbH, zT, zL, zH, zmT, zmL, zmH)
    nEp  = size(intervals, 1);
    flds = {'rawT','rawL','rawH','dbT','dbL','dbH','zT','zL','zH','zmT','zmL','zmH'};
    for f = 1:numel(flds), pow.(flds{f}) = nan(1, nEp); end
    if nEp == 0, return; end
    for jj = 1:nEp
        on  = find(powts >= intervals(jj,1), 1, 'first');
        off = find(powts <= intervals(jj,2), 1, 'last');
        if isempty(on) || isempty(off) || on > off, continue; end
        idx = on:off;
        pow.rawT(jj) = median(rawT(idx)); pow.rawL(jj) = median(rawL(idx)); pow.rawH(jj) = median(rawH(idx));
        pow.dbT(jj)  = median(dbT(idx));  pow.dbL(jj)  = median(dbL(idx));  pow.dbH(jj)  = median(dbH(idx));
        pow.zT(jj)   = median(zT(idx));   pow.zL(jj)   = median(zL(idx));   pow.zH(jj)   = median(zH(idx));
        pow.zmT(jj)  = median(zmT(idx));  pow.zmL(jj)  = median(zmL(idx));  pow.zmH(jj)  = median(zmH(idx));
    end
end