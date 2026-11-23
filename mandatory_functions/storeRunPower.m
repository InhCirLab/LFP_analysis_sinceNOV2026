function ep = storeRunPower(ep, intervals, p)
    ep.thetaRun.newintervals    = intervals;
    ep.thetaRun.thetapowerinuv  = p.rawT;  ep.thetaRun.lgammapowerinuv = p.rawL;  ep.thetaRun.hgammapowerinuv = p.rawH;
    ep.thetaRun.thetapowerindb  = p.dbT;   ep.thetaRun.lgammapowerindb = p.dbL;   ep.thetaRun.hgammapowerindb = p.dbH;
    ep.thetaRun.thetapowerz     = p.zT;    ep.thetaRun.lgammapowerz    = p.zL;    ep.thetaRun.hgammapowerz    = p.zH;
    ep.thetaRun.thetapowerzmod  = p.zmT;   ep.thetaRun.lgammapowerzmod = p.zmL;   ep.thetaRun.hgammapowerzmod = p.zmH;
end