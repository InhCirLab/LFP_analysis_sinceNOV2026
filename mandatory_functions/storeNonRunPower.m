function ep = storeNonRunPower(ep, intervals, p)
    s = 'thetanonrun';
    ep.(s).intervals_stored                     = intervals;
    ep.(s).thetapower_nonTheta_nonRipples_inuv  = p.rawT;  ep.(s).lgammapower_nonTheta_nonRipples_inuv  = p.rawL;  ep.(s).hgammapower_nonTheta_nonRipples_inuv  = p.rawH;
    ep.(s).thetapower_nonTheta_nonRipples_indb  = p.dbT;   ep.(s).lgammapower_nonTheta_nonRipples_indb  = p.dbL;   ep.(s).hgammapower_nonTheta_nonRipples_indb  = p.dbH;
    ep.(s).thetapower_nonTheta_nonRipples_z     = p.zT;    ep.(s).lgammapower_nonTheta_nonRipples_z     = p.zL;    ep.(s).hgammapower_nonTheta_nonRipples_z     = p.zH;
    ep.(s).thetapower_nonTheta_nonRipples_zmod  = p.zmT;   ep.(s).lgammapower_nonTheta_nonRipples_zmod  = p.zmL;   ep.(s).hgammapower_nonTheta_nonRipples_zmod  = p.zmH;
end