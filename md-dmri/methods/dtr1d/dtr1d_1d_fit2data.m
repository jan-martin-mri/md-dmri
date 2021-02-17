function s = dtr1d_1d_fit2data(m, xps)

opt = mdm_opt();
opt = dtr1d_opt(opt);

if m(1)>0
    dtr1d = dtr1d_m2dtr1d(m);
    [dtr1d_nx6,r1,w] = dtr1d_dist2nx6r1w(dtr1d);
    if strcmp(opt.dtr1d.t1_weighting, 'saturation') 
        k = exp(-xps.bt*dtr1d_nx6').*(1 - exp(-xps.tr*r1'));
    elseif strcmp(opt.dtr1d.t1_weighting, 'spin_echo')
        k = exp(-xps.bt*dtr1d_nx6').*(1 - 2*exp((xps.te/2 - xps.tr)*r1') + exp(-xps.tr*r1'));
    end
    s = k*w;
else
    s = zeros(xps.n,1);
end

