function s = dtr2r1d_1d_fit2data(m, xps)

opt = mdm_opt();
opt = dtr2r1d_opt(opt);

if m(1)>0
    dtr2r1d = dtr2r1d_m2dtr2r1d(m);
    [dtr2d_nx6,r2,r1,w] = dtr2r1d_dist2nx6r2r1w(dtr2r1d);
    if strcmp(opt.dtr2r1d.t1_weighting, 'saturation')
        k = exp(-xps.bt*dtr2d_nx6').*exp(-xps.te*r2').*(1-exp(-xps.tr*r1'));
    elseif strcmp(opt.dtr2r1d.t1_weighting, 'spin_echo')
        k = exp(-xps.bt*dtr2d_nx6').*exp(-xps.te*r2').*(1 - 2*exp((xps.te/2 - xps.tr)*r1') + exp(-xps.tr*r1'));
    elseif strcmp(opt.dtr2r1d.t1_weighting, 'wip_inversion')
        k = exp(-xps.bt*dtr2d_nx6') .* exp(-xps.te*r2') .* (1 - 2*exp(-xps.ti*r1') + exp(-xps.tr*r1'));
    end
    s = k*w;
else
    s = zeros(xps.n,1);
end

