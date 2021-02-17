function dtr2d_plot(S, xps, axh, axh2)
% function dtr2d_plot(S, xps, axh, axh2)

if (nargin < 4), axh2 = []; end

opt = mdm_opt();
opt = dtr2d_opt(opt);
opt = mplot_opt(opt);

% Customize options
opt.dtd.dmin = .05/max(xps.b);

S = abs(S);

% Show signal and fit
m = mplot_signal_and_fit(S, xps, @dtr2d_1d_data2fit, @dtr2d_1d_fit2data, axh, opt);

% Plot the tensor distribution
[~,dpar,dperp,theta,phi,r2,w] = dtr2d_dist2par(dtr2d_m2dtr2d(m));

diso = (dpar+2*dperp)/3;
dratio = dpar./dperp;

% Set marker sizes
w = w / sum(w);
ms = 30 * opt.mplot.ms * sqrt(w);
ms = real(ms);

dist_d.x = log10(diso);
dist_d.y = log10(dratio);
dist_d.z = log10(r2);
dist_d.a = ms;
dist_d.r = abs(cos(phi).*sin(theta));
dist_d.g = abs(sin(phi).*sin(theta));
dist_d.b = abs(cos(theta));
dist_d.bright = (abs(dpar-dperp)./max([dpar dperp],[],2)).^2;

axpars.xmin = -10; axpars.xmax = -8; 
axpars.ymin = -2; axpars.ymax = 2; 
axpars.zmin = -.5; axpars.zmax = 2;

contourpars.Nx = 50; 
contourpars.Ny = contourpars.Nx; 
contourpars.Nz = contourpars.Nx;
contourpars.Nlevels = 5;

h = mplot_dtr2d(dist_d, axh2, opt, axpars, contourpars);

set(h,'ZTick',0:.5:2)

    
% dps1d = dtr2d_dist2dps1d(dtr2d);