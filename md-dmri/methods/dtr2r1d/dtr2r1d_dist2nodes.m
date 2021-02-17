function dtr2r1d_nodes = dtr2r1d_dist2nodes(dtr2r1d)

[n,par,perp,theta,phi,r2,r1,~] = dtr2r1d_dist2par(dtr2r1d);

dtr2r1d_nodes = [par'; perp'; theta'; phi'; r2'; r1'];
dtr2r1d_nodes = [n; dtr2r1d_nodes(:)];
