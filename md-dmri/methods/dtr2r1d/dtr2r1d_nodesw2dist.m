function dtr2r1d = dtr2r1d_nodesw2dist(dtr2r1d_nodes,w)

[n,par,perp,theta,phi,r2,r1] = dtr2r1d_nodes2par(dtr2r1d_nodes);

dtr2r1d = [par'; perp'; theta'; phi'; r2'; r1'; w'];
dtr2r1d = [n; dtr2r1d(:)];
