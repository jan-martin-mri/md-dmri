function dtr2d = dtr2d_nodesw2dist(dtr2d_nodes,w)
% Combines dtr2d_nodes and w in a single distribution vector (dtr2d)

[n,par,perp,theta,phi,r2] = dtr2d_nodes2par(dtr2d_nodes);

dtr2d = [par'; perp'; theta'; phi'; r2'; w'];
dtr2d = [n; dtr2d(:)];
