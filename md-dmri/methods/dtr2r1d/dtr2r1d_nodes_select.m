function dtr2r1d_nodes_out = dtr2r1d_nodes_select(dtr2r1d_nodes,ind)

[~,par,perp,theta,phi,r2,r1] = dtr2r1d_nodes2par(dtr2r1d_nodes);

n = numel(ind);
par = par(ind);
perp = perp(ind);
theta = theta(ind);
phi = phi(ind);
r2 = r2(ind);
r1 = r1(ind);

dtr2r1d_nodes_out = [par'; perp'; theta'; phi'; r2'; r1'];
dtr2r1d_nodes_out = [n; dtr2r1d_nodes_out(:)];
