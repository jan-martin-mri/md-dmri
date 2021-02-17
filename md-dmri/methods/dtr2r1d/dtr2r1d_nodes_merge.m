function dtr2r1d_nodes = dtr2r1d_nodes_merge(dtr2r1d_nodes1,dtr2r1d_nodes2)
% Combine the separate nodes of dtr2r1d_nodes1 and dtr2r1d_nodes2 nodes in 
% single parameter vectors

[n1,par1,perp1,theta1,phi1,r21,r11] = dtr2r1d_nodes2par(dtr2r1d_nodes1);
[n2,par2,perp2,theta2,phi2,r22,r12] = dtr2r1d_nodes2par(dtr2r1d_nodes2);

n = n1 + n2;
par = [par1; par2];
perp = [perp1; perp2];
theta = [theta1; theta2];
phi = [phi1; phi2];
r2 = [r21; r22];
r1 = [r11; r12];

dtr2r1d_nodes = [par'; perp'; theta'; phi'; r2'; r1'];
dtr2r1d_nodes = [n; dtr2r1d_nodes(:)];
