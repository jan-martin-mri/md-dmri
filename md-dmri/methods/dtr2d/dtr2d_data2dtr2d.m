function dtr2d = dtr2d_data2dtr2d(stemp, bt_mx6, te, dtr2d_nodes)

% Combine all the various diffusion tensor parameters in a vector
% written in Mandel notation
[dtr2d_nx6,r2] = dtr2d_nodes2nx6r2(dtr2d_nodes); 

% Modify Kernel to match experimental design
k = exp(-bt_mx6*dtr2d_nx6').*exp(-te*r2'); % Typical 5D D-R2 kernel
snorm = max(stemp);
w = snorm*lsqnonneg(k,stemp/snorm); % linear fit to data

dtr2d = dtr2d_nodesw2dist(dtr2d_nodes,w); % couple nodes with weights
% Sort solutions by descending order of weights &
% keep only the nodes with non-zero weight
dtr2d = dtr2d_sort(dtr2d);




