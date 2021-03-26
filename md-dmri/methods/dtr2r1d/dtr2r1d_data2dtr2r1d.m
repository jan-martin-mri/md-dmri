function dtr2r1d = dtr2r1d_data2dtr2r1d(stemp, bt_mx6, te, tr, dtr2r1d_nodes)

[dt_nx6,r2,r1] = dtr2r1d_nodes2nx6r2r1(dtr2r1d_nodes);

% Modify Kernel to match experimental design
% jm, 2021-03-26: Use a single line to speed up computation for now.
k = exp(-bt_mx6 * dt_nx6') .* exp(-te * r2') .* (1 + exp(-tr * r1') - 2 * exp((te/2 - tr) * r1'));
%if strcmp(opt.dtr2r1d.t1_weighting, 'saturation')
%    k = exp(-bt_mx6 * dt_nx6') .* exp(-te * r2') .* (1 - exp(-tr * r1'));
%elseif strcmp(opt.dtr2r1d.t1_weighting, 'spin_echo')
%    k = exp(-bt_mx6 * dt_nx6') .* exp(-te * r2') .* (1 + exp(-tr * r1') - 2 * exp((te/2 - tr) * r1'));
%end
snorm = max(stemp);
w = snorm*lsqnonneg(k,stemp/snorm); % linear fit to data

dtr2r1d = dtr2r1d_nodesw2dist(dtr2r1d_nodes,w); % couple nodes with weights
% Sort solutions by descending order of weights &
% keep only the nodes with non-zero weight
dtr2r1d = dtr2r1d_sort(dtr2r1d); 




