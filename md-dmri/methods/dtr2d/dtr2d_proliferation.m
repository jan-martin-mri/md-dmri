function dtr2d = dtr2d_proliferation(stemp, bt_mx6, te, opt)
% Generates M sets of N nodes &
% sequentially fits them to the signal data

dmin = opt.dtr2d.dmin; % min diffusivity
dmax = opt.dtr2d.dmax; % max diffusivity
r2min = opt.dtr2d.r2min; % max r2-rate
r2max = opt.dtr2d.r2max; % max r2-rate
n_nodes = opt.dtr2d.n_in; % number of nodes (N)
n_proliferation = opt.dtr2d.n_proliferation; % number of node sets (M)

dtr2d_nodes1 = [];
for niter = 1:n_proliferation    
    % generate random nodes
    dtr2d_nodes2 = dtr2d_rand(n_nodes,dmin,dmax,r2min,r2max);
    % merge node sets
    dtr2d_nodes = dtr2d_nodes_merge(dtr2d_nodes1,dtr2d_nodes2);

    % fit to data & get a full distribution
    dtr2d = dtr2d_data2dtr2d(stemp,bt_mx6, te, dtr2d_nodes);
    dtr2d_nodes1 = dtr2d_dist2nodes(dtr2d); % extract nodes from distribution
end

