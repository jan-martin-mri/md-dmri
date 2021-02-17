function dtr2r1d = dtr2r1d_proliferation(stemp, bt_mx6, te, tr, opt)

dmin = opt.dtr2r1d.dmin;
dmax = opt.dtr2r1d.dmax;
r2min = opt.dtr2r1d.r2min;
r2max = opt.dtr2r1d.r2max;
r1min = opt.dtr2r1d.r1min;
r1max = opt.dtr2r1d.r1max;
n_nodes = opt.dtr2r1d.n_in; % number of nodes
n_proliferation = opt.dtr2r1d.n_proliferation; % number of node sets

dtr2r1d_nodes1 = [];
for niter = 1:n_proliferation    
    % generate random nodes
    dtr2r1d_nodes2 = dtr2r1d_rand(n_nodes,dmin,dmax,r2min,r2max,r1min,r1max);
    % merge node sets
    dtr2r1d_nodes = dtr2r1d_nodes_merge(dtr2r1d_nodes1,dtr2r1d_nodes2);

    % fit to data & get a full distribution
    dtr2r1d = dtr2r1d_data2dtr2r1d(stemp,bt_mx6, te, tr, dtr2r1d_nodes);
    
    dtr2r1d_nodes1 = dtr2r1d_dist2nodes(dtr2r1d); % extract nodes from distribution
end

