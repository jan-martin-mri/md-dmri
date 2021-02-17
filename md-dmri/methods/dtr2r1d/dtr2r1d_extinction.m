function dtr2r1d = dtr2r1d_extinction(stemp, bt_mx6, te, tr, dtr2r1d, opt)

n_nodes = opt.dtr2r1d.n_in; % number of nodes
n_extinction = opt.dtr2r1d.n_extinction;

for niter = 1:n_extinction    
    n_in = dtr2r1d(1);
%     dtr2r1d*1e-8
    n_max = min([(n_in-opt.dtr2r1d.n_kill) opt.dtr2r1d.n_out]);

    dtr2r1d_nodes1 = dtr2r1d_dist2nodes(dtr2r1d);
    ind = 1:n_max;
    % pick the n_max sols with higher weights
    dtr2r1d_nodes1 = dtr2r1d_nodes_select(dtr2r1d_nodes1,ind);
    
    ind = 1 + floor((n_max-1)*linspace(0,1,n_nodes).^3);
    ind(ind<1) = 1;
    %dtr2d_nodes1(1)
    %ind
    if dtr2r1d_nodes1(1) == 0
        dtr2r1d = [];
        break
    end
    
    % repeat the dtr2r1d_nodes1 to generate a vector of n_nodes nodes, 
    % biased to contain more nodes corresponding to higher weights 
    dtr2r1d_nodes2 = dtr2r1d_nodes_select(dtr2r1d_nodes1,ind);
    % induce small mutations in the nodes within the dtr2r1d_nodes2 set
    dtr2r1d_nodes2 = dtr2r1d_nodes_mutate(dtr2r1d_nodes2,opt);
    % merge the original nodes with their mutations
    dtr2r1d_nodes = dtr2r1d_nodes_merge(dtr2r1d_nodes1,dtr2r1d_nodes2);
    
    % Fit to data
    dtr2r1d = dtr2r1d_data2dtr2r1d(stemp,bt_mx6,te,tr,dtr2r1d_nodes);    
end

if ~isempty(dtr2r1d)
    dtr2r1d_nodes = dtr2r1d_dist2nodes(dtr2r1d);
    n_max = min([opt.dtr2r1d.n_out dtr2r1d(1)]);
    ind = 1:n_max;
    dtr2r1d_nodes = dtr2r1d_nodes_select(dtr2r1d_nodes,ind);
    dtr2r1d = dtr2r1d_data2dtr2r1d(stemp,bt_mx6,te,tr,dtr2r1d_nodes);    
end
