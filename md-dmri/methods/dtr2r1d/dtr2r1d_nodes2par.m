function [n,par,perp,theta,phi,r2,r1] = dtr2r1d_nodes2par(dtr2r1d_nodes)

if ~isempty(dtr2r1d_nodes)
    n = dtr2r1d_nodes(1);
    if n > 0
        m = numel(dtr2r1d_nodes(2:end))/n;
        dtr2r1d_nodes_array = reshape(dtr2r1d_nodes(2:end),[m n]);
        par = dtr2r1d_nodes_array(1,:)';
        perp = dtr2r1d_nodes_array(2,:)';
        theta = dtr2r1d_nodes_array(3,:)';
        phi = dtr2r1d_nodes_array(4,:)';
        r2 = dtr2r1d_nodes_array(5,:)';
        r1 = dtr2r1d_nodes_array(6,:)';
    else
        n = 0;
        par = [];
        perp = [];
        theta = [];
        phi = [];
        r2 = [];
        r1 = [];
    end
else
    n = 0;
    par = [];
    perp = [];
    theta = [];
    phi = [];
    r2 = [];
    r1 = [];
end
