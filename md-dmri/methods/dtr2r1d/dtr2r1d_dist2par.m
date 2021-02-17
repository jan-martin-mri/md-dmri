function [n,par,perp,theta,phi,r2,r1,w] = dtr2r1d_dist2par(dtr2r1d)

n = dtr2r1d(1);

if n>0
    m = numel(dtr2r1d(2:end))/n;
    dtr2r1d_array = reshape(dtr2r1d(2:end),[m n]);
    par = dtr2r1d_array(1,:)';
    perp = dtr2r1d_array(2,:)';
    theta = dtr2r1d_array(3,:)';
    phi = dtr2r1d_array(4,:)';
    r2 = dtr2r1d_array(5,:)';
    r1 = dtr2r1d_array(6,:)';
    w = dtr2r1d_array(7,:)';
else
    par = [];
    perp = [];
    theta = [];
    phi = [];
    r2 = [];
    r1 = [];
    w = [];
end
