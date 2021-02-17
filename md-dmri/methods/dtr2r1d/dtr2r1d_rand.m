function dtr2r1d_nodes = dtr2r1d_rand(n,dmin,dmax,r2min,r2max,r1min,r1max)


% ratio = 1*(100).^rand(1,n);
% delta = (ratio-1)./(2+ratio);
% iso = dmax./(1+2*delta);
% iso = dmin*(iso/dmin).^rand(1,n);

par = dmin*(dmax/dmin).^rand(1,n);
perp = dmin*(dmax/dmin).^rand(1,n);
perp(round(.9*n):n) = par(round(.9*n):n);
% perp = dmin*(dmax/dmin).^rand(1,n);

% par = iso.*(1+2*delta);
% perp = iso.*(1-delta);

% Good ones to switch back to, replace line 10
% perp = dmin*(par/dmin).^rand(1,n);
% perp(round(.9*n):n) = par(round(.9*n):n);

theta = acos(2*rand(1,n)-1);
phi = 2*pi*rand(1,n);
r2 = r2min*(r2max/r2min).^rand(1,n);
r1 = r1min*(r1max/r1min).^rand(1,n);

dtr2r1d_nodes = [par; perp; theta; phi; r2; r1];
dtr2r1d_nodes = [n; dtr2r1d_nodes(:)];


