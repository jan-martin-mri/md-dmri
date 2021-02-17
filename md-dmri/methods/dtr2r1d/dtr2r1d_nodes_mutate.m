function dtr2r1d_nodes_out = dtr2r1d_nodes_mutate(dtr2r1d_nodes,opt)

[n,par,perp,theta,phi,r2,r1] = dtr2r1d_nodes2par(dtr2r1d_nodes);

par = par.*(1+opt.dtr2r1d.dfuzz*randn(n,1));
perp = perp.*(1+opt.dtr2r1d.dfuzz*randn(n,1));
theta = theta + opt.dtr2r1d.ofuzz*randn(n,1);
phi = phi + opt.dtr2r1d.ofuzz*randn(n,1);
r2 = r2.*(1+opt.dtr2r1d.r2fuzz*randn(n,1));
r1 = r1.*(1+opt.dtr2r1d.r1fuzz*randn(n,1));

par(par>opt.dtr2r1d.dmax) = opt.dtr2r1d.dmax;
perp(perp>opt.dtr2r1d.dmax) = opt.dtr2r1d.dmax;
r2(r2>opt.dtr2r1d.r2max) = opt.dtr2r1d.r2max;
r1(r1>opt.dtr2r1d.r1max) = opt.dtr2r1d.r1max;
par(par<opt.dtr2r1d.dmin) = opt.dtr2r1d.dmin;
perp(perp<opt.dtr2r1d.dmin) = opt.dtr2r1d.dmin;
r2(r2<opt.dtr2r1d.r2min) = opt.dtr2r1d.r2min;
r1(r1<opt.dtr2r1d.r1min) = opt.dtr2r1d.r1min;

dtr2r1d_nodes_out = [par'; perp'; theta'; phi'; r2'; r1'];
dtr2r1d_nodes_out = [n; dtr2r1d_nodes_out(:)];
