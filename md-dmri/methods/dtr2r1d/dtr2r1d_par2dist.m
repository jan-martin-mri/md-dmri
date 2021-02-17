function dtr2r1d = dtr2r1d_par2dist(par,perp,theta,phi,r2,r1,w)

n = numel(par);

if n>0
    dtr2r1d = [par'; perp'; theta'; phi'; r2'; r1'; w'];
    dtr2r1d = [n; dtr2r1d(:)];
else
    dtr2r1d = [];
end