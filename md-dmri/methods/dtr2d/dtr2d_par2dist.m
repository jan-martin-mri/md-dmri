function dtd = dtr2d_par2dist(par,perp,theta,phi,R2,w)

n = numel(par);

if n>0
    dtd = [par'; perp'; theta'; phi'; R2'; w'];
    dtd = [n; dtd(:)];
else
    dtd = [];
end

