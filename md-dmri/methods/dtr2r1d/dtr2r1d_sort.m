function dtr2r1d_out = dtr2r1d_sort(dtr2r1d_in)

[~,par,perp,theta,phi,r2,r1,w] = dtr2r1d_dist2par(dtr2r1d_in);

[wsort,ind] = sort(w,'descend');
ind = ind(wsort>0); % select nodes with non-zero weight

n = numel(ind);
par = par(ind);
perp = perp(ind);
theta = theta(ind);
phi = phi(ind);
r2 = r2(ind);
r1 = r1(ind);
w = w(ind);

dtr2r1d_out = [par'; perp'; theta'; phi'; r2'; r1'; w'];
dtr2r1d_out = [n; dtr2r1d_out(:)];
