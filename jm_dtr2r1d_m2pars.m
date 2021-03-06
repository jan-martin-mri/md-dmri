function [dpar, dperp, theta, phi, r2, r1, w] = jm_dtr2r1d_m2pars( m )

sz = size(m);

ind = false(sz(4), 1);
ind(2:7:end) = 1;

dpar  = m(:,:,:,ind);
dperp = m(:,:,:,circshift(ind,1,1));
theta = m(:,:,:,circshift(ind,2,1));
phi   = m(:,:,:,circshift(ind,3,1));
r2    = m(:,:,:,circshift(ind,4,1));
r1    = m(:,:,:,circshift(ind,5,1));
w     = m(:,:,:,circshift(ind,6,1));

end