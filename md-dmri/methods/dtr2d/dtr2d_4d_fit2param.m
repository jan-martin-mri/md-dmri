function dps = dtr2d_4d_fit2param(mfs_fn, dps_fn, opt)
% function dps = dtr2d_4d_fit2param(mfs_fn, dps_fn, opt)

if (nargin < 2), dps_fn = []; end
if (nargin < 3), opt = []; end

opt = mdm_opt(opt);
opt = dtr2d_opt(opt);
dps = mdm_mfs_load(mfs_fn);

% create parameter maps and save them

m = dps.m;
dps = rmfield(dps,'m');
sz = size(m);

ind = false(sz(4),1);
ind(2:6:end) = 1;
nn = (sz(4)-1)/6; %Number of solution nodes

dpar = m(:,:,:,circshift(ind,0,1));
dperp = m(:,:,:,circshift(ind,1,1));
theta = m(:,:,:,circshift(ind,2,1));
phi = m(:,:,:,circshift(ind,3,1));
r2 = m(:,:,:,circshift(ind,4,1));
w = m(:,:,:,circshift(ind,5,1));

%Calculate derived parameters
diso = (dpar + 2*dperp)/3;
daniso = (dpar - dperp)/3;
ddelta = msf_notfinite2zero(daniso./diso);
sqdaniso = daniso.^2;
sqddelta = msf_notfinite2zero(sqdaniso./diso.^2);
dratio = msf_notfinite2zero(dpar./dperp);
[dxx,dyy,dzz,dxy,dxz,dyz] = dtr2d_pars2elements(dpar,dperp,theta,phi);

dtr2ds = struct('w',w,'dpar',dpar,'dperp',dperp,'theta',theta,'phi',phi,...
    'diso',diso,'daniso',daniso,'ddelta',ddelta,'sqdaniso',sqdaniso,...
    'sqddelta',sqddelta,'dratio',dratio,'dxx',dxx,'dyy',dyy,'dzz',dzz,...
    'dxy',dxy,'dxz',dxz,'dyz',dyz,'r2',r2);

dps = dtr2d_s2dps(dps, dtr2ds, nn);

%Per-bin statistical measures
if opt.dtr2d.do_bin == 1
    for nbin = 1:numel(opt.dtr2d.bin_disomax)
        ind_bin = false([sz(1) sz(2) sz(3) nn 4]);
        ind_bin(:,:,:,:,1) = r2 >= opt.dtr2d.bin_r2min(nbin);
        ind_bin(:,:,:,:,2) = r2 <= opt.dtr2d.bin_r2max(nbin);
        ind_bin(:,:,:,:,3) = diso >= opt.dtr2d.bin_disomin(nbin);
        ind_bin(:,:,:,:,4) = diso <= opt.dtr2d.bin_disomax(nbin);
        ind_bin(:,:,:,:,5) = dratio >= opt.dtr2d.bin_dratiomin(nbin);
        ind_bin(:,:,:,:,6) = dratio <= opt.dtr2d.bin_dratiomax(nbin);
        
        ind_bin = all(ind_bin,5);
        
        dps_bin.no = nbin;
        dtr2ds_temp = dtr2ds;
        dtr2ds_temp.w = dtr2ds.w.*ind_bin;
        dps_bin = dtr2d_s2dps(dps_bin, dtr2ds_temp, nn);
        dps.bin{nbin} = dps_bin;
    end
end

if (~isempty(dps_fn)) mdm_dps_save(dps, dps.s, dps_fn, opt); end

end

