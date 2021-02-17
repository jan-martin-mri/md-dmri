function odf = dtr2d_4d_fit2podf(mfs_fn, opt)
% function odf = dtr2d_4d_fit2odf(mfs_fn, opt)

if (nargin < 2)
    opt = mdm_opt(opt);
    opt = dtr2d_opt(opt);
end

% create parameter maps and save them
mfs = mdm_mfs_load(mfs_fn);
m = mfs.m;
sz = size(m);

ind = false(sz(4),1);
ind(2:6:end) = 1;
nn = (sz(4)-1)/6;

dpar = m(:,:,:,circshift(ind,0,1));
dperp = m(:,:,:,circshift(ind,1,1));
theta = m(:,:,:,circshift(ind,2,1));
phi = m(:,:,:,circshift(ind,3,1));
r2 = m(:,:,:,circshift(ind,4,1));
w = m(:,:,:,circshift(ind,5,1));
if isfield(opt.dtr2d,'r2extrap') == 1
    w = w.*exp(-r2*opt.dtr2d.r2extrap);
end

%Calculate derived parameters
diso = (dpar + 2*dperp)/3;
daniso = (dpar - dperp)/3;
ddelta = msf_notfinite2zero(daniso./diso);
sqdaniso = daniso.^2;
sqddelta = msf_notfinite2zero(sqdaniso./diso.^2);
dratio = msf_notfinite2zero(dpar./dperp);
[dxx,dyy,dzz,dxy,dxz,dyz] = dtr2d_pars2elements(dpar,dperp,theta,phi);

dtds = struct('w',w,'dpar',dpar,'dperp',dperp,'theta',theta,'phi',phi,'diso',diso,'daniso',daniso,'ddelta',ddelta,...
    'sqdaniso',sqdaniso,'sqddelta',sqddelta,'dratio',dratio,'dxx',dxx,'dyy',dyy,'dzz',dzz,'dxy',dxy,'dxz',dxz,'dyz',dyz,'r2',r2);

%ODF nodes
run_path = mfilename('fullpath');
framework_path = fileparts(fileparts(fileparts(run_path)));
angles_path = fullfile(framework_path,'tools','uvec','repulsion_angles_tri');

odf_s.n = opt.dtr2d.odf_nnodes; %250, 350, 500, 1000, 3994, or 15970
angles = load(fullfile(angles_path,num2str(odf_s.n)));
odf_s.x = sin(angles.theta).*cos(angles.phi);
odf_s.y = sin(angles.theta).*sin(angles.phi);
odf_s.z = cos(angles.theta);
odf_s.tri = angles.tri;
ODindex = .05; %Watson distribution smoothing kernel
odf_s.kappa = 1/tan(ODindex*pi/2);

odf = odf_s;
%odf_w = dtds2odf(odf_s, dtds);
%odf.w = odf_w;

%Per-bin ODFs
for nbin = 1:numel(opt.dtr2d.bin_disomax)
    ind_bin = false([sz(1) sz(2) sz(3) nn 4]);
    ind_bin(:,:,:,:,1) = diso >= opt.dtr2d.bin_disomin(nbin);
    ind_bin(:,:,:,:,2) = diso <= opt.dtr2d.bin_disomax(nbin);
    ind_bin(:,:,:,:,3) = sqddelta >= opt.dtr2d.bin_sddeltamin(nbin);
    ind_bin(:,:,:,:,4) = sqddelta <= opt.dtr2d.bin_sddeltamax(nbin);
    ind_bin = all(ind_bin,5);
    
    

    dtds_temp = dtds;
    dtds_temp.w = dtds.w.*ind_bin;
    
    odf_w = zeros(sz(1), sz(2), sz(3), odf_s.n);
    odf_diso = zeros(sz(1), sz(2), sz(3), odf_s.n);
    odf_sqddelta = zeros(sz(1), sz(2), sz(3), odf_s.n);
    odf_r2 = zeros(sz(1), sz(2), sz(3), odf_s.n);
    
    parfor nk = 1:sz(3)
        list_w = zeros(sz(1), sz(2), odf_s.n);
        list_diso = zeros(sz(1), sz(2), odf_s.n);
        list_sqddelta = zeros(sz(1), sz(2), odf_s.n);
        list_r2 = zeros(sz(1), sz(2), odf_s.n);
        for nj = 1:sz(2)
            for ni = 1:sz(1)
                if mfs.mask(ni,nj,nk)
                    theta_vox = squeeze(dtds_temp.theta(ni,nj,nk,:));
                    phi_vox = squeeze(dtds_temp.phi(ni,nj,nk,:));
                    w_vox = squeeze(dtds_temp.w(ni,nj,nk,:));
                    diso_vox = squeeze(dtds_temp.diso(ni,nj,nk,:));
                    sqddelta_vox = squeeze(dtds_temp.sqddelta(ni,nj,nk,:));
                    r2_vox = squeeze(dtds_temp.r2(ni,nj,nk,:));
                    odf_d = struct;
                    odf_d.x = sin(theta_vox).*cos(phi_vox);
                    odf_d.y = sin(theta_vox).*sin(phi_vox);
                    odf_d.z = cos(theta_vox);
                    odf_d.w = w_vox;
                    odf_d.diso = diso_vox;
                    odf_d.sqddelta = sqddelta_vox;
                    odf_d.r2 = r2_vox;
                    odf_vox = dist_podf_discrete2smooth(odf_d,odf_s);
                    
                    list_w(ni,nj,:) = odf_vox.w(:);
                    list_diso(ni,nj,:) = odf_vox.diso(:);
                    list_sqddelta(ni,nj,:) = odf_vox.sqddelta(:);
                    list_r2(ni,nj,:) = odf_vox.r2(:);
                    
%                     odf_w(ni,nj,nk,:) = odf_vox.w(:);
%                     odf_diso(ni,nj,nk,:) = odf_vox.diso(:);
%                     odf_sqddelta(ni,nj,nk,:) = odf_vox.sqddelta(:);
%                     odf_r2(ni,nj,nk,:) = odf_vox.r2(:);
                    
                end
            end
        end
        odf_w(:,:,nk,:) = list_w;
        odf_diso(:,:,nk,:) = list_diso;
        odf_sqddelta(:,:,nk,:) = list_sqddelta;
        odf_r2(:,:,nk,:) = list_r2;
    end
    
    odf.w_bin{nbin} = odf_w;
    odf.diso_bin{nbin} = odf_diso;
    odf.sqddelta_bin{nbin} = odf_sqddelta;
    odf.r2_bin{nbin} = odf_r2;
end

end