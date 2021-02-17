function dps = dtd_mv_gamma_4d_fit2param(mfs_fn, dps_fn, opt)
% function dps = dtd_mv_gamma_4d_fit2param(mfs_fn, dps_fn, opt)
%
% Calculates derived parameters from the primary parameters of the gamma fit.

if (nargin < 2), dps_fn = []; end
if (nargin < 3), opt = []; end

opt = mdm_opt(opt);

%dps = mdm_mfs_load(mfs_fn);

% Hack to allow mgui to access this function
if ischar(mfs_fn)
    dps = mdm_mfs_load(mfs_fn);
else
    m = mfs_fn;
    dps.m = m;
end

s0 = dps.m(:,:,:,1);
kappa = dps.m(:,:,:,2);
psi_1 = dps.m(:,:,:,3);
psi_2 = dps.m(:,:,:,4);
psi_3 = dps.m(:,:,:,5);
h_1   = dps.m(:,:,:,6);
h_2   = dps.m(:,:,:,7);
h_3   = dps.m(:,:,:,8);
alpha = dps.m(:,:,:,9);
beta  = dps.m(:,:,:,10);
gamma = dps.m(:,:,:,11);

Psi = cat(4,psi_1,psi_2,psi_3);
average_D_1 = psi_1./h_1;
average_D_2 = psi_2./h_2;
average_D_3 = psi_3./h_3;
average_D = cat(4,average_D_1,average_D_2,average_D_3);
average_D = msf_notfinite2zero(average_D);

mdiso = sum(average_D,4)/3;
var_eig = var(average_D,1,4)/3;
FA = real(sqrt((3/2)*(1 + mdiso.^2./var_eig).^(-1)));

% [mddelta, mdeta, ~] = dtd_mv_gamma_multidim_Haeberlen(average_D, mdiso);
% mdxx_Haeberlen = mdiso.*(1-mddelta.*(1-mdeta));
% mdyy_Haeberlen = mdiso.*(1-mddelta.*(1+mdeta));
% mdzz_Haeberlen = mdiso.*(1+2.*mddelta);

[mdxx,mdyy,mdzz,mdxy,mdxz,mdyz] = dtd_mv_gamma_multidim_Euler_rotation(average_D,alpha,beta,gamma);

average_D_squared = dtd_mv_gamma_multidim_outer_product(average_D,average_D);
[E_bulk, E_shear, ~] = tm_6x6_iso();
E_bulk = E_bulk(1:3,1:3);
E_shear = E_shear(1:3,1:3);

covariance = dtd_mv_gamma_multidim_outer_product(average_D,Psi)...
             + dtd_mv_gamma_multidim_outer_product(Psi,average_D) ...
             - kappa.*dtd_mv_gamma_multidim_outer_product(Psi,Psi);
covariance = msf_notfinite2zero(covariance);
covariance = dtd_mv_gamma_multidim_ensure_symmetry(covariance); % Make sure that the covariance is numerically symmetric
covariance = dtd_mv_gamma_multidim_nearest_SPD(covariance); % Project to nearest symmetric positive-definite matrix

V_shear2 = dtd_mv_gamma_multidim_inner_product(average_D_squared, E_shear);
V_shear  = dtd_mv_gamma_multidim_inner_product(covariance, E_shear);
V_shear1 = V_shear + V_shear2;

vdiso = dtd_mv_gamma_multidim_inner_product(covariance, E_bulk);
msdaniso = V_shear1/2;
sOP = real(sqrt(V_shear2./V_shear1));

dps.s0 = msf_notfinite2zero(s0);
dps.mdxx = msf_notfinite2zero(mdxx);
dps.mdyy = msf_notfinite2zero(mdyy);
dps.mdzz = msf_notfinite2zero(mdzz);
dps.mdxy = msf_notfinite2zero(mdxy);
dps.mdxz = msf_notfinite2zero(mdxz);
dps.mdyz = msf_notfinite2zero(mdyz);
dps.t1x6 = cat(4, dps.mdxx, dps.mdyy, dps.mdzz, sqrt(2)*dps.mdxy, sqrt(2)*dps.mdxz, sqrt(2)*dps.mdyz);
dps.maxmdii = max(cat(4, dps.mdxx, dps.mdyy, dps.mdzz),[],4);


if ~isfield(opt,'dtd_mv_gamma')
    opt = dtd_mv_gamma_opt(opt);
end
dps.mdiso = dtd_mv_gamma_clamp(msf_notfinite2zero(mdiso),[opt.dtd_mv_gamma.diso_min opt.dtd_mv_gamma.diso_max]);
dps.FA = dtd_mv_gamma_clamp(msf_notfinite2zero(FA),[opt.dtd_mv_gamma.FA_min opt.dtd_mv_gamma.FA_max]);
dps.vdiso = dtd_mv_gamma_clamp(msf_notfinite2zero(vdiso),[opt.dtd_mv_gamma.vdiso_min opt.dtd_mv_gamma.vdiso_max]);
dps.msdaniso = dtd_mv_gamma_clamp(msf_notfinite2zero(msdaniso),[opt.dtd_mv_gamma.msdaniso_min opt.dtd_mv_gamma.msdaniso_max]);
dps.vdison = dtd_mv_gamma_clamp(msf_notfinite2zero(dps.vdiso./dps.mdiso.^2),[opt.dtd_mv_gamma.vdison_min opt.dtd_mv_gamma.vdison_max]); % Normalized
dps.msdanison = dtd_mv_gamma_clamp(msf_notfinite2zero(dps.msdaniso./dps.mdiso.^2),[opt.dtd_mv_gamma.msdanison_min opt.dtd_mv_gamma.msdanison_max]);
dps.sOP = dtd_mv_gamma_clamp(msf_notfinite2zero(sOP),[opt.dtd_mv_gamma.sOP_min opt.dtd_mv_gamma.sOP_max]);

% dps.mddelta = dtd_mv_gamma_clamp(msf_notfinite2zero(mddelta),[-0.5 1]);
% dps.msddelta = dtd_mv_gamma_clamp(msf_notfinite2zero(mddelta.^2),[0 1]);
% dps.mdeta = msf_notfinite2zero(mdeta);
% dps.msdeta = msf_notfinite2zero(mdeta.^2);
% dps.mdxx_Haeberlen = msf_notfinite2zero(mdxx_Haeberlen);
% dps.mdyy_Haeberlen = msf_notfinite2zero(mdyy_Haeberlen);
% dps.mdzz_Haeberlen = msf_notfinite2zero(mdzz_Haeberlen);

% for i = 5:size(dps.m, 4)
%     nam = ['s' num2str(i-4)];
%     dps.(nam) = dps.m(:,:,:,i);
% end

if (~isempty(dps_fn)), mdm_dps_save(dps, dps.s, dps_fn, opt); end
