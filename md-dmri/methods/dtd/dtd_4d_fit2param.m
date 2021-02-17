function dps = dtd_4d_fit2param(mfs_fn, dps_fn, opt)
% function dps = dtd_4d_fit2param(mfs_fn, dps_fn, opt)

if (nargin < 2), dps_fn = []; end
if (nargin < 3), opt = []; end

opt = mdm_opt(opt);
opt = dtd_opt(opt);
%dps = mdm_mfs_load(mfs_fn);

% Hack to allow mgui to access this function
if ischar(mfs_fn)
    dps = mdm_mfs_load(mfs_fn);
else
    m = mfs_fn;
    dps.m = m;
end

sz  = msf_size(dps.m(:,:,:,1), 3);

% reshape help functions
g_reshape = @(a,n) reshape(a, prod(sz(1:3)), n);
f_reshape = @(a,n) reshape(a, sz(1), sz(2), sz(3), n);

% create parameter maps and save them

m = dps.m;
dps = rmfield(dps,'m');
sz = size(m);

ind = false(sz(4),1);
ind(2:5:end) = 1;
nn = (sz(4)-1)/5;

dpar = m(:,:,:,circshift(ind,0,1));
dperp = m(:,:,:,circshift(ind,1,1));
theta = m(:,:,:,circshift(ind,2,1));
phi = m(:,:,:,circshift(ind,3,1));
w = m(:,:,:,circshift(ind,4,1));

%Calculate derived parameters
diso = (dpar + 2*dperp)/3;
daniso = (dpar - dperp)/3;
ddelta = msf_notfinite2zero(daniso./diso);
sdaniso = daniso.^2;
sddelta = msf_notfinite2zero(sdaniso./diso.^2);
dratio = msf_notfinite2zero(dpar./dperp);
[dxx,dyy,dzz,dxy,dxz,dyz] = dtd_pars2elements(dpar,dperp,theta,phi);

dtds = struct('w',w,'dpar',dpar,'dperp',dperp,'theta',theta,'phi',phi,'diso',diso,'daniso',daniso,'ddelta',ddelta,...
    'sdaniso',sdaniso,'sddelta',sddelta,'dratio',dratio,'dxx',dxx,'dyy',dyy,'dzz',dzz,'dxy',dxy,'dxz',dxz,'dyz',dyz);


    function dps = dtds2dps(dps, dtds)

        %Per-voxel statistical measures
        dps.s0 = sum(dtds.w,4);

        %Means
        dps.mdiso = msf_notfinite2zero(sum(dtds.diso.*dtds.w,4)./dps.s0);
        dps.mdaniso = msf_notfinite2zero(sum(dtds.daniso.*dtds.w,4)./dps.s0);
        dps.mddelta = msf_notfinite2zero(sum(dtds.ddelta.*dtds.w,4)./dps.s0);
        dps.msdaniso = msf_notfinite2zero(sum(dtds.sdaniso.*dtds.w,4)./dps.s0);
        dps.msddelta = msf_notfinite2zero(sum(dtds.sddelta.*dtds.w,4)./dps.s0);
        dps.mdxx = msf_notfinite2zero(sum(dtds.dxx.*dtds.w,4)./dps.s0);
        dps.mdyy = msf_notfinite2zero(sum(dtds.dyy.*dtds.w,4)./dps.s0);
        dps.mdzz = msf_notfinite2zero(sum(dtds.dzz.*dtds.w,4)./dps.s0);
        dps.mdxy = msf_notfinite2zero(sum(dtds.dxy.*dtds.w,4)./dps.s0);
        dps.mdxz = msf_notfinite2zero(sum(dtds.dxz.*dtds.w,4)./dps.s0);
        dps.mdyz = msf_notfinite2zero(sum(dtds.dyz.*dtds.w,4)./dps.s0);
        dps.t1x6 = cat(4, dps.mdxx, dps.mdyy, dps.mdzz, sqrt(2)*dps.mdxy, sqrt(2)*dps.mdxz, sqrt(2)*dps.mdyz);
        dps.maxmdii = max(cat(4, dps.mdxx, dps.mdyy, dps.mdzz),[],4);
        dps = tm_dt_to_dps(reshape(dps.t1x6,[prod(sz(1:3)) 6]), dps, f_reshape);

        %Variances
        dps.vdiso = msf_notfinite2zero(sum((dtds.diso-repmat(dps.mdiso,[1 1 1 nn])).^2.*dtds.w,4)./dps.s0);
        dps.vsdaniso = msf_notfinite2zero(sum((dtds.sdaniso-repmat(dps.msdaniso,[1 1 1 nn])).^2.*dtds.w,4)./dps.s0);
        dps.vsddelta = msf_notfinite2zero(sum((dtds.sddelta-repmat(dps.msddelta,[1 1 1 nn])).^2.*dtds.w,4)./dps.s0);

        %Covariances
        dps.cvdisosdaniso = msf_notfinite2zero(sum((dtds.diso-repmat(dps.mdiso,[1 1 1 nn])).*(dtds.sdaniso-repmat(dps.msdaniso,[1 1 1 nn])).*dtds.w,4)./dps.s0);
        dps.cvdisosddelta = msf_notfinite2zero(sum((dtds.diso-repmat(dps.mdiso,[1 1 1 nn])).*(dtds.sddelta-repmat(dps.msddelta,[1 1 1 nn])).*dtds.w,4)./dps.s0);

        %Normalized measures
        dps.mdanison = msf_notfinite2zero(dps.mdaniso./dps.mdiso);
        dps.msdanison = msf_notfinite2zero(dps.msdaniso./dps.mdiso.^2);
        dps.vdison = msf_notfinite2zero(dps.vdiso./dps.mdiso.^2);
        dps.vsdanison = msf_notfinite2zero(dps.vsdaniso./dps.mdiso.^4);
%         dps.vsdanison = msf_notfinite2zero(dps.vsdaniso./dps.msdaniso.^2);
        dps.vsddeltan = msf_notfinite2zero(dps.vsddelta./dps.msddelta.^2);
        
        dps.Vl = 5/2 * 4/5*dps.msdaniso;
        dps.MKi = 3 * dps.vdison; % Multiply by 3 to get kurtosis
        dps.MKa = 3 * 4/5*dps.msdanison;

        % Calculate uFA. Take real component to avoid complex values due to
        % sqrt of negative variances.
        dps.ufa_old = real(sqrt(3/2) * sqrt(1./(dps.mdiso.^2./dps.Vl+1))); % Lasic (2014)
        dps.ufa     = real(sqrt(3/2) * sqrt( dps.Vl ./ (dps.Vl + dps.vdiso + dps.mdiso.^2) )); % Szczepankiewicz (2016)
        dps.uFA = dps.ufa;
        dps.MD = dps.mdiso*1e9;
        
        dps.signaniso = sign(dps.mdanison);
    end

dps = dtds2dps(dps, dtds);

%Per-bin statistical measures
for nbin = 1:numel(opt.dtd.bin_disomax)
    ind_bin = false([sz(1) sz(2) sz(3) nn 4]);
    ind_bin(:,:,:,:,1) = diso >= opt.dtd.bin_disomin(nbin);
    ind_bin(:,:,:,:,2) = diso <= opt.dtd.bin_disomax(nbin);
    ind_bin(:,:,:,:,3) = dratio >= opt.dtd.bin_dratiomin(nbin);
    ind_bin(:,:,:,:,4) = dratio <= opt.dtd.bin_dratiomax(nbin);
    ind_bin = all(ind_bin,5);

    dps_bin.no = nbin;
    dtds_temp = dtds;
    dtds_temp.w = dtds.w.*ind_bin;
    dps_bin = dtds2dps(dps_bin, dtds_temp);
    dps.bin{nbin} = dps_bin;
end

% %DTI params
% dps.lambdazzvec = zeros([sz(1) sz(2) sz(3) 3]);
% dps.lambdaxxvec = zeros([sz(1) sz(2) sz(3) 3]);
% dps.lambdayyvec = zeros([sz(1) sz(2) sz(3) 3]);
% dps.lambda11vec = zeros([sz(1) sz(2) sz(3) 3]);
% dps.lambda22vec = zeros([sz(1) sz(2) sz(3) 3]);
% dps.lambda33vec = zeros([sz(1) sz(2) sz(3) 3]);
% dtiparam = {'trace','iso','lambda33','lambda22','lambda11','lambdazz','lambdaxx','lambdayy','vlambda',...
%     'delta','eta','s','p','l','fa','cs','cp','cl','cm'};
% param = {dtiparam{:}};
% for nparam = 1:numel(param)
%     eval(['dps.' param{nparam} ' = zeros([sz(1) sz(2) sz(3)]);']);
% end
% 
% Nbin = numel(opt.dtd.bin_isomax);
% for nbin = 1:Nbin
%     eval(['dps.bin' num2str(nbin) ' = zeros([sz(1) sz(2) sz(3)]);']);
% end
% 
% for nk = 1:sz(3)
%     for nj = 1:sz(2)
%         for ni = 1:sz(1)
%             %ni = 11; nj = 9; nk = 1;
%             if dps.mask(ni,nj,nk)
% %                     for nbin = 1:Nbin
% %                         eval(['ind = all([iso_v > opt.dtd.bin_isomin(' num2str(nbin) '), iso_v < opt.dtd.bin_isomax(' num2str(nbin) '), ratio_v > opt.dtd.bin_ratiomin(' num2str(nbin) '), ratio_v < opt.dtd.bin_ratiomax(' num2str(nbin) ')],2);'])
% %                         bin = ind'*w;
% %                         eval(['dps.bin' num2str(nbin) '(ni,nj,nk) = bin;']);
% %                     end
% 
%                 dt1x6 = [dps.mdxx(ni,nj,nk) dps.mdyy(ni,nj,nk) dps.mdzz(ni,nj,nk) sqrt(2)*[dps.mdxy(ni,nj,nk) dps.mdxz(ni,nj,nk) dps.mdyz(ni,nj,nk)]];
%                 dt3x3 = tm_1x6_to_3x3(dt1x6);
% 
%                 dt = tm_3x3_to_tpars(dt3x3);
% 
%                 dps.t1x6(ni,nj,nk,:) = dt.t1x6;
%                 dps.lambdazzvec(ni,nj,nk,:) = dt.lambdazzvec;
%                 dps.lambdaxxvec(ni,nj,nk,:) = dt.lambdaxxvec;
%                 dps.lambdayyvec(ni,nj,nk,:) = dt.lambdayyvec;
%                 dps.lambda11vec(ni,nj,nk,:) = dt.lambda11vec;
%                 dps.lambda22vec(ni,nj,nk,:) = dt.lambda22vec;
%                 dps.lambda33vec(ni,nj,nk,:) = dt.lambda33vec;
% 
%                 for nparam = 1:numel(dtiparam)
%                     eval(['dps.' dtiparam{nparam} '(ni,nj,nk) = dt.' dtiparam{nparam} ';']);
%                 end
%                     
%                                 
%             end
%         end
%     end
% end

% %2nd moments of P(Deff)
% dps.mu2iso = dps.vdiso;
% dps.mu2aniso = 4/5*dps.msdaniso;
% dps.mu2macro = 4/45*((dps.lambda33-dps.lambda11).^2 + (dps.lambda22-dps.lambda11).*(dps.lambda22-dps.lambda33));
% 
% % dps.ufa = sqrt(3/2)*sqrt(1./((dps.miso.^2+dps.mu2iso)./(5/2*dps.mu2aniso)+1));
% % dps.ufa(isnan(dps.ufa)) = 0;
% % dps.ufa(dps.ufa>1) = 1;
% 
% kronecker = permute(repmat([1 1 1 0 0 0]',[1 size(dps.s0,3) size(dps.s0,2) size(dps.s0,1)]),[4 3 2 1]);
% dps.s1x6 = (dps.t1x6./repmat(dps.miso,[1 1 1 6]) - kronecker)./repmat(dps.mdelta,[1 1 1 6])/2;
% dps.s1x6prim = (2*dps.s1x6 + kronecker)/3;
% dps.s1x6(isnan(dps.s1x6)) = 0;
% dps.s1x6prim(isnan(dps.s1x6prim)) = 0;
% 
% dps.slambdaxx = (dps.lambdaxx./dps.miso - 1)./dps.mdelta/2;
% dps.slambdaxx(isnan(dps.slambdaxx)) = 0;
% dps.slambdayy = (dps.lambdayy./dps.miso - 1)./dps.mdelta/2;
% dps.slambdayy(isnan(dps.slambdayy)) = 0;
% dps.slambdazz = (dps.lambdazz./dps.miso - 1)./dps.mdelta/2;
% dps.slambdazz(isnan(dps.slambdazz)) = 0;
% 
% dps.slambdaxxprim = (2*dps.slambdaxx + 1)/3;
% dps.slambdayyprim = (2*dps.slambdayy + 1)/3;
% dps.slambdazzprim = (2*dps.slambdazz + 1)/3;
% 
% dps.slambdas = zeros([sz(1) sz(2) sz(3) 3]);
% dps.slambdas(:,:,:,1) = dps.slambdaxxprim;
% dps.slambdas(:,:,:,2) = dps.slambdayyprim;
% dps.slambdas(:,:,:,3) = dps.slambdazzprim;
% dps.slambda11prim = min(dps.slambdas,[],4);
% dps.slambda33prim = max(dps.slambdas,[],4);
% dps = rmfield(dps,'slambdas');
% 
% dps.op = sqrt(dps.mu2macro./dps.mu2aniso);
% dps.op(dps.op>1) = 1;
% dps.op(isnan(dps.op)) = 0;

if (~isempty(dps_fn)) mdm_dps_save(dps, dps.s, dps_fn, opt); end

end

