function dps = dtr2d_s2dps(dps, dtr2ds, nn)

%Per-voxel statistical measures
dps.s0 = sum(dtr2ds.w,4);

%Means
dps.mdiso = msf_notfinite2zero(sum(dtr2ds.diso.*dtr2ds.w,4)./dps.s0);
dps.msqdaniso = msf_notfinite2zero(sum(dtr2ds.sqdaniso.*dtr2ds.w,4)./dps.s0);
dps.msqddelta = msf_notfinite2zero(sum(dtr2ds.sqddelta.*dtr2ds.w,4)./dps.s0);
dps.mdxx = msf_notfinite2zero(sum(dtr2ds.dxx.*dtr2ds.w,4)./dps.s0);
dps.mdyy = msf_notfinite2zero(sum(dtr2ds.dyy.*dtr2ds.w,4)./dps.s0);
dps.mdzz = msf_notfinite2zero(sum(dtr2ds.dzz.*dtr2ds.w,4)./dps.s0);
dps.mdxy = msf_notfinite2zero(sum(dtr2ds.dxy.*dtr2ds.w,4)./dps.s0);
dps.mdxz = msf_notfinite2zero(sum(dtr2ds.dxz.*dtr2ds.w,4)./dps.s0);
dps.mdyz = msf_notfinite2zero(sum(dtr2ds.dyz.*dtr2ds.w,4)./dps.s0);
dps.t1x6 = [dps.mdxx dps.mdyy dps.mdzz sqrt(2)*[dps.mdxy dps.mdxz dps.mdyz]];
dps.mr2 = msf_notfinite2zero(sum(dtr2ds.r2.*dtr2ds.w,4)./dps.s0);

%Variances
dps.vdiso = msf_notfinite2zero(sum((dtr2ds.diso-repmat(dps.mdiso,[1 1 1 nn])).^2.*dtr2ds.w,4)./dps.s0);
dps.vsqdaniso = msf_notfinite2zero(sum((dtr2ds.sqdaniso-repmat(dps.msqdaniso,[1 1 1 nn])).^2.*dtr2ds.w,4)./dps.s0);
dps.vsqddelta = msf_notfinite2zero(sum((dtr2ds.sqddelta-repmat(dps.msqddelta,[1 1 1 nn])).^2.*dtr2ds.w,4)./dps.s0);
dps.vr2 = msf_notfinite2zero(sum((dtr2ds.r2-repmat(dps.mr2,[1 1 1 nn])).^2.*dtr2ds.w,4)./dps.s0);

%Covariances
dps.cvdisosqdaniso = msf_notfinite2zero(sum((dtr2ds.diso-repmat(dps.mdiso,[1 1 1 nn])).*(dtr2ds.sqdaniso-repmat(dps.msqdaniso,[1 1 1 nn])).*dtr2ds.w,4)./dps.s0);
dps.cvdisosqddelta = msf_notfinite2zero(sum((dtr2ds.diso-repmat(dps.mdiso,[1 1 1 nn])).*(dtr2ds.sqddelta-repmat(dps.msqddelta,[1 1 1 nn])).*dtr2ds.w,4)./dps.s0);
dps.cvdisor2 = msf_notfinite2zero(sum((dtr2ds.diso-repmat(dps.mdiso,[1 1 1 nn])).*(dtr2ds.r2-repmat(dps.mr2,[1 1 1 nn])).*dtr2ds.w,4)./dps.s0);
dps.cvsqddeltar2 = msf_notfinite2zero(sum((dtr2ds.r2-repmat(dps.mr2,[1 1 1 nn])).*(dtr2ds.sqddelta-repmat(dps.msqddelta,[1 1 1 nn])).*dtr2ds.w,4)./dps.s0);

%Normalized measures
dps.vdison = msf_notfinite2zero(dps.vdiso./dps.mdiso.^2);
dps.vsqdanison = msf_notfinite2zero(dps.vsqdaniso./dps.msqdaniso.^2);
dps.vsqddeltan = msf_notfinite2zero(dps.vsqddelta./dps.msqddelta.^2);
dps.vr2n = msf_notfinite2zero(dps.vr2./dps.mr2.^2);


end