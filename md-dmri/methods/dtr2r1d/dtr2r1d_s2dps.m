function dps = dtr2r1d_s2dps(dps, dtr2r1ds, nn)

%Per-voxel statistical measures
dps.s0 = sum(dtr2r1ds.w,4);

%Means
dps.mdiso = msf_notfinite2zero(sum(dtr2r1ds.diso.*dtr2r1ds.w,4)./dps.s0);
dps.msqdaniso = msf_notfinite2zero(sum(dtr2r1ds.sqdaniso.*dtr2r1ds.w,4)./dps.s0);
dps.msqddelta = msf_notfinite2zero(sum(dtr2r1ds.sqddelta.*dtr2r1ds.w,4)./dps.s0);
dps.mdxx = msf_notfinite2zero(sum(dtr2r1ds.dxx.*dtr2r1ds.w,4)./dps.s0);
dps.mdyy = msf_notfinite2zero(sum(dtr2r1ds.dyy.*dtr2r1ds.w,4)./dps.s0);
dps.mdzz = msf_notfinite2zero(sum(dtr2r1ds.dzz.*dtr2r1ds.w,4)./dps.s0);
dps.mdxy = msf_notfinite2zero(sum(dtr2r1ds.dxy.*dtr2r1ds.w,4)./dps.s0);
dps.mdxz = msf_notfinite2zero(sum(dtr2r1ds.dxz.*dtr2r1ds.w,4)./dps.s0);
dps.mdyz = msf_notfinite2zero(sum(dtr2r1ds.dyz.*dtr2r1ds.w,4)./dps.s0);
dps.mr2 = msf_notfinite2zero(sum(dtr2r1ds.r2.*dtr2r1ds.w,4)./dps.s0);
dps.mr1 = msf_notfinite2zero(sum(dtr2r1ds.r1.*dtr2r1ds.w,4)./dps.s0);

%Variances
dps.vdiso = msf_notfinite2zero(sum((dtr2r1ds.diso-repmat(dps.mdiso,[1 1 1 nn])).^2.*dtr2r1ds.w,4)./dps.s0);
dps.vsqdaniso = msf_notfinite2zero(sum((dtr2r1ds.sqdaniso-repmat(dps.msqdaniso,[1 1 1 nn])).^2.*dtr2r1ds.w,4)./dps.s0);
dps.vsqddelta = msf_notfinite2zero(sum((dtr2r1ds.sqddelta-repmat(dps.msqddelta,[1 1 1 nn])).^2.*dtr2r1ds.w,4)./dps.s0);
dps.vr2 = msf_notfinite2zero(sum((dtr2r1ds.r2-repmat(dps.mr2,[1 1 1 nn])).^2.*dtr2r1ds.w,4)./dps.s0);
dps.vr1 = msf_notfinite2zero(sum((dtr2r1ds.r1-repmat(dps.mr1,[1 1 1 nn])).^2.*dtr2r1ds.w,4)./dps.s0);

%Covariances
dps.cvdisosqdaniso = msf_notfinite2zero(sum((dtr2r1ds.diso-repmat(dps.mdiso,[1 1 1 nn])).*(dtr2r1ds.sqdaniso-repmat(dps.msqdaniso,[1 1 1 nn])).*dtr2r1ds.w,4)./dps.s0);
dps.cvdisosqddelta = msf_notfinite2zero(sum((dtr2r1ds.diso-repmat(dps.mdiso,[1 1 1 nn])).*(dtr2r1ds.sqddelta-repmat(dps.msqddelta,[1 1 1 nn])).*dtr2r1ds.w,4)./dps.s0);
dps.cvdisor2 = msf_notfinite2zero(sum((dtr2r1ds.diso-repmat(dps.mdiso,[1 1 1 nn])).*(dtr2r1ds.r2-repmat(dps.mr2,[1 1 1 nn])).*dtr2r1ds.w,4)./dps.s0);
dps.cvsqddeltar2 = msf_notfinite2zero(sum((dtr2r1ds.r2-repmat(dps.mr2,[1 1 1 nn])).*(dtr2r1ds.sqddelta-repmat(dps.msqddelta,[1 1 1 nn])).*dtr2r1ds.w,4)./dps.s0);
dps.cvdisor1 = msf_notfinite2zero(sum((dtr2r1ds.diso-repmat(dps.mdiso,[1 1 1 nn])).*(dtr2r1ds.r1-repmat(dps.mr1,[1 1 1 nn])).*dtr2r1ds.w,4)./dps.s0);
dps.cvsqddeltar1 = msf_notfinite2zero(sum((dtr2r1ds.r1-repmat(dps.mr1,[1 1 1 nn])).*(dtr2r1ds.sqddelta-repmat(dps.msqddelta,[1 1 1 nn])).*dtr2r1ds.w,4)./dps.s0);
dps.cvr2r1 = msf_notfinite2zero(sum((dtr2r1ds.r2-repmat(dps.mr2,[1 1 1 nn])).*(dtr2r1ds.r1-repmat(dps.mr1,[1 1 1 nn])).*dtr2r1ds.w,4)./dps.s0);

%Normalized measures
dps.vdison = msf_notfinite2zero(dps.vdiso./dps.mdiso.^2);
dps.vsqdanison = msf_notfinite2zero(dps.vsqdaniso./dps.msqdaniso.^2);
dps.vsqddeltan = msf_notfinite2zero(dps.vsqddelta./dps.msqddelta.^2);
dps.vr2n = msf_notfinite2zero(dps.vr2./dps.mr2.^2);
dps.vr1n = msf_notfinite2zero(dps.vr1./dps.mr1.^2);

end