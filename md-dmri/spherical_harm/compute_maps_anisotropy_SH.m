function compute_maps_anisotropy_SH(input_directory, output_directory)

if nargin < 1
    input_directory = pwd;
end

if nargin < 2
    output_directory = pwd;
end

[diso_SH_coeffs, h] = mdm_nii_read(fullfile(input_directory, 'ODF_SH_coeffs_diso.nii.gz'));
sqddelta_SH_coeffs = mdm_nii_read(fullfile(input_directory, 'ODF_SH_coeffs_sqddelta.nii.gz'));

squared_diso_SH_coeffs = abs(diso_SH_coeffs).^2;
squared_sqddelta_SH_coeffs = abs(sqddelta_SH_coeffs).^2;

aniso_SH_coeffs_diso_map = msf_notfinite2zero(sum(squared_diso_SH_coeffs(:,:,:,2:end),4)./sum(squared_diso_SH_coeffs,4));
aniso_SH_coeffs_sqddelta_map = msf_notfinite2zero(sum(squared_sqddelta_SH_coeffs(:,:,:,2:end),4)./sum(squared_sqddelta_SH_coeffs,4));

aniso_SH_coeffs_sqddelta_map(aniso_SH_coeffs_sqddelta_map > 0.05) = 0.05;

mdm_nii_write(aniso_SH_coeffs_diso_map, fullfile(output_directory, 'ODF_aniso_SH_coeffs_diso.nii.gz'), h);
mdm_nii_write(aniso_SH_coeffs_sqddelta_map, fullfile(output_directory, 'ODF_aniso_SH_coeffs_sqddelta.nii.gz'), h);