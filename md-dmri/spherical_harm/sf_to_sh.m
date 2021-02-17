function sf_to_sh(struct_sf_to_sh, SH_order, flip, mask_file, output_directory, method)
% Spherical function to spherical harmonics

FLAG_relaxation = false;
if strcmp(method, 'dtr2d') || strcmp(method, 'dtr1d')
    FLAG_relaxation = true;
end

FLAG_double_relaxation = strcmp(method, 'dtr2r1d');

x = struct_sf_to_sh.x;
y = struct_sf_to_sh.y;
z = struct_sf_to_sh.z;
w = struct_sf_to_sh.w;
diso = struct_sf_to_sh.diso;
sqddelta = struct_sf_to_sh.sqddelta;
if FLAG_relaxation
    r = struct_sf_to_sh.r;
    t = struct_sf_to_sh.t;
elseif FLAG_double_relaxation
    r1 = struct_sf_to_sh.r1;
    t1 = struct_sf_to_sh.t1;
    r2 = struct_sf_to_sh.r2;
    t2 = struct_sf_to_sh.t2;
end

odf_norm = 0.15*max(w(:));
w = w/odf_norm;
[Nx, Ny, Nz, N_nodes] = size(w);

if nargin < 3
    flip = [0 0 0];
end

if nargin < 4
    mask = true([Nx Ny Nz]);
    h = mdm_nii_h_empty();
else
    mask = mdm_nii_read(mask_file);
    h = mdm_nii_read_header(mask_file);
end

if nargin < 5
    output_directory = pwd;
end

if nargin < 6
    method = 'dtd';
end

%%%%%%%%%

% Bounds for colormaps
disomin = 1e-10; disomax = 2.5e-9;
sddeltamin = 0.5^2; sddeltamax = 1;

if FLAG_relaxation || FLAG_double_relaxation
    if strcmp(method, 'dtr2d')
        rmin = 0; rmax = 30;
        tmin = 1/rmax; tmax = 0.11;
    elseif strcmp(method, 'dtr1d')
        rmin = 0.1; rmax = 0.9;
        tmin = 1/rmax; tmax = 1/rmin;
    elseif strcmp(method, 'dtr2r1d')
        r2min = 0; r2max = 30;
        t2min = 1/r2max; t2max = 0.11;
        
        r1min = 0.1; r1max = 0.9;
        t1min = 1/r1max; t1max = 1/r1min;
    end
end

% Normalization for colormap
diso = (diso - disomin)/(disomax - disomin);
sqddelta = (sqddelta - sddeltamin)/(sddeltamax - sddeltamin);
diso(diso < 0) = 0;
diso(diso > 1) = 1;
sqddelta(sqddelta < 0) = 0;
sqddelta(sqddelta > 1) = 1;
if FLAG_relaxation 
    r = (r - rmin)/(rmax - rmin);
    t = (t - tmin)/(tmax - tmin);
    r(r < 0) = 0;
    r(r > 1) = 1;
    t(t < 0) = 0;
    t(t > 1) = 1;
elseif FLAG_double_relaxation
    r1 = (r1 - r1min)/(r1max - r1min);
    t1 = (t1 - t1min)/(t1max - t1min);
    r1(r1 < 0) = 0;
    r1(r1 > 1) = 1;
    t1(t1 < 0) = 0;
    t1(t1 > 1) = 1;
    
    r2 = (r2 - r2min)/(r2max - r2min);
    t2 = (t2 - t2min)/(t2max - t2min);
    r2(r2 < 0) = 0;
    r2(r2 > 1) = 1;
    t2(t2 < 0) = 0;
    t2(t2 > 1) = 1;
end

% Plan for a direction flip
x = (-1)^flip(1)*x;
y = (-1)^flip(2)*y;
z = (-1)^flip(3)*z;
[theta_list, phi_list] = cartesian2spherical_unit_sphere(x,y,z);
phi_list = phi_list';

% Prepare SH matrix
[m_list, n_list] = SH_ind_list(SH_order);
N_SH_coeff = length(m_list);
SH_matrix = compute_SH_matrix(m_list, n_list, theta_list, phi_list);
pinv_SH_matrix = pinv(SH_matrix);

% Decrease dimensionality to lower memory cost
w = reshape(w, Nx*Ny*Nz, N_nodes);
diso = reshape(diso, Nx*Ny*Nz, N_nodes);
sqddelta = reshape(sqddelta, Nx*Ny*Nz, N_nodes);
if FLAG_relaxation 
    r = reshape(r, Nx*Ny*Nz, N_nodes);
    t = reshape(t, Nx*Ny*Nz, N_nodes);
elseif FLAG_double_relaxation
    r1 = reshape(r1, Nx*Ny*Nz, N_nodes);
    t1 = reshape(t1, Nx*Ny*Nz, N_nodes);
    r2 = reshape(r2, Nx*Ny*Nz, N_nodes);
    t2 = reshape(t2, Nx*Ny*Nz, N_nodes);
end
mask = reshape(mask, Nx*Ny*Nz, 1)';

% Retain only non-masked-out voxels
w = w(mask>0, :);
diso = diso(mask>0, :);
sqddelta = sqddelta(mask>0, :);
if FLAG_relaxation 
    r = r(mask>0, :);
    t = t(mask>0, :);
elseif FLAG_double_relaxation
    r1 = r1(mask>0, :);
    t1 = t1(mask>0, :);
    r2 = r2(mask>0, :);
    t2 = t2(mask>0, :);
end

[N, ~] = size(w);
SH_coeffs_map = zeros([N N_SH_coeff]);
SH_coeffs_map_normalized = zeros([N N_SH_coeff]);
SH_coeffs_diso_map = zeros([N N_SH_coeff]);
SH_coeffs_sqddelta_map = zeros([N N_SH_coeff]);
if FLAG_relaxation
    SH_coeffs_r_map = zeros([N N_SH_coeff]);
    SH_coeffs_t_map = zeros([N N_SH_coeff]);
elseif FLAG_double_relaxation
    SH_coeffs_r1_map = zeros([N N_SH_coeff]);
    SH_coeffs_t1_map = zeros([N N_SH_coeff]);
    SH_coeffs_r2_map = zeros([N N_SH_coeff]);
    SH_coeffs_t2_map = zeros([N N_SH_coeff]);
end

if FLAG_double_relaxation
    parfor v = 1:N
    % for v = 1:N
        w_vox = squeeze(w(v,:))';
        w_vox_normalized = msf_notfinite2zero(w_vox/nanmax(w_vox));
        diso_vox = squeeze(diso(v,:))';
        sqddelta_vox = squeeze(sqddelta(v,:))';
        
        SH_coeffs_map(v,:) = pinv_SH_matrix*w_vox;
        SH_coeffs_map_normalized(v,:) = pinv_SH_matrix*w_vox_normalized;
        SH_coeffs_diso_map(v,:) = pinv_SH_matrix*diso_vox;
        SH_coeffs_sqddelta_map(v,:) = pinv_SH_matrix*sqddelta_vox;
        
        r1_vox = squeeze(r1(v,:))';
        t1_vox = squeeze(t1(v,:))';
        SH_coeffs_r1_map(v,:) = pinv_SH_matrix*r1_vox;
        SH_coeffs_t1_map(v,:) = pinv_SH_matrix*t1_vox;
        
        r2_vox = squeeze(r2(v,:))';
        t2_vox = squeeze(t2(v,:))';
        SH_coeffs_r2_map(v,:) = pinv_SH_matrix*r2_vox;
        SH_coeffs_t2_map(v,:) = pinv_SH_matrix*t2_vox;
    end
elseif FLAG_relaxation
    parfor v = 1:N
        w_vox = squeeze(w(v,:))';
        w_vox_normalized = msf_notfinite2zero(w_vox/nanmax(w_vox));
        diso_vox = squeeze(diso(v,:))';
        sqddelta_vox = squeeze(sqddelta(v,:))';
        
        SH_coeffs_map(v,:) = pinv_SH_matrix*w_vox;
        SH_coeffs_map_normalized(v,:) = pinv_SH_matrix*w_vox_normalized;
        SH_coeffs_diso_map(v,:) = pinv_SH_matrix*diso_vox;
        SH_coeffs_sqddelta_map(v,:) = pinv_SH_matrix*sqddelta_vox;
        
        r_vox = squeeze(r(v,:))';
        t_vox = squeeze(t(v,:))';
        SH_coeffs_r_map(v,:) = pinv_SH_matrix*r_vox;
        SH_coeffs_t_map(v,:) = pinv_SH_matrix*t_vox;    
    end
else
    parfor v = 1:N
        w_vox = squeeze(w(v,:))';
        w_vox_normalized = msf_notfinite2zero(w_vox/nanmax(w_vox));
        diso_vox = squeeze(diso(v,:))';
        sqddelta_vox = squeeze(sqddelta(v,:))';
        
        SH_coeffs_map(v,:) = pinv_SH_matrix*w_vox;
        SH_coeffs_map_normalized(v,:) = pinv_SH_matrix*w_vox_normalized;
        SH_coeffs_diso_map(v,:) = pinv_SH_matrix*diso_vox;
        SH_coeffs_sqddelta_map(v,:) = pinv_SH_matrix*sqddelta_vox;   
    end
end

dummy_SH = zeros(Nx*Ny*Nz, N_SH_coeff);
dummy_SH(mask>0, :) = SH_coeffs_map;
SH_coeffs_map = single(reshape(dummy_SH, Nx, Ny, Nz, N_SH_coeff));

dummy_SH = zeros(Nx*Ny*Nz, N_SH_coeff);
dummy_SH(mask>0, :) = SH_coeffs_map_normalized;
SH_coeffs_map_normalized = single(reshape(dummy_SH, Nx, Ny, Nz, N_SH_coeff));

dummy_SH = zeros(Nx*Ny*Nz, N_SH_coeff);
dummy_SH(mask>0, :) = SH_coeffs_diso_map;
SH_coeffs_diso_map = single(reshape(dummy_SH, Nx, Ny, Nz, N_SH_coeff));

dummy_SH = zeros(Nx*Ny*Nz, N_SH_coeff);
dummy_SH(mask>0, :) = SH_coeffs_sqddelta_map;
SH_coeffs_sqddelta_map = single(reshape(dummy_SH, Nx, Ny, Nz, N_SH_coeff));

if FLAG_relaxation
    dummy_SH = zeros(Nx*Ny*Nz, N_SH_coeff);
    dummy_SH(mask>0, :) = SH_coeffs_r_map;
    SH_coeffs_r_map = single(reshape(dummy_SH, Nx, Ny, Nz, N_SH_coeff));
    
    dummy_SH = zeros(Nx*Ny*Nz, N_SH_coeff);
    dummy_SH(mask>0, :) = SH_coeffs_t_map;
    SH_coeffs_t_map = single(reshape(dummy_SH, Nx, Ny, Nz, N_SH_coeff));
elseif FLAG_double_relaxation
    dummy_SH = zeros(Nx*Ny*Nz, N_SH_coeff);
    dummy_SH(mask>0, :) = SH_coeffs_r1_map;
    SH_coeffs_r1_map = single(reshape(dummy_SH, Nx, Ny, Nz, N_SH_coeff));
    
    dummy_SH = zeros(Nx*Ny*Nz, N_SH_coeff);
    dummy_SH(mask>0, :) = SH_coeffs_t1_map;
    SH_coeffs_t1_map = single(reshape(dummy_SH, Nx, Ny, Nz, N_SH_coeff));
    
    dummy_SH = zeros(Nx*Ny*Nz, N_SH_coeff);
    dummy_SH(mask>0, :) = SH_coeffs_r2_map;
    SH_coeffs_r2_map = single(reshape(dummy_SH, Nx, Ny, Nz, N_SH_coeff));
    
    dummy_SH = zeros(Nx*Ny*Nz, N_SH_coeff);
    dummy_SH(mask>0, :) = SH_coeffs_t2_map;
    SH_coeffs_t2_map = single(reshape(dummy_SH, Nx, Ny, Nz, N_SH_coeff));
end

SH_coeffs_map = msf_notfinite2zero(SH_coeffs_map);
SH_coeffs_map_normalized = msf_notfinite2zero(SH_coeffs_map_normalized);
SH_coeffs_diso_map = msf_notfinite2zero(SH_coeffs_diso_map);
SH_coeffs_sqddelta_map = msf_notfinite2zero(SH_coeffs_sqddelta_map);

mdm_nii_write(single(SH_coeffs_map), fullfile(output_directory, 'SH_coeffs_ODF.nii.gz'), h);
mdm_nii_write(single(SH_coeffs_map_normalized), fullfile(output_directory, 'SH_coeffs_ODF_normalized.nii.gz'), h);
mdm_nii_write(single(SH_coeffs_diso_map), fullfile(output_directory, 'SH_coeffs_ODF_diso.nii.gz'), h);
mdm_nii_write(single(SH_coeffs_sqddelta_map), fullfile(output_directory, 'SH_coeffs_ODF_sqddelta.nii.gz'), h);

if FLAG_relaxation
    SH_coeffs_r_map = msf_notfinite2zero(SH_coeffs_r_map);
    SH_coeffs_t_map = msf_notfinite2zero(SH_coeffs_t_map);
    if strcmp(method, 'dtr2d')
        mdm_nii_write(single(SH_coeffs_r_map), fullfile(output_directory, 'SH_coeffs_ODF_r2.nii.gz'), h);
        mdm_nii_write(single(SH_coeffs_t_map), fullfile(output_directory, 'SH_coeffs_ODF_t2.nii.gz'), h);
    elseif strcmp(method, 'dtr1d')
        mdm_nii_write(single(SH_coeffs_r_map), fullfile(output_directory, 'SH_coeffs_ODF_r1.nii.gz'), h);
        mdm_nii_write(single(SH_coeffs_t_map), fullfile(output_directory, 'SH_coeffs_ODF_t1.nii.gz'), h);
    end
elseif FLAG_double_relaxation
    mdm_nii_write(single(SH_coeffs_r1_map), fullfile(output_directory, 'SH_coeffs_ODF_r1.nii.gz'), h);
    mdm_nii_write(single(SH_coeffs_t1_map), fullfile(output_directory, 'SH_coeffs_ODF_t1.nii.gz'), h);
    
    mdm_nii_write(single(SH_coeffs_r2_map), fullfile(output_directory, 'SH_coeffs_ODF_r2.nii.gz'), h);
    mdm_nii_write(single(SH_coeffs_t2_map), fullfile(output_directory, 'SH_coeffs_ODF_t2.nii.gz'), h);
end
