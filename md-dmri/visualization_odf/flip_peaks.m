directory = fullfile(pwd, 'Glioma_Spectrum', '20200504083444', '3_odfs');
data_file = fullfile(directory, 'SH_peaks.nii.gz');
flip = [1 0 0];

[I, h] = mdm_nii_read(data_file);
[~, ~, ~, N] = size(I);
nb_peaks = N/3;

ind_x = 1:nb_peaks-1:N;
ind_y = 2:nb_peaks-1:N;
ind_z = 3:nb_peaks-1:N;

I(:,:,:,ind_x) = (-1)^flip(1)*I(:,:,:,ind_x);
I(:,:,:,ind_y) = (-1)^flip(2)*I(:,:,:,ind_y);
I(:,:,:,ind_z) = (-1)^flip(3)*I(:,:,:,ind_z);

mdm_nii_write(I, data_file, h);