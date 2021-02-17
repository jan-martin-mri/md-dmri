function exitflag = run_inversion_pipeline(fn_data, fn_mask, fn_xps)
%% run_inversion_pipeline
% Settings for md-dmri data reconstruction using nonparametric Monte Carlo inversion.
% This script currently supports the following methods:
%     dtd
%     dtr1d
%     dtr2d
%     dtr2r1d
% Input
%     fn_data: Path to the data in .nii.gz format
%     fn_mask: Path to data mask in .nii.gz format. Replace with '' if no
%              mask is available.
%     fn_xps: Path to xps that corresponds to data in .mat format

exitflag = 0;
tStart = tic;

%% Input files
% Transform input to match Alexis way of handling files.
% Long-term this should be done with absolute paths and no fullfile during
% individual steps.
[pn_d, fn_d, ext_d] = msf_fileparts(fn_data);
input_parameters.initial_directory = pn_d; % Where the data can be found
input_parameters.data_file = [fn_d ext_d]; % Diffusion/relaxation data
[pn_x, fn_x, ext_x] = msf_fileparts(fn_xps);
input_parameters.xps_file = [fn_x ext_x]; % xps matfile
cond1 = ~strcmp(pn_d, pn_x);
cond2 = 0;
if ~isempty(fn_mask)
    [pn_m, fn_m, ext_m] = msf_fileparts(fn_mask);
    input_parameters.mask_file = [fn_m ext_m]; % Data mask - If no mask, simply write '', one will be generated from mdm_s_mask
    cond2 = ~strcmp(pn_d, pn_m);
else
    input_parameters.mask_file = '';
end
% Check if data, mask, and xps are located in the same folder
if cond1 || cond2
    error('run_inversion_pipeline: Data, mask, and xps need to be at the same location.');
end

%% Inversion method
input_parameters.inversion_method = 'dtr2r1d'; % 'dtr2r1d', 'dtr2d', 'dtr1d', 'dtd'
input_parameters.do_covariance = 0;

%% Data correction (extrapolated references + Elastix)
input_parameters.do_data_correction = 0;
input_parameters.do_rotate_bvec = 0;

%% Global parameters
input_parameters.nb_MC_inversions = 100; % Number of Monte Carlo realizations
input_parameters.nb_parallel_workers = inf; % Number of parallel workers, put Inf to have the largest number of workers on your machine
input_parameters.threshold_d_iso_big = 2.5; % Threshold on d_iso to define the big bin (in {\mu}m^2/ms)
input_parameters.threshold_d_delta_thin = 0.5; % Threshold on anisotropy to define the thin bin
input_parameters.threshold_weight_thin = 0.05; % Threshold on the thin-bin total weight to create a thin-bin data mask (useless for now)
input_parameters.dir_flips = [0 0 0]; % Orientation filps for ODFs and clustering in case of gradient flips in the data [x y z] = [0/1 0/1 0/1]

%% Colormap bounds
input_parameters.clim.diso_clim = 3.5e-9*[0 1];
input_parameters.clim.msqddelta_clim = 1*[0 1];
input_parameters.clim.vdiso_clim = .3*3e-9^2*[0 1];
input_parameters.clim.vsqddelta_clim = .1*[0 1];
input_parameters.clim.cvdisosqddelta_clim = .1*3e-9*1*[-1 1];

if strcmp(input_parameters.inversion_method, 'dtr2d')
    mr_clim = 30*[0 1];
    input_parameters.clim.mr_clim = mr_clim;
    input_parameters.clim.vr_clim = .2*max(mr_clim)^2*[0 1];
    input_parameters.clim.cvdisor_clim = 0.1*max(mr_clim)*3e-9*[-1 1];
    input_parameters.clim.cvsqddeltar_clim = 0.1*max(mr_clim)*[-1 1];
    
    input_parameters.clim.mt_clim = [1/max(mr_clim) 0.11];
    input_parameters.clim.mt_clim_blackandwhite = [1/max(mr_clim) 0.15];
    input_parameters.clim.vt_clim = 0.5*[0 1];
    input_parameters.clim.cvdisot_clim = 1e-9*[-1 1];
    input_parameters.clim.cvsqddeltat_clim = 0.05*[-1 1];  
    
elseif strcmp(input_parameters.inversion_method, 'dtr1d')
    mr_clim = 0.8*[0 1];
    input_parameters.clim.mr_clim = mr_clim;  
    input_parameters.clim.vr_clim = .2*max(mr_clim)^2*[0 1];
    input_parameters.clim.cvdisor_clim = 0.1*max(mr_clim)*3e-9*[-1 1];
    input_parameters.clim.cvsqddeltar_clim = 0.1*max(mr_clim)*1*[-1 1];
    
    mt_clim = [0.5 5];
    input_parameters.clim.mt_clim = mt_clim;
    input_parameters.clim.mt_clim_blackandwhite = [1 10];
    input_parameters.clim.vt_clim = 0.25*[0 max(mt_clim)^2];
    input_parameters.clim.cvdisot_clim = 0.1*max(mt_clim)*3e-9*[-1 1];
    input_parameters.clim.cvsqddeltat_clim = 0.1*max(mt_clim)*1*[-1 1]; 
    
elseif strcmp(input_parameters.inversion_method, 'dtr2r1d')
    mr2_clim = 30*[0 1];
    input_parameters.clim.mr2_clim = mr2_clim;
    input_parameters.clim.vr2_clim = .2*max(mr2_clim)^2*[0 1];
    input_parameters.clim.cvdisor2_clim = 0.1*max(mr2_clim)*3e-9*[-1 1];
    input_parameters.clim.cvsqddeltar2_clim = 0.1*max(mr2_clim)*[-1 1];
    
    mt2_clim = [1/max(mr2_clim) 0.11];
    input_parameters.clim.mt2_clim = mt2_clim;
    input_parameters.clim.mt2_clim_blackandwhite = [1/max(mr2_clim) 0.15];
    input_parameters.clim.vt2_clim = 0.5*[0 1];
    input_parameters.clim.cvdisot2_clim = 1e-9*[-1 1];
    input_parameters.clim.cvsqddeltat2_clim = 0.05*[-1 1];  

    mr1_clim = 0.8*[0 1];
    input_parameters.clim.mr1_clim = mr1_clim;  
    input_parameters.clim.vr1_clim = .2*max(mr1_clim)^2*[0 1];
    input_parameters.clim.cvdisor1_clim = 0.1*max(mr1_clim)*3e-9*[-1 1];
    input_parameters.clim.cvsqddeltar1_clim = 0.1*max(mr1_clim)*1*[-1 1];
    
    mt1_clim = [0.5 5];
    input_parameters.clim.mt1_clim = mt1_clim;
    input_parameters.clim.mt1_clim_blackandwhite = [1 10];
    input_parameters.clim.vt1_clim = 0.25*[0 max(mt1_clim)^2];
    input_parameters.clim.cvdisot1_clim = 0.1*max(mt1_clim)*3e-9*[-1 1];
    input_parameters.clim.cvsqddeltat1_clim = 0.1*max(mt1_clim)*1*[-1 1]; 
    
    input_parameters.clim.cvr1r2_clim = 0.1*max(mr1_clim)*max(mr2_clim)*[-1 1];
    input_parameters.clim.cvt1t2_clim = 0.1*max(mt1_clim)*max(mt2_clim)*[-1 1];
end

% If one wants to optimize the colormap bounds without deleting input_parameters.parameter_maps_directory everytime
input_parameters.repeat_colormaps = 0; 

%% ODF parameters
input_parameters.do_maps = 0;
input_parameters.framework_directory = 'inversion_pipeline/md-dmri'; % Directory wherein 'methods' and 'tools' are located
input_parameters.mrtrix_directory = ''; % Directory wherein the binaries of mrtrix are located
input_parameters.max_nb_odf_peaks = 4; % Maximal number of ODF peaks per voxel
input_parameters.threshold_ODF_w = 0.1; % 0.05, 0.1, 0.15, take out the core of the ODF for peak calculation
input_parameters.nb_mesh_nodes = 1000; % 250, 350, 500, 1000, 3994, or 15970, to serve as projecting mesh for the discrete ODFs

%% Monte-Carlo density-peak clustering (MC-DPC)
input_parameters.do_clustering = 0; % Do you want clustering to be performed?
struct_clustering.max_nb_clusters = input_parameters.max_nb_odf_peaks; % Maximal number of clusters
struct_clustering.do_plot = 0; % Print clustering plots, do not use for big volumes
struct_clustering.min_weight_per_cluster = input_parameters.threshold_ODF_w; % Repeat clustering at the bootstrap level if one cluster does not weigh at least this
input_parameters.structure_clustering = struct_clustering;

%% Output directories
if strcmp(input_parameters.inversion_method, 'dtd') && input_parameters.do_covariance
    input_parameters.covariance_directory = fullfile(input_parameters.initial_directory, '0_covariance');
end
input_parameters.data_directory = fullfile(input_parameters.initial_directory, '0_data');
input_parameters.bootstrap_directory = fullfile(input_parameters.initial_directory, '1_bootstrap');
input_parameters.parameter_maps_directory = fullfile(input_parameters.initial_directory, '2_nii_pdf_maps');
input_parameters.odf_directory = fullfile(input_parameters.initial_directory, '3_odfs');
input_parameters.clustering_directory = fullfile(input_parameters.initial_directory, '4_clustering');

%% Inversion pipeline
global_inversion_pipeline(input_parameters)

exitflag = 1;
pipeline_runtime = toc(tStart);
disp(['run_inversion_pipeline finished in ' num2str(pipeline_runtime) ' s with exitflag ' num2str(exitflag)]);