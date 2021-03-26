function fn = dtr2r1d_pipe(s, paths, opt)
% function fn = dtr2r1d_pipe(s, paths, opt)
%
% s     - input structure
% paths - either a pathname or a path structure (see mdm_paths)
% opt   - (optional) options that drive the pipeline
%
% fn    - a cell array with filenames to generated nii files

fn = '';

if (nargin < 3), opt.present = 1; end

opt = mdm_opt(opt);
opt = dtr2r1d_opt(opt);
paths = mdm_paths(paths);
msf_log(['Starting ' mfilename], opt);

% Prepare: mask etc
if (opt.do_mask)
    s = mdm_s_mask(s, @mio_mask_threshold, [], opt);
end

% jm, 2021-03-26: Note that bootstrapping with replacement might no be
% ideal here. We could instead use Monte Carlo cross-validation by randomly
% sampling a fixed percentage of the acquisition scheme and not replacing
% any data points. Code:
% Siehe randperm(n, k), berechne k aus n mit neuem opt Parameter der den
% Anteil an unique Messungen pro bootstrap festsetzt.
% Run the analysis
if (opt.do_data2fit)
    if opt.do_bootstrap
        % Create bootstrap directory
        msf_mkdir(fileparts(paths.ind_fn));
        % Randomly choose signal indices for bootstrapping
        ind = (opt.dtr2r1d.ind_start-0) + round(rand([s.xps.n-(opt.dtr2r1d.ind_start-1),1])*(s.xps.n-(opt.dtr2r1d.ind_start-0)));
        % Save bootstrap signal indices
        save(paths.ind_fn, 'ind');
        ind_fn = mdm_ind_save(ind, paths.ind_fn);
        load(paths.ind_fn)
        % Set bootstrap index in opt struct
        opt.bootstrap.ind = ind;
    end
    % Fit data
    mdm_data2fit(@dtr2r1d_4d_data2fit, s, paths.mfs_fn, opt);
    % Convert mfs.m to single to save disk space
    mfs = mdm_mfs_load(paths.mfs_fn);
    mfs.m = single(mfs.m);
    save(paths.mfs_fn, 'mfs');
end

if (opt.do_fit2param)
    mdm_fit2param(@dtr2r1d_4d_fit2param, paths.mfs_fn, paths.dps_fn, opt);
end

% Save nifti parameter maps    
if (opt.do_param2nii)
    fn = mdm_param2nii(paths.dps_fn, paths.nii_path, opt.dtr2r1d, opt);
    % Convert all .nii.gz files to .pdf
    if (opt.do_nii2pdf)
        mdm_nii2pdf(fn, [], opt);
    end
end

% Save dtr2r1d pdf   
if (opt.do_m2pdf)
    dtr2r1d_m2pdf(paths.dps_fn, paths.nii_path, opt);
end



