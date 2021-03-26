function res = dtr2r1d_4d_data2fit(s, mfs_fn, opt)
% function mfs_fn = dtr2d_4d_data2fit(s, o, opt)

if (nargin < 3), opt = []; end

res = -1;

% jm, 2021-03-25: Reintroduced bootstrapping and indexing
% Note that b-values are checked and removed during 
% global_inversion_pipeline before transferring the xps to the data 
% directory.
if opt.do_bootstrap
    disp('DEBUG: opt.do_bootstrap = 1 / Using bootstrap indices.');
    ind = opt.bootstrap.ind;
else
    disp('DEBUG: opt.do_bootstrap = 0 / Using all measurements.');
    ind = opt.dtr2r1d.ind_start : s.xps.n;
end

%Verify the xps
%dti_euler_mic_check_xps(s.xps);

% Loop over the volume and fit the model
xps = s.xps; % this appears to improve parallel performance
f = @(signal) dtr2r1d_1d_data2fit(signal, xps, opt, ind);
dummy = mio_fit_model(f, s, mfs_fn, opt);

res = 1;
