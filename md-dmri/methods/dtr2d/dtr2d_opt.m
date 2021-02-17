function opt = dtr2d_opt(opt)
% function opt = dtr2d_opt(opt)
%
% Makes sure that all needed fields in the options structure are present

opt.dtr2d.present = 1;

opt.dtr2d = msf_ensure_field(opt.dtr2d, 'tmp', 1); 
opt.dtr2d = msf_ensure_field(opt.dtr2d, 'do_plot', 0);
opt.dtr2d = msf_ensure_field(opt.dtr2d, 'ind_start', 1);
opt.dtr2d = msf_ensure_field(opt.dtr2d, 'dmin', 5e-11);
opt.dtr2d = msf_ensure_field(opt.dtr2d, 'dmax', 5e-9);
opt.dtr2d = msf_ensure_field(opt.dtr2d, 'r2min', 1); % 1
opt.dtr2d = msf_ensure_field(opt.dtr2d, 'r2max', 30); % 30
opt.dtr2d = msf_ensure_field(opt.dtr2d, 'n_in', 2e2); % n_in: Number of nodes in NNLS inversion. [100 - 1000]
opt.dtr2d = msf_ensure_field(opt.dtr2d, 'n_out', 50); % 20
opt.dtr2d = msf_ensure_field(opt.dtr2d, 'n_kill', 0);
opt.dtr2d = msf_ensure_field(opt.dtr2d, 'n_proliferation', 30); % 20
opt.dtr2d = msf_ensure_field(opt.dtr2d, 'n_extinction', 30); % 20
opt.dtr2d = msf_ensure_field(opt.dtr2d, 'ofuzz', .1*2*pi);
opt.dtr2d = msf_ensure_field(opt.dtr2d, 'dfuzz', .1);
opt.dtr2d = msf_ensure_field(opt.dtr2d, 'r2fuzz', .1);
 
opt.dtr2d = msf_ensure_field(opt.dtr2d, 'odf_nnodes', 1000); %250, 500, 1000, 3994, or 15970

% Statistical measures
opt.dtr2d = msf_ensure_field(opt.dtr2d, 'fig_maps', ...
   {'s0','mdiso','msqddelta','mr2','vdiso','vsqddelta','vr2','cvdisosqddelta','cvdisor2','cvsqddeltar2'});
opt.dtr2d = msf_ensure_field(opt.dtr2d, 'fig_prefix', 'dtr2d');
opt.dtr2d = msf_ensure_field(opt.dtr2d, 'fig_cmaps',{});
opt.dtr2d = msf_ensure_field(opt.dtr2d, 'fig_ccol',{});
opt.dtr2d = msf_ensure_field(opt.dtr2d, 'fig_ccolnorm',{});
opt.dtr2d = msf_ensure_field(opt.dtr2d, 'do_dtr2dpdf', 0);


% Binning options
opt.dtr2d = msf_ensure_field(opt.dtr2d, 'do_bin', 0);
opt.dtr2d = msf_ensure_field(opt.dtr2d, 'bin_disomax', [opt.dtr2d.dmax 2*1e-9 2*1e-9]);
opt.dtr2d = msf_ensure_field(opt.dtr2d, 'bin_disomin', [2*1e-9 opt.dtr2d.dmin  opt.dtr2d.dmin]);
opt.dtr2d = msf_ensure_field(opt.dtr2d, 'bin_dratiomax', [4 4 opt.dtr2d.dmax/opt.dtr2d.dmin]);
opt.dtr2d = msf_ensure_field(opt.dtr2d, 'bin_dratiomin', [opt.dtr2d.dmin/opt.dtr2d.dmax opt.dtr2d.dmin/opt.dtr2d.dmax 4]);
opt.dtr2d = msf_ensure_field(opt.dtr2d, 'bin_r2max', [opt.dtr2d.r2max opt.dtr2d.r2max opt.dtr2d.r2max]);
opt.dtr2d = msf_ensure_field(opt.dtr2d, 'bin_r2min', [opt.dtr2d.r2min opt.dtr2d.r2min opt.dtr2d.r2min]);




