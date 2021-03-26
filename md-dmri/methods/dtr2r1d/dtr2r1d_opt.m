function opt = dtr2r1d_opt(opt)
% function opt = dtr2r1d_opt(opt)
%
% Makes sure that all needed fields in the options structure are present

opt.dtr2r1d.present = 1;

opt.dtr2r1d = msf_ensure_field(opt.dtr2r1d, 'tmp', 1); 
opt.dtr2r1d = msf_ensure_field(opt.dtr2r1d, 'do_plot', 0);
opt.dtr2r1d = msf_ensure_field(opt.dtr2r1d, 'ind_start', 1);
opt.dtr2r1d = msf_ensure_field(opt.dtr2r1d, 'dmin', 5e-11);
opt.dtr2r1d = msf_ensure_field(opt.dtr2r1d, 'dmax', 5e-9);
opt.dtr2r1d = msf_ensure_field(opt.dtr2r1d, 'r2min', 1);
opt.dtr2r1d = msf_ensure_field(opt.dtr2r1d, 'r2max', 30);
opt.dtr2r1d = msf_ensure_field(opt.dtr2r1d, 'r1min', 0.2);
opt.dtr2r1d = msf_ensure_field(opt.dtr2r1d, 'r1max', 2);
opt.dtr2r1d = msf_ensure_field(opt.dtr2r1d, 'n_in', 2e2); % number of nodes in NNLS inversion. [100 - 1000]
opt.dtr2r1d = msf_ensure_field(opt.dtr2r1d, 'n_out', 50); % number of nodes constituting a solution
opt.dtr2r1d = msf_ensure_field(opt.dtr2r1d, 'n_kill', 0); 
opt.dtr2r1d = msf_ensure_field(opt.dtr2r1d, 'n_proliferation', 20); % number of initial iterations
opt.dtr2r1d = msf_ensure_field(opt.dtr2r1d, 'n_extinction', 20); % number of mutation rounds
opt.dtr2r1d = msf_ensure_field(opt.dtr2r1d, 'ofuzz', .1*2*pi);
opt.dtr2r1d = msf_ensure_field(opt.dtr2r1d, 'dfuzz', .1);
opt.dtr2r1d = msf_ensure_field(opt.dtr2r1d, 'r2fuzz', .1);
opt.dtr2r1d = msf_ensure_field(opt.dtr2r1d, 'r1fuzz', .1);

opt.dtr2r1d = msf_ensure_field(opt.dtr2r1d, 't1_weighting', 'spin_echo'); % 'spin_echo', 'saturation', ('wip_inversion' only for 1d_fit2data)

opt.dtr2r1d = msf_ensure_field(opt.dtr2r1d, 'do_bin', 0);
opt.dtr2r1d = msf_ensure_field(opt.dtr2r1d, 'bin_r1min', [0 0 0]);
opt.dtr2r1d = msf_ensure_field(opt.dtr2r1d, 'bin_r1max', [20 20 20]);
opt.dtr2r1d = msf_ensure_field(opt.dtr2r1d, 'bin_r2min', [1 1 1]);
opt.dtr2r1d = msf_ensure_field(opt.dtr2r1d, 'bin_r2max', [30 30 30]);
opt.dtr2r1d = msf_ensure_field(opt.dtr2r1d, 'bin_disomin', [2 .1 .1]*1e-9);
opt.dtr2r1d = msf_ensure_field(opt.dtr2r1d, 'bin_disomax', [3.5 2 2]*1e-9);
opt.dtr2r1d = msf_ensure_field(opt.dtr2r1d, 'bin_dratiomin', [0 0 4]);
opt.dtr2r1d = msf_ensure_field(opt.dtr2r1d, 'bin_dratiomax', [4 4 200]);


% 's_fw','s_st','s_lt','s_pt','s_plt', 's_splt','s_residue'
opt.dtr2r1d = msf_ensure_field(opt.dtr2r1d, 'fig_maps', {'s0', 'mdiso','vdiso_n','mdaniso_n','vdaniso_n','msdaniso_n','vsdaniso_n','mr2','vr2','mr1','vr1','ufa','fa','op'});
opt.dtr2r1d = msf_ensure_field(opt.dtr2r1d, 'fig_prefix', 'dtr2r1d');
opt.dtr2r1d = msf_ensure_field(opt.dtr2r1d, 'fig_cmaps',{'fa','cl','cp','ufa','ucl','ucp'});
opt.dtr2r1d = msf_ensure_field(opt.dtr2r1d, 'fig_ccol',{'t1x6','lambda33vec','lambda11vec','s1x6prim','s1x6prim','s1x6prim'});
opt.dtr2r1d = msf_ensure_field(opt.dtr2r1d, 'fig_ccolnorm',{'lambda33','mask','mask','slambda33prim','slambda33prim','slambda33prim'});
opt.dtr2r1d = msf_ensure_field(opt.dtr2r1d, 'do_dtr2r1dpdf', 0);

