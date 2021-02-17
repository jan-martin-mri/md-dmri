function opt = dtd_mv_gamma_opt(opt)
% function opt = dtd_gamma_opt(opt)
%
% Makes sure that all needed fields in the options structure are present

opt.dtd_mv_gamma.present = 1;

opt.dtd_mv_gamma = msf_ensure_field(opt.dtd_mv_gamma, 'tmp', 1); 
opt.dtd_mv_gamma = msf_ensure_field(opt.dtd_mv_gamma, 'lsq_opts', optimoptions('lsqcurvefit', 'display', 'off','MaxFunEvals',1e3));
opt.dtd_mv_gamma = msf_ensure_field(opt.dtd_mv_gamma, 'do_plot', 0);
opt.dtd_mv_gamma = msf_ensure_field(opt.dtd_mv_gamma, 'do_weight', 1);
opt.dtd_mv_gamma = msf_ensure_field(opt.dtd_mv_gamma, 'weight_sthresh', .07); % 0.07
opt.dtd_mv_gamma = msf_ensure_field(opt.dtd_mv_gamma, 'weight_wthresh', 2); % 2
opt.dtd_mv_gamma = msf_ensure_field(opt.dtd_mv_gamma, 'weight_mdthresh', 1e-9);
opt.dtd_mv_gamma = msf_ensure_field(opt.dtd_mv_gamma, 'fig_maps', {'s0', 'mdiso', 'FA', 'vdiso', 'msdaniso', 'vdison', 'msdanison', 'sOP'});
opt.dtd_mv_gamma = msf_ensure_field(opt.dtd_mv_gamma, 'fig_prefix', 'dtd_mv_gamma');
opt.dtd_mv_gamma = msf_ensure_field(opt.dtd_mv_gamma, 'fig_cmaps',    {'FA', 'msdanison'});
opt.dtd_mv_gamma = msf_ensure_field(opt.dtd_mv_gamma, 'fig_ccol',     {'t1x6','t1x6'});
opt.dtd_mv_gamma = msf_ensure_field(opt.dtd_mv_gamma, 'fig_ccolnorm', {'maxmdii','maxmdii'});
%     {'s0','MD', 'MKi', 'MKa', 'MKt', 'ufa', 'miso', 'viso_n', 'msaniso_n'});

opt.dtd_mv_gamma = msf_ensure_field(opt.dtd_mv_gamma, 'diso_min', 0);
opt.dtd_mv_gamma = msf_ensure_field(opt.dtd_mv_gamma, 'diso_max', 3.5*1e-9);
opt.dtd_mv_gamma = msf_ensure_field(opt.dtd_mv_gamma, 'vdiso_min', 0);
opt.dtd_mv_gamma = msf_ensure_field(opt.dtd_mv_gamma, 'vdiso_max', 8e-18); % 2e-18 for in vivo
opt.dtd_mv_gamma = msf_ensure_field(opt.dtd_mv_gamma, 'vdison_min', 0);
opt.dtd_mv_gamma = msf_ensure_field(opt.dtd_mv_gamma, 'vdison_max', 1);
opt.dtd_mv_gamma = msf_ensure_field(opt.dtd_mv_gamma, 'msdaniso_min', 0);
opt.dtd_mv_gamma = msf_ensure_field(opt.dtd_mv_gamma, 'msdaniso_max', 1e-18);
opt.dtd_mv_gamma = msf_ensure_field(opt.dtd_mv_gamma, 'msdanison_min', 0);
opt.dtd_mv_gamma = msf_ensure_field(opt.dtd_mv_gamma, 'msdanison_max', 0.7);
opt.dtd_mv_gamma = msf_ensure_field(opt.dtd_mv_gamma, 'FA_min', 0);
opt.dtd_mv_gamma = msf_ensure_field(opt.dtd_mv_gamma, 'FA_max', 0.7);
opt.dtd_mv_gamma = msf_ensure_field(opt.dtd_mv_gamma, 'sOP_min', 0.6);
opt.dtd_mv_gamma = msf_ensure_field(opt.dtd_mv_gamma, 'sOP_max', 1);

% Number of identical fitting repetitions (keeping smallest residual)
% Very expensive
opt.dtd_mv_gamma = msf_ensure_field(opt.dtd_mv_gamma, 'fit_iters', 1);

% Number of random guesses to start from (keeping guess with smallest residual)
% Not expensive
opt.dtd_mv_gamma = msf_ensure_field(opt.dtd_mv_gamma, 'guess_iters', 200);

% Decide if multiple series should be assumed to have different baseline
% signal (s_ind has multiple unique values).
opt.dtd_mv_gamma = msf_ensure_field(opt.dtd_mv_gamma, 'do_multiple_s0', 0);

% Bounds for initial guess and fitting (not including relative signal for
% multiple series. This is done in the data2fit function.
% Note that the first value is S(b=0)/max(signal), i.e., a relative signal.
% <D> bounded between 1e-12 and 4e-9 and H = 1/(kappa*I + Theta) bounded between 1e-12 and 1e3 (Psi = <D>H)
lb_kappa = 1+1e-12;
ub_kappa = 10;
lb_averageD = 1e-12;
ub_averageD = 3.5e-9;
lb_H = 1e-12;
ub_H = 1e2;
lb_psi = lb_averageD*lb_H;
ub_psi = ub_averageD*ub_H;
opt.dtd_mv_gamma = msf_ensure_field(opt.dtd_mv_gamma, 'fit_lb', [ 0 lb_kappa [1 1 1]*lb_psi [1 1 1]*lb_H [0 0 0]]);
opt.dtd_mv_gamma = msf_ensure_field(opt.dtd_mv_gamma, 'fit_ub', [10 ub_kappa [1 1 1]*ub_psi [1 1 1]*ub_H [2 1 2]*pi]);
