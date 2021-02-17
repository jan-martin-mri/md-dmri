function m = dtr2d_1d_data2fit(signal, xps, opt, ind)
% function m = dtr2d_1d_data2fit(signal, xps, opt, ind)
%
% 5D diffusion tensor & R2 distribution

if (nargin < 4), ind = ones(size(signal)) > 0; end

bt_mx6 = xps.bt(ind,:); % b-tensor (Mandel notation) 
te = xps.te(ind); % echo time
stemp = signal(ind);

% Generate opt.n_proliferation sets of opt.n_in nodes &
% sequentially fit them to the data
dtr2d = dtr2d_proliferation(stemp, bt_mx6, te, opt);
%dtr2d
%pause


dtr2d = dtr2d_extinction(stemp, bt_mx6, te, dtr2d, opt);
m = dtr2d_dtr2d2m(dtr2d,opt);


% size(m)
if (opt.dtr2d.do_plot)
    figure(1), clf
    signal_fit = dtr2d_1d_fit2data(m, xps);
    %[~,s_ind] = sort(signal_fit,'descend');
    %semilogy(xps.b,signal,'.',xps.b,signal_fit,'o',xps.b,m(1)*weight,'x');
    plot(1:xps.n,signal,'o',1:xps.n,signal_fit,'x');
    %set(gca,'YLim',m(1)*[.01 1.2])
    pause(0.05);
end