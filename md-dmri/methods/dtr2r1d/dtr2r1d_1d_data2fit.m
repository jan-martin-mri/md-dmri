function m = dtr2r1d_1d_data2fit(signal, xps, opt, ind)
% function m = dtr2r1d_1d_data2fit(signal, xps, opt, ind)
%
% Size-shape-orientation diffusion tensor distribution

if (nargin < 4), ind = ones(size(signal)) > 0; end

% Filter out signal points with 0 or negative value
% Implemented by Jan Martin on 2020-12-09
ind(signal(ind) <= 0) = [];

% bt_mx6 = xps.bt(find(ind),:);
bt_mx6 = xps.bt(ind,:); % b-tensor (Mandel notation) 
te = xps.te(ind); % echo time
tr = xps.tr(ind); % saturation/recovery time
% tr(1) = 0;
stemp = signal(ind);

% Generate n_proliferation sets of n_nodes nodes &
% sequentially fit them to the data
dtr2r1d = dtr2r1d_proliferation(stemp, bt_mx6, te, tr, opt);

%dtr2d
%pause

dtr2r1d = dtr2r1d_extinction(stemp, bt_mx6, te, tr, dtr2r1d, opt);

m = dtr2r1d_dtr2r1d2m(dtr2r1d,opt);
% size(m)
if (opt.dtr2r1d.do_plot)
    figure(1), clf
    signal_fit = dtr2r1d_1d_fit2data(m, xps);
    %[~,s_ind] = sort(signal_fit,'descend');
    %semilogy(xps.b,signal,'.',xps.b,signal_fit,'o',xps.b,m(1)*weight,'x');
    plot(1:xps.n,signal,'o',1:xps.n,signal_fit,'x');
    %set(gca,'YLim',m(1)*[.01 1.2])
    pause(0.05);
end