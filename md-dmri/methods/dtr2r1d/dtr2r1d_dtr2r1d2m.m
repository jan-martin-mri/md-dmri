function m = dtr2r1d_dtr2r1d2m(dtr2r1d,opt)

m = zeros(1 + 7*opt.dtr2r1d.n_out,1);
if ~isempty(dtr2r1d)
    m(1:numel(dtr2r1d)) = dtr2r1d;
%     m = m(1:(1 + 7*opt.dtr2r1d.n_out),1);
    m(1) = min([dtr2r1d(1) opt.dtr2r1d.n_out]);
end
m = m';
