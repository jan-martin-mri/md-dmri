function list = dtd_mv_gamma_clamp(list, bounds, new_values)
% function opt = dtd_gamma_clamp(opt, bounds, new_values)

if nargin < 3
    new_values = bounds;
end

list(list > max(bounds)) = max(new_values);
list(list < min(bounds)) = min(new_values);