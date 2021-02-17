function [m_list, n_list] = SH_ind_list(SH_order)
% n is equivalent to l in quantum physics
% m remains unchanged from quantum physics

n_range = linspace(0, SH_order, SH_order/2+1);
n_list = repelem(n_range, 2*n_range+1); % For each n, 2n+1 values of m
nb_coef = (SH_order + 2)*(SH_order + 1)/ 2; % Number of integer element between 0 and SH_order (inclusive)
offset = 0;
m_list = zeros([1, nb_coef]);
for i = n_range
    m_list(offset+1:(offset + 2*i + 1)) = -i:i;
    offset = offset + 2*i + 1;
end