function SH_matrix = compute_SH_matrix(m_list, n_list,theta_list,phi_list)
% Tournier convention

N_SH_coeff = length(m_list);
N_odf = length(theta_list);
SH_matrix = zeros([N_odf, N_SH_coeff]);

for i = 1:N_SH_coeff
    n = n_list(i);
    m = m_list(i);
    if m > 0
        SH_matrix(:,i) = (-1)^m*compute_SH(n, m, theta_list, phi_list, 'type', 'real')';
    elseif m < 0
        SH_matrix(:,i) = (-1)^m*compute_SH(n, abs(m), theta_list, phi_list, 'type', 'imag')';
    else
        SH_matrix(:,i) = compute_SH(n, 0, theta_list, phi_list)';
    end
end