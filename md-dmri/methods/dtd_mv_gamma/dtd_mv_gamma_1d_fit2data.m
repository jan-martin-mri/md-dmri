function s = dtd_mv_gamma_1d_fit2data(m, xps)
% function s = dtd_mv_gamma_1d_fit2data(m, xps)

% Convert to readable parameters
s0    = m(1);         % Signal at no diff encoding (baseline signal)
kappa = m(2);
psi_1 = m(3);
psi_2 = m(4);
psi_3 = m(5);
h_1   = m(6);
h_2   = m(7);
h_3   = m(8);
alpha = m(9);
beta  = m(10);
gamma = m(11);

% Define Psi and H, with their rotations
Psi_diag = diag([psi_1, psi_2, psi_3]);
% H_diag = diag([h_1, h_2, h_3]);
H_inv_diag = diag([1/h_1, 1/h_2, 1/h_3]);

Rz_alpha = [[cos(alpha) -sin(alpha) 0]; ...
           [sin(alpha) cos(alpha)  0]; ...
           [0          0           1]];
Ry_beta =  [[cos(beta)  0 sin(beta)]; ...
           [0          1         0]; ...
           [-sin(beta) 0 cos(beta)]];
Rz_gamma = [[cos(gamma) -sin(gamma) 0]; ...
           [sin(gamma) cos(gamma)  0]; ...
           [0          0           1]];

R_euler = Rz_gamma*Ry_beta*Rz_alpha;
Psi = R_euler*Psi_diag/R_euler;
% H = R_euler*H_diag/R_euler; 
H_inv = R_euler*H_inv_diag/R_euler; 

% Signal equation
s = zeros([xps.n,1]);
for i = 1:xps.n
    bt = xps.bt(i,:);
    product_Psi_B = Psi*tm_1x6_to_3x3(bt); % 3x3
    product_tensors = tm_3x3_to_1x6((eye(3) + product_Psi_B)\Psi*(H_inv - kappa*eye(3))); % 1x6
    s(i) = s0*det(eye(3) + product_Psi_B)^(-kappa)*exp(-tm_inner(bt,product_tensors));
%     s0
%     det(eye(3) + product_Psi_B)
%     kappa
%     det(eye(3) + product_Psi_B)^(-kappa)
%     exp(-tm_inner(bt,product_tensors))
end

% Force signal to be real. This is necessary when allowing negative variances.
s = real(s);

