function [tensor_xx,tensor_yy,tensor_zz,tensor_xy,tensor_xz,tensor_yz] = dtd_mv_gamma_multidim_Euler_rotation(tensor,alpha,beta,gamma)
% dtd_mv_gamma_multidim_Euler_rotation(tensor,alpha,beta,gamma)

% [nx,ny,nz,~] = size(tensor);
% tensor_rotated = zeros([nx,ny,nz,6]);

tensor_X = tensor(:,:,:,1);
tensor_Y = tensor(:,:,:,2);
tensor_Z = tensor(:,:,:,3);

R_euler_11 = cos(alpha).*cos(beta).*cos(gamma) - sin(alpha).*sin(gamma);
R_euler_21 = cos(alpha).*cos(beta).*sin(gamma) + sin(alpha).*cos(gamma);
R_euler_22 = -sin(alpha).*cos(beta).*sin(gamma) + cos(alpha).*cos(gamma);
R_euler_12 = -sin(alpha).*cos(beta).*cos(gamma) - cos(alpha).*sin(gamma);
R_euler_13 = sin(beta).*cos(gamma);
R_euler_23 = sin(beta).*sin(gamma);
R_euler_31 = -cos(alpha).*sin(beta);
R_euler_32 = sin(alpha).*sin(beta);
R_euler_33 = cos(beta);


tensor_xx = tensor_X.*R_euler_11.^2 + tensor_Y.*R_euler_12.^2 + tensor_Z.*R_euler_13.^2;
tensor_yy = tensor_X.*R_euler_21.^2 + tensor_Y.*R_euler_22.^2 + tensor_Z.*R_euler_23.^2;
tensor_zz = tensor_X.*R_euler_31.^2 + tensor_Y.*R_euler_32.^2 + tensor_Z.*R_euler_33.^2;
tensor_xy = tensor_X.*R_euler_11.*R_euler_21 + tensor_Y.*R_euler_12.*R_euler_22 + tensor_Z.*R_euler_13.*R_euler_23;
tensor_xz = tensor_X.*R_euler_11.*R_euler_31 + tensor_Y.*R_euler_12.*R_euler_32 + tensor_Z.*R_euler_13.*R_euler_33;
tensor_yz = tensor_X.*R_euler_21.*R_euler_31 + tensor_Y.*R_euler_22.*R_euler_32 + tensor_Z.*R_euler_23.*R_euler_33;