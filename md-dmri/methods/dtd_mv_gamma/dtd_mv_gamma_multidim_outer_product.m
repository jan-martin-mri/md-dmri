function tensor_outer_product = dtd_mv_gamma_multidim_outer_product(tensor_1, tensor_2)
% function dtd_mv_gamma_multidim_kronecker_product(tensor_1, tensor_2)

% Outputs 3x3 tensors (assuming tensor_1 and tensor_2 share the same eigenvectors)

[nx,ny,nz,~] = size(tensor_1);
tensor_outer_product = zeros([nx,ny,nz,3,3]);

tensor_outer_product(:,:,:,1,1) = tensor_1(:,:,:,1).*tensor_2(:,:,:,1);
tensor_outer_product(:,:,:,1,2) = tensor_1(:,:,:,1).*tensor_2(:,:,:,2);
tensor_outer_product(:,:,:,1,3) = tensor_1(:,:,:,1).*tensor_2(:,:,:,3);
tensor_outer_product(:,:,:,2,1) = tensor_1(:,:,:,2).*tensor_2(:,:,:,1);
tensor_outer_product(:,:,:,2,2) = tensor_1(:,:,:,2).*tensor_2(:,:,:,2);
tensor_outer_product(:,:,:,2,3) = tensor_1(:,:,:,2).*tensor_2(:,:,:,3);
tensor_outer_product(:,:,:,3,1) = tensor_1(:,:,:,3).*tensor_2(:,:,:,1);
tensor_outer_product(:,:,:,3,2) = tensor_1(:,:,:,3).*tensor_2(:,:,:,2);
tensor_outer_product(:,:,:,3,3) = tensor_1(:,:,:,3).*tensor_2(:,:,:,3);