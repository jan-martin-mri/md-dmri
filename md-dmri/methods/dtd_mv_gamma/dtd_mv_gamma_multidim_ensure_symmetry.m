function tensor = dtd_mv_gamma_multidim_ensure_symmetry(tensor)
% function dtd_mv_gamma_multidim_kronecker_product(tensor_1, tensor_2)

tensor(:,:,:,1,2) = (tensor(:,:,:,1,2)+tensor(:,:,:,2,1))/2;
tensor(:,:,:,1,3) = (tensor(:,:,:,1,3)+tensor(:,:,:,3,1))/2;
tensor(:,:,:,2,3) = (tensor(:,:,:,2,3)+tensor(:,:,:,3,2))/2;

tensor(:,:,:,2,1) = tensor(:,:,:,1,2);
tensor(:,:,:,3,1) = tensor(:,:,:,1,3);
tensor(:,:,:,3,2) = tensor(:,:,:,2,3);
