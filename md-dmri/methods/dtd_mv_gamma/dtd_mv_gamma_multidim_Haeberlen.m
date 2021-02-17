function [tensor_delta, tensor_eta, ind_Haeberlen] = dtd_mv_gamma_multidim_Haeberlen(tensor,tensor_iso)
% function dtd_mv_gamma_multidim_Haeberlen(tensor,tensor_iso)

% Outputs 3x3 tensors (assuming tensor and tensor_iso are nx x ny x nz x 3)

nnz(tensor)

[~, ind_Haeberlen] = sort(abs(tensor - tensor_iso),4);
tensor_Haeberlen = tensor(ind_Haeberlen);

nnz(tensor_Haeberlen)

tensor_delta = 1./(3*tensor_iso).*(tensor_Haeberlen(:,:,:,3)-(tensor_Haeberlen(:,:,:,2)+tensor_Haeberlen(:,:,:,1))/2);
tensor_eta = (tensor_Haeberlen(:,:,:,2)-tensor_Haeberlen(:,:,:,1))./(3*tensor_iso.*tensor_delta);