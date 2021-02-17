function tensor_SPD = dtd_mv_gamma_multidim_nearest_SPD(tensor)

% Here tensor is sized [nx,ny,nz,3,3]

[nx,ny,nz,~,~] = size(tensor);
tensor_SPD = zeros(size(tensor));

for vx = 1:nx
    for vy = 1:ny
        for vz = 1:nz
            [~,p] = chol(squeeze(tensor(vx,vy,vz,:,:)));
            if p ~= 0
                tensor_SPD(vx,vy,vz,:,:) = dtd_mv_gamma_nearest_SPD(squeeze(tensor(vx,vy,vz,:,:)));
            end
        end
    end
end










