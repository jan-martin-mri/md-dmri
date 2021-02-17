function [x_coord,y_coord,z_coord,color] = dist_draw_polygon(x_v,y_v,z_v,color)

if nargin > 3
    color = reshape(color,size(color,1),1,3);
end

x_coord = [];
y_coord = [];
z_coord = [];

for i = 1:numel(z_v)
    
    x_basis = [x_v(1) x_v(2) x_v(2) x_v(1)];
    y_basis = [y_v(1) y_v(1) y_v(2) y_v(2)];
    z_basis = z_v(i)*ones(1,4);% y_v(1)];
    
%     [x_grid,y_grid,z_grid] = ndgrid(x_basis,y_basis,z_basis);
%     
%     x_coord = x_basis]; %reshape(x_grid,numel(x_grid),1)];
%     y_coord = [y_coord y_basis]; %reshape(y_grid,numel(y_grid),1)];
%     z_coord = [z_coord z_basis]; %reshape(z_grid,numel(z_grid),1)];
    
    x_coord = [x_coord; x_basis];
    y_coord = [y_coord; y_basis];
    z_coord = [z_coord; z_basis];

end