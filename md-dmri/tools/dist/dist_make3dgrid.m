function [x_coord,y_coord,z_coord] = dist_make3dgrid(x_v,y_v,z_v)

% Draw x lines
x_basis = linspace(min(x_v),max(x_v),10);
y_basis = y_v;
z_basis = z_v;

[x_grid,y_grid,z_grid] = ndgrid(x_basis,y_basis,z_basis);
sz = size(x_grid);
x_grid = reshape(x_grid,sz(1),sz(2)*sz(3));
y_grid = reshape(y_grid,sz(1),sz(2)*sz(3));
z_grid = reshape(z_grid,sz(1),sz(2)*sz(3));

x_coord = x_grid;
y_coord = y_grid;
z_coord = z_grid;

% hold on
% plot3(x_coord,y_coord,z_coord,'Color',[0 0 0],'LineWidth',lw_axes)

% Draw y lines
x_basis = x_v;
y_basis = linspace(min(y_v),max(y_v),10);
z_basis = z_v;

[y_grid,x_grid,z_grid] = ndgrid(y_basis,x_basis,z_basis);
sz = size(x_grid);
x_grid = reshape(x_grid,sz(1),sz(2)*sz(3));
y_grid = reshape(y_grid,sz(1),sz(2)*sz(3));
z_grid = reshape(z_grid,sz(1),sz(2)*sz(3));

% hold on
% plot3(x_coord,y_coord,z_coord,'Color',[0 0 0],'LineWidth',lw_axes)

x_coord = [x_coord x_grid];
y_coord = [y_coord y_grid];
z_coord = [z_coord z_grid];

% Draw z lines
x_basis = x_v;
y_basis = y_v;
z_basis = linspace(min(z_v),max(z_v),10);

[z_grid,x_grid,y_grid] = ndgrid(z_basis,x_basis,y_basis);
sz = size(x_grid);
x_grid = reshape(x_grid,sz(1),sz(2)*sz(3));
y_grid = reshape(y_grid,sz(1),sz(2)*sz(3));
z_grid = reshape(z_grid,sz(1),sz(2)*sz(3));

% hold on
% plot3(x_coord,y_coord,z_coord,'Color',[0 0 0],'LineWidth',lw_axes)

x_coord = [x_coord x_grid];
y_coord = [y_coord y_grid];
z_coord = [z_coord z_grid];