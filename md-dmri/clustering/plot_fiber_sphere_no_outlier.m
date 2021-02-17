function plot_fiber_sphere_no_outlier(x_fibers, y_fibers, z_fibers, w_fibers, optimal_K, colors_cluster, figure_path, vx, vy, vz, do_HD, suffix)

set(0, 'defaultLegendInterpreter','latex');

if nargin < 11
    do_HD = 1;
end

if nargin < 12
   suffix = '_fiber_clusters'; 
end

all_w = [];
for m = 1:optimal_K
    all_w = [all_w w_fibers{m}];
end
max_w = max(all_w);
nb_fibers = length(all_w);

radius_plotted_sphere = 0.98;

if do_HD
    %f = figure('Position',[10 10 300 650], 'Units', 'pixels','visible','off','renderer','Painters');
    f = figure('Position',[10 10 300 650], 'Units', 'pixels','renderer','Painters');
else
    %f = figure('Position',[10 10 300 650], 'Units', 'pixels','visible','off');
    f = figure('Position',[10 10 300 650], 'Units', 'pixels');
end

sub = subplot(2,1,1);
lines_axes = plot3([0, 0; -radius_plotted_sphere, radius_plotted_sphere; 0, 0]',[0, 0; 0, 0; -radius_plotted_sphere, radius_plotted_sphere]',[-radius_plotted_sphere, radius_plotted_sphere; 0, 0; 0, 0]', '-k', 'LineWidth',1);
hold on
[x_sph,y_sph,z_sph] = sphere(50);
sphere_plot = surf(radius_plotted_sphere*x_sph, radius_plotted_sphere*y_sph, radius_plotted_sphere*z_sph, 'FaceColor', 'black', 'EdgeColor', 'black', 'FaceAlpha', 0.1,'EdgeAlpha', 0.2);
hold on
points_cluster = cell([1,nb_fibers]);
count = 1;
for m = 1:optimal_K
    %color_to_plot = colors_cluster((3*(m-1)+1):((3*(m-1)+1)+2));
    color_to_plot = colors_cluster{m};
    x_fibers_to_plot = x_fibers{m};
    y_fibers_to_plot = y_fibers{m};
    z_fibers_to_plot = z_fibers{m};
    w_fibers_to_plot = w_fibers{m};
    for i = 1:length(x_fibers_to_plot)
        points_cluster{count} = scatter3(x_fibers_to_plot(i), y_fibers_to_plot(i), z_fibers_to_plot(i), 20, 'filled', 'MarkerFaceAlpha', w_fibers_to_plot(i)/max_w, 'MarkerEdgeAlpha', w_fibers_to_plot(i)/max_w, 'MarkerFaceColor', color_to_plot, 'MarkerEdgeColor', color_to_plot);
        hold on
        count = count + 1;
    end
end

axis equal
sub.XAxis.TickLabelInterpreter = 'latex';
sub.YAxis.TickLabelInterpreter = 'latex';
sub.ZAxis.TickLabelInterpreter = 'latex';
set(gca,'FontSize',16)
xlabel('$x$', 'FontSize', 18, 'Interpreter', 'latex')
ylabel('$y$', 'FontSize', 18, 'Interpreter', 'latex')
zlabel('$z$', 'FontSize', 18, 'Interpreter', 'latex')
xlim([-1 1])
ylim([-1 1])
zlim([-1 1])
view(142.5, 90);

sub = subplot(2,1,2);
copyobj(lines_axes,sub);
hold on
copyobj(sphere_plot,sub);
hold on
for m = 1:length(points_cluster)
    copyobj(points_cluster{m},sub);
    hold on
end

axis equal
sub.XAxis.TickLabelInterpreter = 'latex';
sub.YAxis.TickLabelInterpreter = 'latex';
sub.ZAxis.TickLabelInterpreter = 'latex';
set(gca,'FontSize',16)
xlabel('$x$', 'FontSize', 18, 'Interpreter', 'latex')
ylabel('$y$', 'FontSize', 18, 'Interpreter', 'latex')
zlabel('$z$', 'FontSize', 18, 'Interpreter', 'latex')
xlim([-1 1])
ylim([-1 1])
zlim([-1 1])
view(142.5, 20);

% print(f, [figure_path '/' num2str(vx) '_' num2str(vy) '_' num2str(vz) '_fiber_clusters.pdf'], '-dpdf', '-bestfit', '-r500')
saveas(f, [figure_path '/' num2str(vx) '_' num2str(vy) '_' num2str(vz) , suffix, '.pdf'])
hold off