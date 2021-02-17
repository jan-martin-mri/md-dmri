function plot_fiber_sphere(x_fibers, y_fibers, z_fibers, w_fibers, ind_fibers_in_relevant_cluster, ind_fibers_not_in_relevant_cluster, optimal_K, colors_cluster, figure_path, vx, vy, vz)

set(0, 'defaultLegendInterpreter','latex');

max_w = max(w_fibers);
nb_fibers = length(w_fibers);
nb_fibers_out_of_clustering = length(ind_fibers_not_in_relevant_cluster);

%color_out_of_clustering = [round(184/255,1), 0, round(230/255,1)]; % Purple
color_out_of_clustering = [0, 0, 0];
radius_plotted_sphere = 0.98;

f = figure('Position',[10 10 300 650],'visible','off', 'renderer','Painters');

sub = subplot(2,1,1);
lines_axes = plot3([0, 0; -radius_plotted_sphere, radius_plotted_sphere; 0, 0]',[0, 0; 0, 0; -radius_plotted_sphere, radius_plotted_sphere]',[-radius_plotted_sphere, radius_plotted_sphere; 0, 0; 0, 0]', '-k', 'LineWidth',1);
hold on
[x_sph,y_sph,z_sph] = sphere(50);
sphere_plot = surf(radius_plotted_sphere*x_sph, radius_plotted_sphere*y_sph, radius_plotted_sphere*z_sph, 'FaceColor', 'black', 'EdgeColor', 'black', 'FaceAlpha', 0.1,'EdgeAlpha', 0.2);
hold on
points_cluster = cell([1,nb_fibers]);
count = 1;
for m = 1:optimal_K
    color_to_plot = colors_cluster{m};
    for i = ind_fibers_in_relevant_cluster{m}
        points_cluster{count} = scatter3(x_fibers(i), y_fibers(i), z_fibers(i), 20, 'filled', 'MarkerFaceAlpha', w_fibers(i)/max_w, 'MarkerEdgeAlpha', w_fibers(i)/max_w, 'MarkerFaceColor', color_to_plot, 'MarkerEdgeColor', color_to_plot);
        hold on
        count = count + 1;
    end
end

points_out_of_cluster = cell([1,nb_fibers_out_of_clustering]);
count = 1;
for i = ind_fibers_not_in_relevant_cluster
    points_out_of_cluster{count} = scatter3(x_fibers(i), y_fibers(i), z_fibers(i), 20, 'filled', 'MarkerFaceAlpha', w_fibers(i)/max_w, 'MarkerEdgeAlpha', w_fibers(i)/max_w, 'MarkerFaceColor', color_out_of_clustering, 'MarkerEdgeColor', color_out_of_clustering);
    hold on
    count = count + 1;
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
for m = 1:length(points_out_of_cluster)
    copyobj(points_out_of_cluster{m},sub);
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
saveas(f, [figure_path '/' num2str(vx) '_' num2str(vy) '_' num2str(vz) '_fiber_clusters_all_bootstrap.pdf'])
hold off