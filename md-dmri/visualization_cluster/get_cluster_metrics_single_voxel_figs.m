function get_cluster_metrics_single_voxel_figs()

close all

% CC body
% char_mask = 'CC_body';
% vx = 39;
% vy = 44;
% vz = 14;
% ylim_diso = [0.6 1]*1e-9;
% ylim_sqddelta = [0.7 0.9];
% ylim_r = [15 25];
% ylim_t = [50 100]*1e-3;
% ylim_w = [0 1];
% unit_factor = [1e9, 1, 1, 1e3, 1];
% FLAG_ylim = 1;

% CC splenium (back)
% char_mask = 'CC_splenium';
% vx = 39;
% vy = 32;
% vz = 9;
% ylim_diso = [0.6 1]*1e-9;
% ylim_sqddelta = [0.7 0.9];
% ylim_r = [15 25];
% ylim_t = [50 100]*1e-3;
% ylim_w = [0 1];
% unit_factor = [1e9, 1, 1, 1e3, 1];
% FLAG_ylim = 1;

% CC genu (front)
% char_mask = 'CC_genu';
% vx = 39;
% vy = 53;
% vz = 8;
% ylim_diso = [0.6 1]*1e-9;
% ylim_sqddelta = [0.7 0.9];
% ylim_r = [15 25];
% ylim_t = [50 100]*1e-3;
% ylim_w = [0 1];
% unit_factor = [1e9, 1, 1, 1e3, 1];
% FLAG_ylim = 1;

% CS
% char_mask = 'CS';
% vx = 48;
% vy = 39;
% vz = 16;
% ylim_diso = [0.7 1]*1e-9;
% ylim_sqddelta = [0.8 0.9];
% ylim_r = [10 20];
% ylim_t = [55 110]*1e-3;
% ylim_w = [0 1];
% unit_factor = [1e9, 1, 1, 1e3, 1];
% FLAG_ylim = 1;

% CC_CING
char_mask = 'CC_CING';
vx = 37;
vy = 50;
vz = 13;
ylim_diso = [0.5 1]*1e-9;
ylim_sqddelta = [0.6 0.9];
ylim_r = [15 30];
ylim_t = [40 70]*1e-3;
ylim_w = [0 1];
unit_factor = [1e9, 1, 1, 1e3, 1];
FLAG_ylim = 1;


initial_directory = 'Data_Chantal';
method = 'dtr2d';

% initial_directory = 'Data_T1_D';
% method = 'dtr1d';

directory_odf = fullfile(pwd, initial_directory, '3_odfs');

directory_input = fullfile(pwd, initial_directory, '4_clustering');
directory_output = fullfile(pwd, initial_directory, 'cluster_metrics_single_voxel');

if ~exist(directory_output, 'dir')
    msf_mkdir(directory_output)
end

load(fullfile(directory_input, 'cluster_metrics.mat'), 'map_cluster_metrics')
load(fullfile(directory_odf, 'ODF_peaks_metrics.mat'), 'cell_odf_peaks')

mask = zeros(size(cell_odf_peaks));
mask(vx,vy,vz) = 1;
mdm_nii_write(mask, fullfile(directory_output, strcat('mask_', char_mask, '.nii.gz')));

odf_metrics = cell_odf_peaks{vx,vy,vz};
orientation_peaks = odf_metrics.orientation_peaks;
diso_peaks = odf_metrics.diso_peaks*1e9;
sddelta_peaks = odf_metrics.sddelta_peaks;
if strcmp(method,'dtr2d')
    r_peaks = odf_metrics.r2_peaks;
    t_peaks = odf_metrics.t2_peaks*1e3;
elseif strcmp(method,'dtr1d')
    r_peaks = odf_metrics.r1_peaks;
    t_peaks = odf_metrics.t1_peaks*1e3;
end

metrics_voxel = map_cluster_metrics{vx,vy,vz};
N = metrics_voxel.n;
x = metrics_voxel.x;
y = metrics_voxel.y;
z = metrics_voxel.z;
dpar = metrics_voxel.dpar;
dperp = metrics_voxel.dperp;
w = metrics_voxel.w;
if strcmp(method,'dtr2d')
    r = metrics_voxel.r2;
    t = metrics_voxel.t2;
    char_r = '$\mathring{\mathrm{E}}[R_2]\,(\mathrm{s}^{-1})$';
    char_t = '$\mathring{\mathrm{E}}[T_2]\,(\mathrm{ms})$';
elseif strcmp(method,'dtr1d')
    r = metrics_voxel.r1;
    t = metrics_voxel.t1;
    char_r = '$\mathring{\mathrm{E}}[R_1]\,(\mathrm{s}^{-1})$';
    char_t = '$\mathring{\mathrm{E}}[T_1]\,(\mathrm{s})$';
end
median_orientations = metrics_voxel.median_orientations;

w_sum = zeros([1 N]);
for m = 1:N
    w_sum(m) = sum(w{m});
end

diso = cell([1 N]);
sqddelta = cell([1 N]);
colors = cell([1 N]);
relevant_orientation_peaks = zeros([N, 3]);
relevant_diso_peaks = zeros([1 N]);
relevant_sqddelta_peaks = zeros([1 N]);
relevant_r_peaks = zeros([1 N]);
relevant_t_peaks = zeros([1 N]);

max_w = 0;
for m = 1:N
   diso{m} = (dpar{m}+2*dperp{m})/3;
   sqddelta{m} = ((dpar{m}-dperp{m})./(dpar{m}+2*dperp{m})).^2;
   median_xyz = median_orientations(3*(m-1)+1:3*(m-1)+3);
   color_cell = cluster_orientational_colors(median_xyz(1), median_xyz(2), median_xyz(3));
   colors{m} = color_cell{1};
   potential_max_w = max(w{m});
   if potential_max_w > max_w
       max_w = potential_max_w;
   end
   
   relevant_dots = (orientation_peaks*median_xyz')';
   [max_dot, ind_max] = max(abs(relevant_dots));
   if acos(abs(max_dot))*180/pi < 30
       relevant_orientation_peaks(m,:) = orientation_peaks(ind_max,:);
       relevant_diso_peaks(m) = diso_peaks(ind_max);
       relevant_sqddelta_peaks(m) = sddelta_peaks(ind_max);
       relevant_r_peaks(m) = r_peaks(ind_max);
       relevant_t_peaks(m) = t_peaks(ind_max);
   end
   
end

%% PLOT
param_names = {'diso'; 'sqddelta'; 'r'; 't' ; 'w'};
param_names_plot = {'$\mathring{\mathrm{E}}[D_\mathrm{iso}]\,(\mu\mathrm{m}^2/\mathrm{ms})$'; '$\mathring{\mathrm{E}}[D_\Delta^2]$'; char_r ; char_t ; '$\mathring{w}$'};
Nparam = numel(param_names);

set(0, 'defaultLegendInterpreter','latex');
lw = 2;
global_font_size = 16;
global_font_size_labels = 17;

width_sphere = 0.25;
height_sphere = 0.25;
width_cluster_metrics = 0.17;
height_cluster_metrics = 0.17;

inter_sphere = 0.15;
inter_h = 0.2;
inter_v = 0.01;
inter_v_2 = 0.04;

radius_plotted_sphere = 0.98;

error_cap = 10;

offset_bottom = 0.03;

y_bottom = offset_bottom + (1 - 2*inter_v - inter_v_2 - 3*height_cluster_metrics - height_sphere)/2;
y1 = y_bottom;
y2 = y1 + inter_v + height_cluster_metrics;
y3 = y2 + inter_v + height_cluster_metrics;
y_top = y3 + inter_v_2 + height_cluster_metrics;

x_top_left = (1 - 2*width_sphere - inter_sphere)/2;
x_top_right = x_top_left + inter_sphere + width_sphere;
x_left = (1 - 2*width_cluster_metrics - inter_h)/2;
x_right = x_left + width_cluster_metrics + inter_h;
x_mid = (1 - width_cluster_metrics)/2;
    
f = figure('Position', [0,0,1500,1500], 'Units', 'pixels', 'visible', 'off');
% f.PaperOrientation = 'landscape';
f.PaperUnits = 'normalized';
f.PaperPosition = [0 0 1 1];
set(gcf, 'Renderer', 'painters')

axh1 = axes('position',[x_top_left y_top width_sphere height_sphere]);
axh2 = axes('position',[x_top_right y_top width_sphere height_sphere]);
axh_diso = axes('position',[x_left y3 width_cluster_metrics height_cluster_metrics]);
axh_sqddelta = axes('position',[x_left y2 width_cluster_metrics height_cluster_metrics]);
axh_r = axes('position',[x_right y3 width_cluster_metrics height_cluster_metrics]);
axh_t = axes('position',[x_right y2 width_cluster_metrics height_cluster_metrics]);
axh_w = axes('position',[x_mid y1 width_cluster_metrics height_cluster_metrics]);

%% Orientations
hold(axh1, 'on')
plot3(axh1, [0, 0; -radius_plotted_sphere, radius_plotted_sphere; 0, 0]',[0, 0; 0, 0; -radius_plotted_sphere, radius_plotted_sphere]',[-radius_plotted_sphere, radius_plotted_sphere; 0, 0; 0, 0]', '-k', 'LineWidth',1);
[x_sph,y_sph,z_sph] = sphere(50);
hold(axh1, 'on')
surf(axh1, radius_plotted_sphere*x_sph, radius_plotted_sphere*y_sph, radius_plotted_sphere*z_sph, 'FaceColor', 'black', 'EdgeColor', 'black', 'FaceAlpha', 0.1,'EdgeAlpha', 0.2);
hold(axh1, 'on')

hold(axh2, 'on')
plot3(axh2, [0, 0; -radius_plotted_sphere, radius_plotted_sphere; 0, 0]',[0, 0; 0, 0; -radius_plotted_sphere, radius_plotted_sphere]',[-radius_plotted_sphere, radius_plotted_sphere; 0, 0; 0, 0]', '-k', 'LineWidth',1);
[x_sph,y_sph,z_sph] = sphere(50);
hold(axh2, 'on')
surf(axh2, radius_plotted_sphere*x_sph, radius_plotted_sphere*y_sph, radius_plotted_sphere*z_sph, 'FaceColor', 'black', 'EdgeColor', 'black', 'FaceAlpha', 0.1,'EdgeAlpha', 0.2);
hold(axh2, 'on')

for m = 1:N
    color_to_plot = colors{m};
    x_fibers_to_plot = x{m};
    y_fibers_to_plot = y{m};
    z_fibers_to_plot = z{m};
    w_fibers_to_plot = w{m};
    for i = 1:length(x_fibers_to_plot)
        scatter3(axh1, x_fibers_to_plot(i), y_fibers_to_plot(i), z_fibers_to_plot(i), 20, 'filled', 'MarkerFaceAlpha', w_fibers_to_plot(i)/max_w, 'MarkerEdgeAlpha', w_fibers_to_plot(i)/max_w, 'MarkerFaceColor', color_to_plot, 'MarkerEdgeColor', color_to_plot);
        hold(axh1, 'on')
        scatter3(axh2, x_fibers_to_plot(i), y_fibers_to_plot(i), z_fibers_to_plot(i), 20, 'filled', 'MarkerFaceAlpha', w_fibers_to_plot(i)/max_w, 'MarkerEdgeAlpha', w_fibers_to_plot(i)/max_w, 'MarkerFaceColor', color_to_plot, 'MarkerEdgeColor', color_to_plot);
        hold(axh2, 'on')
    end
end

axis(axh1, 'equal')
set(axh1, 'LineWidth', lw)
set(axh1, 'TickLength', [0.04 0.04], 'TickDir','out');
set(axh1,'FontSize', global_font_size, 'TickLabelInterpreter', 'latex')
xlabel(axh1, '$x$', 'FontSize', global_font_size_labels, 'Interpreter', 'latex')
ylabel(axh1, '$y$', 'FontSize', global_font_size_labels, 'Interpreter', 'latex')
zlabel(axh1, '$z$', 'FontSize', global_font_size_labels, 'Interpreter', 'latex')
xlim(axh1, [-1 1])
ylim(axh1, [-1 1])
zlim(axh1, [-1 1])
view(axh1, 142.5, 20);

axis(axh2, 'equal')
set(axh2, 'LineWidth', lw)
set(axh2, 'TickLength', [0.04 0.04], 'TickDir','out');
set(axh2,'FontSize', global_font_size, 'TickLabelInterpreter', 'latex')
xlabel(axh2, '$x$', 'FontSize', global_font_size_labels, 'Interpreter', 'latex')
ylabel(axh2, '$y$', 'FontSize', global_font_size_labels, 'Interpreter', 'latex')
zlabel(axh2, '$z$', 'FontSize', global_font_size_labels, 'Interpreter', 'latex')
xlim(axh2, [-1 1])
ylim(axh2, [-1 1])
zlim(axh2, [-1 1])
view(axh2, 142.5, 90);

%% Metrics
for nparam = 1:Nparam
    eval(['axh = axh_' param_names{nparam} ';'])
    eval(strcat('metric =', param_names{nparam}, ';'))
    hold(axh, 'on')
    
    if ~strcmp(param_names{nparam}, 'w')
        eval(['peak_metric = relevant_' param_names{nparam} '_peaks;'])
    end
    
    medians = zeros(1,N);
    Q1s = zeros(1,N);
    Q3s = zeros(1,N);
    lower = zeros(1,N);
    upper = zeros(1,N);
    
    max_length = 0;
    repeated_metric = cell([1 N]);
    for m = 1:N
        color = colors{m};
        if strcmp(param_names{nparam}, 'w')
            weighted_quantiles = weighted_quantile(unit_factor(nparam)*metric{m}, ones(size(w{m})), [0.25 0.5 0.75]);
            display(['Total weight = ' num2str(sum(unit_factor(nparam)*metric{m}))]);
        else
            weighted_quantiles = weighted_quantile(unit_factor(nparam)*metric{m}, w{m}, [0.25 0.5 0.75]);
        end
        Q1s(m) = weighted_quantiles(1);
        medians(m) = weighted_quantiles(2);
        Q3s(m) = weighted_quantiles(3);
        lower(m) = medians(m) - Q1s(m);
        upper(m) = Q3s(m) - medians(m);
        errorbar(axh, m, medians(m), zeros(size(lower(m))), zeros(size(upper(m))), 'o', 'MarkerFaceColor', color, 'Capsize', 0, 'Linewidth', 3*lw/4, 'Color', color)
        errorbar(axh, m, medians(m), lower(m), upper(m), '', 'Capsize', error_cap, 'Linewidth', 3*lw/4, 'Color', color)
        
        if ~strcmp(param_names{nparam}, 'w') && peak_metric(m)
            errorbar(axh, m, peak_metric(m), zeros(size(median(m))), zeros(size(median(m))), 's', 'MarkerSize', 5, 'MarkerFaceColor', color, 'MarkerEdgeColor', 'k', 'Capsize', 0, 'Linewidth', lw/4)
        end
        
        w_integer = w{m};
        while ~all(floor(w_integer) >= 1)
            w_integer = 10*w_integer; 
        end
        repeated_metric{m} = repelem(unit_factor(nparam)*metric{m}, floor(w_integer));
        if length(repeated_metric{m}) > max_length
            max_length = length(repeated_metric{m});
        end
    end

    max_length = max(length(repeated_metric{m}));
    matrix_Kruskal_Wallis = NaN*zeros(max_length, N);
    for m = 1:N
        matrix_Kruskal_Wallis(1:length(repeated_metric{m}),m) = repeated_metric{m};
    end
    
    fprintf([param_names{nparam} ':\n'])
    p = kruskalwallis(matrix_Kruskal_Wallis,[],'off');
    fprintf(['Kruskal_Wallis: ' num2str(p) '\n'])
    for m = 1:N
        for m2 = (m+1):N
            max_length = max([length(repeated_metric{m}),length(repeated_metric{m2})]);
            matrix_Kruskal_Wallis = NaN*zeros(max_length, 2);
            matrix_Kruskal_Wallis(1:length(repeated_metric{m}),1) = repeated_metric{m};
            matrix_Kruskal_Wallis(1:length(repeated_metric{m2}),2) = repeated_metric{m2};
            p = kruskalwallis(matrix_Kruskal_Wallis,[],'off');
            
            fprintf(['(' num2str(m) ',' num2str(m2) '): \n'])
            fprintf(['Kruskal_Wallis: ' num2str(p) '\n'])
            p = ranksum(repeated_metric{m},repeated_metric{m2});
            fprintf(['Mann-Whitney: ' num2str(p) '\n'])
        end
    end
    fprintf('\n')
    
    box(axh, 'on')
    grid(axh, 'on')
    set(axh, 'XTick', 0:N+1, 'XTickLabel', {' '})
    set(axh,'FontSize', global_font_size, 'TickLabelInterpreter', 'latex')
    set(axh, 'LineWidth', lw)
    set(axh, 'TickLength', [0.04 0.035], 'TickDir', 'out');
    xlim(axh, [0, N+1])
    
    if FLAG_ylim
        eval(['y_limits =ylim_' param_names{nparam} ';'])
        ylim(axh, unit_factor(nparam)*y_limits)
    end
    
    axis(axh, 'square')
    
    ylab = ylabel(axh, param_names_plot{nparam}, 'Fontsize', global_font_size_labels, 'Interpreter', 'latex');
    set(ylab, 'Units', 'Normalized', 'Position', [-0.5, 0.5, 0]);  

    if nparam == Nparam
        xlab = xlabel(axh, 'Clusters', 'Fontsize', global_font_size_labels, 'Interpreter', 'latex');
        set(xlab, 'Units', 'Normalized', 'Position', [0.5, -0.2, 0]);
    end
    
    
%     
%     if FLAG_ylim
% %         if strcmp(param_names{nparam},"R2")
% %             ylim(axh, [min_ylim_T2 max_ylim_T2]);
% %         else
% %             eval(strcat('ylim(axh, [min_ylim_', param_names{nparam}, ' max_ylim_', param_names{nparam},']);'))
% %         end
%         
%         eval(strcat('ylim(axh, [min_ylim_', param_names{nparam}, ' max_ylim_', param_names{nparam},']);'))
%         
%     end
%     
% %     if nparam == 1 && it_SNR == 2
% %         [h_legend, hObj] = legend(axh, L, legend_names);
% %         pos = get(h_legend,'Position');
% %         set(h_legend,'Position', [x1+0.017, inter_v_top_bottom+height_polar/2, pos(3:4)]);
% %         hl = findobj(hObj,'type','line');
% %         set(hl,'LineWidth',1.5);
% %     end
% end
% 
% 

end



saveas(gcf, fullfile(directory_output, strcat('cluster_metrics_vx', num2str(vx), '_vy', num2str(vy), '_vz', num2str(vz), '.pdf')))
clf(f)




end