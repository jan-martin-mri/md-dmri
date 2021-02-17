close all

set(0, 'defaultLegendInterpreter','latex');

%% General
% method = 'dtr2d';
% initial_directory = 'Data_test';
% peak_map = mdm_nii_read(fullfile(pwd, initial_directory, '4_clustering', 'cluster_median_orientations_weighted_by_w.nii.gz'));
% metric_map = mdm_nii_read(fullfile(pwd, initial_directory, '2_nii_pdf_maps', strcat(method, '_fraction_not_thin.nii.gz'))); 
% fraction_map = mdm_nii_read(fullfile(pwd, initial_directory, '2_nii_pdf_maps', strcat(method, '_fraction_thin.nii.gz'))); 
% mask = mdm_nii_read(fullfile(pwd, initial_directory, '0_data', 'b0_bet_mask_data.nii.gz')); 
% output_directory = fullfile(pwd, initial_directory, 'clustering_peak_maps');
% nb_parallel_workers = 6;

method = 'dtr1d';
initial_directory = 'Data_T1_D';
peak_map = mdm_nii_read(fullfile(pwd, initial_directory, '4_clustering', 'cluster_median_orientations_weighted_by_w.nii.gz'));
metric_map = mdm_nii_read(fullfile(pwd, initial_directory, '2_nii_pdf_maps', strcat(method, '_fraction_not_thin.nii.gz'))); 
fraction_map = mdm_nii_read(fullfile(pwd, initial_directory, '2_nii_pdf_maps', strcat(method, '_fraction_thin.nii.gz'))); 
mask = mdm_nii_read(fullfile(pwd, initial_directory, '0_data', 'data_mask.nii.gz')); 
output_directory = fullfile(pwd, initial_directory, 'clustering_peak_maps');
nb_parallel_workers = 6;

% method = 'dtd';
% initial_directory = 'Data_tumor_Spectrum';
% peak_map = mdm_nii_read(fullfile(pwd, initial_directory, '4_clustering', 'cluster_median_orientations_weighted_by_w.nii.gz'));
% metric_map = mdm_nii_read(fullfile(pwd, initial_directory, '2_nii_pdf_maps', strcat(method, '_fraction_not_thin.nii.gz'))); 
% fraction_map = mdm_nii_read(fullfile(pwd, initial_directory, '2_nii_pdf_maps', strcat(method, '_fraction_thin.nii.gz'))); 
% mask = mdm_nii_read(fullfile(pwd, initial_directory, '0_data', 'data_mc_low_b_bet_mask_dilated.nii.gz')); 
% output_directory = fullfile(pwd, initial_directory, 'clustering_peak_maps');
% nb_parallel_workers = 6;

do_background = 1;
use_other_colored_peak_metric = 1;
do_IQR = 0;
colored_metrics_char = 'r1';

do_coronal = 1;
do_axial = 1;
do_sagittal = 1;

mask = double(mask);

if do_background
    metric_map = metric_map.*mask;
else
    metric_map = zeros(size(metric_map));
end

[nb_voxel_x, nb_voxel_y, nb_voxel_z, N] = size(peak_map);
max_nb_peaks = N/3;


%% Other colored metrics along peaks?
if use_other_colored_peak_metric
    load(fullfile(pwd, initial_directory, '4_clustering', 'cluster_metrics.mat')); % Creates map_cluster_metrics
    
    if do_IQR
        colored_metrics_char_plot = ['IQR_' colored_metrics_char '_'];
        if strcmp(colored_metrics_char, 'r2')
            bound_min = 0;
            bound_max = 5;
        elseif strcmp(colored_metrics_char, 't2')
            bound_min = 0;
            bound_max = 0.02;
        elseif strcmp(colored_metrics_char, 'diso')
            bound_min = 0*1e-9;
            bound_max = 0.5*1e-9;
        elseif strcmp(colored_metrics_char, 'sqddelta')
            bound_min = 0;
            bound_max = 0.5;
        end
    else
        colored_metrics_char_plot = [colored_metrics_char '_'];
        if strcmp(colored_metrics_char, 'r2')
            bound_min = 10;
            bound_max = 25;
        elseif strcmp(colored_metrics_char, 't2')
            bound_min = 0.04;
            bound_max = 0.09;
        elseif strcmp(colored_metrics_char, 'r1')
            bound_min = 0;
            bound_max = 0.8;
        elseif strcmp(colored_metrics_char, 't1')
            bound_min = 1;
            bound_max = 2;
        elseif strcmp(colored_metrics_char, 'diso')
            bound_min = 0.1*1e-9;
            bound_max = 1.5*1e-9;
        elseif strcmp(colored_metrics_char, 'sqddelta')
            bound_min = 0.5;
            bound_max = 1;
        end
    end
else
    colored_metrics_char_plot = '';
end

%% Preparation
normalized_peak_map = zeros([nb_voxel_x, nb_voxel_y, nb_voxel_z N]);
color_peak_map = zeros([nb_voxel_x, nb_voxel_y, nb_voxel_z max_nb_peaks]);
norm_peaks_map = zeros([nb_voxel_x, nb_voxel_y, nb_voxel_z max_nb_peaks]);
for vx = 1:nb_voxel_x
    for vy = 1:nb_voxel_y
        for vz = 1:nb_voxel_z
            if mask(vx,vy,vz) && (sum(abs(peak_map(vx,vy,vz,:))) ~= 0)
                peaks = squeeze(peak_map(vx,vy,vz,:));
                
                if use_other_colored_peak_metric
                    colored_metric_voxel = map_cluster_metrics{vx,vy,vz};
                    colored_metric = zeros([1, max_nb_peaks]);
                    w = colored_metric_voxel.w;
                    if strcmp(colored_metrics_char, 'diso')
                        dpar = colored_metric_voxel.dpar;
                        dperp = colored_metric_voxel.dperp;
                        nb_peaks = length(dpar);
                        for i = 1:nb_peaks
                            diso = (dpar{i}+2*dperp{i})/3;
                            if do_IQR
                                colored_metric(i) = abs(diff(weighted_quantile(diso, w{i}, [0.25 0.75])));
                            else
                                colored_metric(i) = weighted_quantile(diso, w{i}, 0.5);
                            end
                        end          
                    elseif strcmp(colored_metrics_char, 'sqddelta')
                        dpar = colored_metric_voxel.dpar;
                        dperp = colored_metric_voxel.dperp;
                        nb_peaks = length(dpar);
                        for i = 1:nb_peaks
                            sqddelta = msf_notfinite2zero(((dpar{i} - dperp{i})/(dpar{i} + 2*dperp{i})).^2);
                            if do_IQR
                                colored_metric(i) = abs(diff(weighted_quantile(sqddelta, w{i}, [0.25 0.75])));
                            else
                                colored_metric(i) = weighted_quantile(sqddelta, w{i}, 0.5);
                            end
                        end
                    else
                        eval(strcat('cm = colored_metric_voxel.', colored_metrics_char, ';'))
                        nb_peaks = length(cm);
                        for i = 1:nb_peaks
                            if do_IQR
                                colored_metric(i) = abs(diff(weighted_quantile(cm{i}, w{i}, [0.25 0.75])));
                            else
                                colored_metric(i) = weighted_quantile(cm{i}, w{i}, 0.5);
                            end
                        end
                    end
                    colored_metric = (colored_metric - bound_min)/(bound_max - bound_min);
                    colored_metric(colored_metric > 1) = 1;
                    colored_metric(colored_metric < 0) = 0;
                end
                
                for m = 1:max_nb_peaks
                    if sum(abs(peaks(3*(m-1)+1:3*(m-1)+3))) ~= 0
                        norm_peak = norm(peaks(3*(m-1)+1:3*(m-1)+3));
                        norm_peaks_map(vx,vy,vz,m) = norm_peak;
                        if use_other_colored_peak_metric
                            rgb_colored_metric = dist_cind2rgb_jet(colored_metric(m));
                            color_peak_map(vx,vy,vz,3*(m-1)+1:3*(m-1)+3) = [rgb_colored_metric.r rgb_colored_metric.g rgb_colored_metric.b];
                        else
                            color_peak_map(vx,vy,vz,3*(m-1)+1:3*(m-1)+3) = abs(peaks(3*(m-1)+1:3*(m-1)+3))/norm_peak;
                        end
                    end
                    normalized_peak_map(vx,vy,vz,:) = fraction_map(vx,vy,vz)*peaks/max(norm_peaks_map(vx,vy,vz,:)); 
                end
            end
        end
    end
end
% max_norm = max(norm_peaks_map, [], 'all');
% peak_map = peak_map/max_norm;

peak_map = 1.05*normalized_peak_map/nanmax(fraction_map, [], 'all');

% Because of the strides
metric_map = flip(metric_map,2); 
metric_map = flip(metric_map,3);
peak_map = flip(peak_map,2);
peak_map = flip(peak_map,3);
color_peak_map = flip(color_peak_map,2);
color_peak_map = flip(color_peak_map,3);
peak_map(:,:,:,1:3:N) = -peak_map(:,:,:,1:3:N); 

% Because of the way Matlab reads the voxels
metric_map_axial = permute(metric_map, [2 1 3]);
peak_map_axial = permute(peak_map, [2 1 3 4]);
color_peak_map_axial = permute(color_peak_map, [2 1 3 4]);

metric_map_sagittal = permute(metric_map, [1 3 2]);
peak_map_sagittal = permute(peak_map, [1 3 2 4]);
color_peak_map_sagittal = permute(color_peak_map, [1 3 2 4]);

metric_map_coronal = permute(metric_map, [3 2 1]);
peak_map_coronal = permute(peak_map, [3 2 1 4]);
color_peak_map_coronal = permute(color_peak_map, [3 2 1 4]);

%% Axial figures
if do_axial
    output_fig_directory = fullfile(output_directory, 'axial'); %#ok<*UNRCH>
    if ~exist(output_fig_directory, 'dir')
        msf_mkdir(output_fig_directory);
    end
    
    scale = 1.35*nb_voxel_y/nb_voxel_x;
    width = 0.8;
    height = scale*width;
    x_left = (1-width)/2;
    y_bottom = (1-height)/2;
    
    parfor (vz = 1:nb_voxel_z, nb_parallel_workers)
        if sum(sum(mask(:,:,vz),1))
            
            metric = squeeze(metric_map_axial(:,:,vz));
            peaks = squeeze(peak_map_axial(:,:,vz,:));
            colors = squeeze(color_peak_map_axial(:,:,vz,:));
            
            f = figure('Position', [0,0,900,706], 'Units', 'pixels', 'visible', 'off');
            f.PaperOrientation = 'landscape';
            f.PaperUnits = 'normalized';
            set(f,'renderer','painters');
            f.PaperPosition = [0 0 1 1];
            
            axh = axes('position',[x_left y_bottom width height]);
            
            %% Maps
            colormap gray
            
            imagesc(axh, metric);
            
            for vx = 1:nb_voxel_x
                for vy = 1:nb_voxel_y
                    peaks_voxel = squeeze(peaks(vy,vx,:));
                    colors_voxel = squeeze(colors(vy,vx,:));
                    for m = 1:max_nb_peaks
                        peak = peaks_voxel(3*(m-1)+1:3*(m-1)+3);
                        color = colors_voxel(3*(m-1)+1:3*(m-1)+3);
                        if sum(color)
                            hold(axh, 'on')
                            quiver(axh, vx-0.5*peak(1), vy-0.5*peak(2), peak(1), peak(2), 'Color', color, 'LineWidth', 1.5, 'ShowArrowHead', 'off')
                        end
                    end
                end
            end
            
            %         originalSize_LS = get(axh_LS, 'Position');
            %         c = colorbar(axh_LS, 'Location', 'eastoutside', 'Fontsize', 16, 'TickLabelInterpreter', 'latex');
            %         set(axh_LS, 'Position', originalSize_LS)
            %         clims = get(c, 'Limits');
            %         cpos = get(c, 'Position');
            %         set(c, 'Ticks', [clims(1) clims(2)], 'TickLabels', {'0', num2str(round(max_median,2))}, 'Position', [cpos(1)+0.01 cpos(2) cpos(3) 2*cpos(4)+inter_v])
            %         text(axh_LS, 1.15, 1+0.5*inter_v/height, title_latex_colorbar, 'Fontsize', 16, 'Rotation', 270, 'Color', [0 0 0], 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'Units', 'Normalized', 'Interpreter', 'latex')
            
            %         text(axh, 0.96, 0.97, 'LL', 'Fontsize', 16, 'Color', [1 1 1], 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'Units', 'Normalized', 'Interpreter', 'latex')
            %         text(axh, 0.04, 0.97, '(a)', 'Fontsize', 16, 'Color', [1 1 1], 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'Units', 'Normalized', 'Interpreter', 'latex')
            
            set(axh,'xtick',[],'ytick',[],'xticklabel',[],'yticklabel',[])
            %         axh.XColor = [1 0 0];
            %         axh.YColor = [1 0 0];
            %         axh.LineWidth = 2.5;
            
            saveas(gcf, fullfile(output_fig_directory, strcat(colored_metrics_char_plot, 'peaks_axial_vz', num2str(vz), '.pdf')))
            clf(f)
            
        end
    end
end

clear vz

%% Coronal figures
if do_coronal
    output_fig_directory = fullfile(output_directory, 'coronal');
    if ~exist(output_fig_directory, 'dir')
        msf_mkdir(output_fig_directory);
    end
    
    scale = 1.35*nb_voxel_z/nb_voxel_x;
    width = 0.8;
    height = scale*width;
    x_left = (1-width)/2;
    y_bottom = (1-height)/2;
    
    parfor (vy = 1:nb_voxel_y, nb_parallel_workers)
        if sum(sum(mask(:,vy,:),1))
                      
            metric = squeeze(metric_map_coronal(:,vy,:));
            peaks = squeeze(peak_map_coronal(:,vy,:,:));
            colors = squeeze(color_peak_map_coronal(:,vy,:,:));
            
            f = figure('Position', [0,0,900,706], 'Units', 'pixels', 'visible', 'off');
            f.PaperOrientation = 'landscape';
            f.PaperUnits = 'normalized';
            set(f,'renderer','painters');
            f.PaperPosition = [0 0 1 1];
            
            axh = axes('position',[x_left y_bottom width height]);
            
            %% Maps
            colormap gray
            
            imagesc(axh, metric);
            
            for vx = 1:nb_voxel_x
                for vz = 1:nb_voxel_z
                    peaks_voxel = squeeze(peaks(vz,vx,:));
                    colors_voxel = squeeze(colors(vz,vx,:));
                    for m = 1:max_nb_peaks
                        peak = peaks_voxel(3*(m-1)+1:3*(m-1)+3);
                        color = colors_voxel(3*(m-1)+1:3*(m-1)+3);
                        if sum(color)
                            hold(axh, 'on')
                            quiver(axh, vx-0.5*peak(1), vz-0.5*peak(3), peak(1), peak(3), 'Color', color, 'LineWidth', 1.5, 'ShowArrowHead', 'off')
                        end
                    end
                end
            end
            
            %         originalSize_LS = get(axh_LS, 'Position');
            %         c = colorbar(axh_LS, 'Location', 'eastoutside', 'Fontsize', 16, 'TickLabelInterpreter', 'latex');
            %         set(axh_LS, 'Position', originalSize_LS)
            %         clims = get(c, 'Limits');
            %         cpos = get(c, 'Position');
            %         set(c, 'Ticks', [clims(1) clims(2)], 'TickLabels', {'0', num2str(round(max_median,2))}, 'Position', [cpos(1)+0.01 cpos(2) cpos(3) 2*cpos(4)+inter_v])
            %         text(axh_LS, 1.15, 1+0.5*inter_v/height, title_latex_colorbar, 'Fontsize', 16, 'Rotation', 270, 'Color', [0 0 0], 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'Units', 'Normalized', 'Interpreter', 'latex')
            
            %         text(axh, 0.96, 0.97, 'LL', 'Fontsize', 16, 'Color', [1 1 1], 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'Units', 'Normalized', 'Interpreter', 'latex')
            %         text(axh, 0.04, 0.97, '(a)', 'Fontsize', 16, 'Color', [1 1 1], 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'Units', 'Normalized', 'Interpreter', 'latex')
            
            set(axh,'xtick',[],'ytick',[],'xticklabel',[],'yticklabel',[])
            %         axh.XColor = [1 0 0];
            %         axh.YColor = [1 0 0];
            %         axh.LineWidth = 2.5;
            
            saveas(gcf, fullfile(output_fig_directory, strcat(colored_metrics_char_plot, 'peaks_coronal_vy', num2str(vy), '.pdf')))
            clf(f)
            
        end
    end
end

clear vy

%% Sagittal figures
if do_sagittal
    output_fig_directory = fullfile(output_directory, 'sagittal');
    if ~exist(output_fig_directory, 'dir')
        msf_mkdir(output_fig_directory);
    end
    
    scale = 1.35*nb_voxel_z/nb_voxel_y;
    width = 0.8;
    height = scale*width;
    x_left = (1-width)/2;
    y_bottom = (1-height)/2;
    
    parfor (vx = 1:nb_voxel_x, nb_parallel_workers)
        if sum(sum(mask(vx,:,:),1))
            
            metric = squeeze(metric_map_sagittal(vx,:,:));
            peaks = squeeze(peak_map_sagittal(vx,:,:,:));
            colors = squeeze(color_peak_map_sagittal(vx,:,:,:));
            
            f = figure('Position', [0,0,900,706], 'Units', 'pixels', 'visible', 'off');
            f.PaperOrientation = 'landscape';
            f.PaperUnits = 'normalized';
            set(f,'renderer','painters');
            f.PaperPosition = [0 0 1 1];
            
            axh = axes('position',[x_left y_bottom width height]);
            
            %% Maps
            colormap gray
            
            imagesc(axh, metric);
            
            for vy = 1:nb_voxel_y
                for vz = 1:nb_voxel_z
                    peaks_voxel = squeeze(peaks(vz,vy,:));
                    colors_voxel = squeeze(colors(vz,vy,:));
                    for m = 1:max_nb_peaks
                        peak = peaks_voxel(3*(m-1)+1:3*(m-1)+3);
                        color = colors_voxel(3*(m-1)+1:3*(m-1)+3);
                        if sum(color)
                            hold(axh, 'on')
                            quiver(axh, vy-0.5*peak(2), vz-0.5*peak(3), peak(2), peak(3), 'Color', color, 'LineWidth', 1.5, 'ShowArrowHead', 'off')
                        end
                    end
                end
            end
            
            
            %         originalSize_LS = get(axh_LS, 'Position');
            %         c = colorbar(axh_LS, 'Location', 'eastoutside', 'Fontsize', 16, 'TickLabelInterpreter', 'latex');
            %         set(axh_LS, 'Position', originalSize_LS)
            %         clims = get(c, 'Limits');
            %         cpos = get(c, 'Position');
            %         set(c, 'Ticks', [clims(1) clims(2)], 'TickLabels', {'0', num2str(round(max_median,2))}, 'Position', [cpos(1)+0.01 cpos(2) cpos(3) 2*cpos(4)+inter_v])
            %         text(axh_LS, 1.15, 1+0.5*inter_v/height, title_latex_colorbar, 'Fontsize', 16, 'Rotation', 270, 'Color', [0 0 0], 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'Units', 'Normalized', 'Interpreter', 'latex')
            
            %         text(axh, 0.96, 0.97, 'LL', 'Fontsize', 16, 'Color', [1 1 1], 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'Units', 'Normalized', 'Interpreter', 'latex')
            %         text(axh, 0.04, 0.97, '(a)', 'Fontsize', 16, 'Color', [1 1 1], 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'Units', 'Normalized', 'Interpreter', 'latex')
            
            set(axh,'xtick',[],'ytick',[],'xticklabel',[],'yticklabel',[])
            %         axh.XColor = [1 0 0];
            %         axh.YColor = [1 0 0];
            %         axh.LineWidth = 2.5;
            
            saveas(gcf, fullfile(output_fig_directory, strcat(colored_metrics_char_plot, 'peaks_sagittal_vx', num2str(vx), '.pdf')))
            clf(f)
            
        end
    end
end

clearvars