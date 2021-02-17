function get_odf_metrics_ROI()

close all

method = 'dtr1d';
initial_directory = 'Data_T1_D';
maps_directory = fullfile(pwd, initial_directory, '2_nii_pdf_maps');
output_directory = fullfile(pwd, initial_directory, 'ODF_metrics_ROI');

% input_directory = fullfile(pwd, initial_directory, '3_odfs');
% load(fullfile(input_directory, 'ODF_peaks_metrics.mat'), 'cell_odf_peaks')
% peak_metrics = cell_odf_peaks;

input_directory = fullfile(pwd, initial_directory, '4_clustering');
load(fullfile(input_directory, 'cluster_metrics.mat'), 'map_cluster_metrics')
peak_metrics = map_cluster_metrics;

ROI_name = 'three-way';
% ROI_name = 'two-way';

if strcmp(ROI_name, 'three-way')
    ylim_diso = [0.75 1.2];
    ylim_sddelta = [0.92 1];
    ylim_r = [0.5 0.8];
    ylim_t = [1.3 2.1];
elseif strcmp(ROI_name, 'two-way')
    ylim_diso = [0.75 1.3];
    ylim_sddelta = [0.7 1];
    ylim_r = [0.35 0.85];
    ylim_t = [1 3.5];
end

[mask, h] = mdm_nii_read(fullfile(maps_directory, [method '_fraction_thin.nii.gz']));
mask = mask > 0;
[nb_vx, nb_vy, nb_vz] = size(mask);

if ~exist(output_directory, 'dir')
    msf_mkdir(output_directory)
end

FLAG_relaxation = false;
if strcmp(method, 'dtr2d') || strcmp(method, 'dtr1d')
   FLAG_relaxation = true;     
end

if strcmp(method, 'dtr2d')
    char_r = 'r2';
    char_t = 't2';
    char_r_2 = 'R_2';
    char_t_2 = 'T_2';
elseif strcmp(method, 'dtr1d')
    char_r = 'r1';
    char_t = 't1';
    char_r_2 = 'R_1';
    char_t_2 = 'T_1';
end

% if FLAG_relaxation
%     param_names = {'diso'; 'sddelta'; 'r'; 't'};
%     if strcmp(method, 'dtr2d')
%         param_names_plot = {'$\hat{\mathrm{E}}[D_\mathrm{iso}]\,(\mu\mathrm{m}^2/\mathrm{ms})$'; '$\hat{\mathrm{E}}[D_\Delta^2]$'; ['$\hat{\mathrm{E}}[' char_r_2 ']\,(\mathrm{s}^{-1})$']; ['$\hat{\mathrm{E}}[' char_t_2 ']\,(\mathrm{ms})$']};
%     elseif strcmp(method, 'dtr1d')
%         param_names_plot = {'$\hat{\mathrm{E}}[D_\mathrm{iso}]\,(\mu\mathrm{m}^2/\mathrm{ms})$'; '$\hat{\mathrm{E}}[D_\Delta^2]$'; ['$\hat{\mathrm{E}}[' char_r_2 ']\,(\mathrm{s}^{-1})$']; ['$\hat{\mathrm{E}}[' char_t_2 ']\,(\mathrm{s})$']};
%     end
% else
%     param_names = {'diso'; 'sddelta'};
%     param_names_plot = {'$\hat{\mathrm{E}}[D_\mathrm{iso}]\,(\mu\mathrm{m}^2/\mathrm{ms})$'; '$\hat{\mathrm{E}}[D_\Delta^2]$'};    
% end

if FLAG_relaxation
    param_names = {'diso'; 'sddelta'; 'r'; 't'};
    if strcmp(method, 'dtr2d')
        param_names_plot = {'$\mathring{\mathrm{E}}[D_\mathrm{iso}]\,(\mu\mathrm{m}^2/\mathrm{ms})$'; '$\mathring{\mathrm{E}}[D_\Delta^2]$'; ['$\mathring{\mathrm{E}}[' char_r_2 ']\,(\mathrm{s}^{-1})$']; ['$\mathring{\mathrm{E}}[' char_t_2 ']\,(\mathrm{ms})$']};
    elseif strcmp(method, 'dtr1d')
        param_names_plot = {'$\mathring{\mathrm{E}}[D_\mathrm{iso}]\,(\mu\mathrm{m}^2/\mathrm{ms})$'; '$\mathring{\mathrm{E}}[D_\Delta^2]$'; ['$\mathring{\mathrm{E}}[' char_r_2 ']\,(\mathrm{s}^{-1})$']; ['$\mathring{\mathrm{E}}[' char_t_2 ']\,(\mathrm{s})$']};
    end
else
    param_names = {'diso'; 'sddelta'};
    param_names_plot = {'$\mathring{\mathrm{E}}[D_\mathrm{iso}]\,(\mu\mathrm{m}^2/\mathrm{ms})$'; '$\mathring{\mathrm{E}}[D_\Delta^2]$'};    
end

Nparam = numel(param_names);

%% Define ROI

if strcmp(ROI_name, 'three-way')
    ROI_vx = 28:36;
    mirror_ROI_vx = nb_vx - ROI_vx + 1;
    ROI_vx = [ROI_vx mirror_ROI_vx];
    ROI_vy = 33:41;
    ROI_vz = 3:4;
elseif strcmp(ROI_name, 'two-way')
    ROI_vx = 38:39;
    mirror_ROI_vx = nb_vx - ROI_vx + 2;
    ROI_vx = [ROI_vx mirror_ROI_vx];
    ROI_vy = 50:55;
    mirror_ROI_vy = nb_vy - ROI_vy + 7;
    ROI_vy = [ROI_vy mirror_ROI_vy];
    ROI_vz = 3:4;
end

mask_ROI = zeros([nb_vx nb_vy nb_vz]);
[ROI_vx, ROI_vy, ROI_vz] = meshgrid(ROI_vx, ROI_vy, ROI_vz);
for vx = ROI_vx
    for vy = ROI_vy
        for vz = ROI_vz
            mask_ROI(vx,vy,vz) = 1;
        end
    end
end
mask_ROI = mask_ROI.*mask;
mdm_nii_write(mask_ROI, fullfile(output_directory, ['mask_ROI_' ROI_name '.nii.gz']), h);

%%
if strcmp(ROI_name, 'three-way')
    char_fibers = {'CC'; 'AF'; 'CST'};
    orientations_fibers = [[1 0 0]; [0 1 0]; [0 0 1]];
    color_group = {[1 0 0]; [0 0.8 0]; [0 0 1]};
elseif strcmp(ROI_name, 'two-way')
    char_fibers = {'CC'; 'CING'};
    orientations_fibers = [[1 0 0]; [0 1 0]];
    color_group = {[1 0 0]; [0 0.8 0]};
end

for p = 1:length(char_fibers)
   eval(['w_' char_fibers{p} ' = [];']); 
   eval(['diso_' char_fibers{p} ' = [];']);
   eval(['sddelta_' char_fibers{p} ' = [];']);
   if FLAG_relaxation
       eval(['r_' char_fibers{p} ' = [];']);
       eval(['t_' char_fibers{p} ' = [];']);
   end
end

for vx = 1:nb_vx
    for vy = 1:nb_vy
        for vz = 1:nb_vz
            metrics = peak_metrics{vx,vy,vz};
            if mask_ROI(vx, vy, vz) && ~isempty(metrics)
                
%                 n = metrics.n_peaks;
%                 orientations = metrics.orientation_peaks;
%                 w = metrics.w_peaks;
%                 diso = metrics.diso_peaks;
%                 sddelta = metrics.sddelta_peaks; 
%                 if strcmp(method, 'dtr1d')
%                     r = metrics.r1_peaks;
%                     t = metrics.t1_peaks;
%                 elseif strcmp(method, 'dtr1d')
%                     r = metrics.r2_peaks;
%                     t = metrics.t2_peaks;
%                 end
                
                n = metrics.n;
                orientations = metrics.median_orientations;
                w = zeros([1 n]);
                diso = zeros([1 n]);
                sddelta = zeros([1 n]);
                r = zeros([1 n]);
                t = zeros([1 n]);
                
                for m = 1:n
                    ind = (3*(m-1)+1):(3*(m-1)+3);
                    orientation = orientations(ind);
                    
                    relevant_dots = zeros([1 length(char_fibers)]);
                    for p = 1:length(char_fibers)
                        relevant_dots(p) = abs(dot(orientation, orientations_fibers(p,:)));
                    end
                    
                    [~, ind_max_dot] = max(relevant_dots);
                    char_fiber = char_fibers{ind_max_dot};
                    
                    w(m) = median(metrics.w{m});
                    diso(m) = weighted_median((metrics.dpar{m} + 2*metrics.dperp{m})./3, metrics.w{m});
                    sddelta(m) = weighted_median(((metrics.dpar{m} - metrics.dperp{m})./(metrics.dpar{m} + 2*metrics.dperp{m})).^2, metrics.w{m});
                    r(m) = weighted_median(metrics.r1{m}, metrics.w{m});
                    t(m) = weighted_median(metrics.t1{m}, metrics.w{m});
                    
                    eval(['w_' char_fiber ' = [w_' char_fiber ' w(m)];'])
                    eval(['diso_' char_fiber ' = [diso_' char_fiber ' diso(m)];'])
                    eval(['sddelta_' char_fiber ' = [sddelta_' char_fiber ' sddelta(m)];'])
                    if FLAG_relaxation
                        eval(['r_' char_fiber ' = [r_' char_fiber ' r(m)];'])
                        eval(['t_' char_fiber ' = [t_' char_fiber ' t(m)];'])
                    end
                end
            end
        end
    end
end

if strcmp(ROI_name, 'three-way')
    [ranksum(diso_CC, diso_AF) ranksum(diso_CC, diso_CST) ranksum(diso_AF, diso_CST)]
    [ranksum(sddelta_CC, sddelta_AF) ranksum(sddelta_CC, sddelta_CST) ranksum(sddelta_AF, sddelta_CST)]
    [ranksum(r_CC, r_AF) ranksum(r_CC, r_CST) ranksum(r_AF, r_CST)]
    [ranksum(t_CC, t_AF) ranksum(t_CC, t_CST) ranksum(t_AF, t_CST)]
elseif strcmp(ROI_name, 'two-way')
    ranksum(diso_CC, diso_CING)
    ranksum(sddelta_CC, sddelta_CING)
    ranksum(r_CC, r_CING)
    ranksum(t_CC, t_CING)
end

%% PLOT

set(0, 'defaultLegendInterpreter','latex');
lw = 2;
global_font_size = 16;
global_font_size_labels = 18;

width = 0.3;
height = 0.3;

offset_left = 0.037;
inter_v = 0.15;
inter_h = 0.15;

box_width = 0.4;

y_bottom = (1 - inter_v - 2*height)/2;
y_top = y_bottom + inter_v + height;

x_left = offset_left + (1 - inter_h - 2*width)/2;
x_right = x_left + width + inter_h;

f = figure('Position', [0,0,900,706], 'Units', 'pixels', 'visible', 'off');
f.PaperOrientation = 'landscape';
f.PaperUnits = 'normalized';
f.PaperPosition = [0 0 1 1];

axh_diso = axes('position',[x_left y_top width height]);
axh_sddelta = axes('position',[x_left y_bottom width height]);

if FLAG_relaxation
    axh_r = axes('position',[x_right y_top width height],'YAxisLocation','left');
    axh_t = axes('position',[x_right y_bottom width height],'YAxisLocation','left');
end

for nparam = 1:Nparam
    eval(['axh = axh_' param_names{nparam} ';'])
    hold(axh,'on')
    
    for m = 1:length(char_fibers) 
        eval(['quantity = ' param_names{nparam} '_' char_fibers{m} ';'])
        if strcmp(param_names{nparam},'diso')
            quantity = quantity*1e9;
        elseif strcmp(param_names{nparam},'t') && strcmp(method, 'dtr2d')
            quantity = quantity*1e3;
        end
        
        boxplot(axh, quantity, 'Color', color_group{m}, 'Whisker', 0, 'boxstyle', 'filled', 'Symbol', '', 'OutlierSize', 3, 'widths', box_width, 'position', m)     
    end
    
    axh.XAxis.TickLabelInterpreter = 'latex';
    axh.YAxis.TickLabelInterpreter = 'latex';
    
    set(findobj(axh,'tag','Median'),'linewidth',2);
    set(findobj(axh,'tag','Whisker'),'linewidth',1.5);
    set(findobj(axh,'tag','Box'),'linewidth',4.5);
    
    box(axh,'off')
    grid(axh,'off')
    set(axh,'XTick',0:length(char_fibers)+1)
    
    XTickCell{1} = '';
    for p = 1:length(char_fibers)
        XTickCell{p+1} = char_fibers{p};
    end
    XTickCell{length(char_fibers)+2} = '';
    set(axh,'XTickLabel', XTickCell)
    
    set(axh, 'Fontsize', global_font_size)
    xlim(axh, [0, length(char_fibers)+1])
    
    ylab = ylabel(axh, param_names_plot{nparam}, 'Fontsize', global_font_size_labels, 'Interpreter', 'latex');
    set(ylab, 'Units', 'Normalized', 'Position', [-0.2, 0.5, 0]);
        
    eval(['ylim(axh_' param_names{nparam} ', ylim_' param_names{nparam} ')'])
    
    set(axh, 'LineWidth', lw)
    set(axh, 'TickDir','out');
    
    if strcmp(ROI_name, 'three-way')
        if strcmp(param_names{nparam}, 'diso')
            plot_statistical_significance_line(axh, [1 2], 1.02+0.1, 2)
            plot_statistical_significance_line(axh, [1 3], 1.06+0.1, 1)
        elseif strcmp(param_names{nparam}, 'sddelta')
            plot_statistical_significance_line(axh, [1 3], 0.986, 1)
            plot_statistical_significance_line(axh, [2 3], 0.978, 3)
        elseif strcmp(param_names{nparam}, 'r')
            plot_statistical_significance_line(axh, [1 3], 0.79, 3)
            plot_statistical_significance_line(axh, [2 3], 0.765, 3)
        elseif strcmp(param_names{nparam}, 't')
            plot_statistical_significance_line(axh, [1 3], 1.85+0.21, 3)
            plot_statistical_significance_line(axh, [2 3], 1.78+0.21, 3)
        end
    elseif strcmp(ROI_name, 'two-way')
        if strcmp(param_names{nparam}, 'diso')
            plot_statistical_significance_line(axh, [1 2], 1.25, 1)
        elseif strcmp(param_names{nparam}, 'sddelta')
            plot_statistical_significance_line(axh, [1 2], 0.98, 1)
        elseif strcmp(param_names{nparam}, 'r')
            plot_statistical_significance_line(axh, [1 2], 0.78, 2)
        end
    end
end

saveas(gcf,fullfile(output_directory, strcat('ODF_metrics_ROI_', ROI_name, '.pdf')))
clf(f)


end