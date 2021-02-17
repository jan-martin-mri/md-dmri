function get_odf_metrics_single_voxel_figs()

close all

% vx, vy, vz = 38, 56, 2 - CC+cingulum

% CC
% vx = 40;
% vy = 44;
% vz = 14;

% CS
% vx = 48;
% vy = 39;
% vz = 16;

vx = 38;
vy = 56;
vz = 2;

% MDDMRI - tools - uvec - repulsion_angles_tree -> numbers that are allowed
nb_MeshNodes = 1000; % info needs to match odfs max. nodes odf 3994

method = 'dtr1d';
initial_directory = '';
directory_input = fullfile(pwd, initial_directory, '3_odfs');
directory_output = fullfile(pwd, initial_directory, 'ODF_peaks_single_voxel');

if ~exist(directory_output, 'dir')
    msf_mkdir(directory_output)
end

FLAG_relaxation = false;
if strcmp(method, 'dtr2d') || strcmp(method, 'dtr1d')
   FLAG_relaxation = true;     
end

max_n_peaks = 4;
threshold_w = 0.1;

disomin = 5e-10; disomax = 1.5e-9;
sqddeltamin = 0.25; sqddeltamax = 1;

ylim_diso = [0.5, 1.5];
ylim_sqddelta = [0.25, 1];

if strcmp(method, 'dtr2d') % t2 in ms
    rmin = 5; rmax = 20;
    tmin = 50; tmax = 150;
    ylim_r = [8, 15];
    ylim_t = [70, 110];
    char_r = 'r2';
    char_t = 't2';
    char_r_2 = 'R_2';
    char_t_2 = 'T_2';
elseif strcmp(method, 'dtr1d')
    rmin = 0; rmax = 1;
    tmin = 0; tmax = 4;
    ylim_r = [0 1];
    ylim_t = [0 4];
    char_r = 'r1';
    char_t = 't1';
    char_r_2 = 'R_1';
    char_t_2 = 'T_1';
end

% disomin = 1e-10; disomax = 1.5e-9;
% sqddeltamin = 0.25; sqddeltamax = 1;
% 
% ylim_diso = [0.75, 0.9];
% ylim_sqddelta = [0.85, 0.87];
% 
% if strcmp(method, 'dtr2d') % t2 in ms
%     rmin = 5; rmax = 20;
%     tmin = 50; tmax = 150;
%     ylim_r = [8, 15];
%     ylim_t = [70, 110];
%     char_r = 'r2';
%     char_t = 't2';
%     char_r_2 = 'R_2';
%     char_t_2 = 'T_2';
% elseif strcmp(method, 'dtr1d')
%     rmin = 0; rmax = 0.8;
%     tmin = 1; tmax = 2;
%     ylim_r = [8, 15];
%     ylim_t = [70, 110];
%     char_r = 'r1';
%     char_t = 't1';
%     char_r_2 = 'R_1';
%     char_t_2 = 'T_1';
% end
    
length_peaks = 1.22;

%% Load ODFs
odf_bsmedian_fn = fullfile(directory_input, ['odf_bsmedian_' num2str(nb_MeshNodes)]);
load(odf_bsmedian_fn, 'odf_bsmedian'); 
odf = odf_bsmedian;
odf.minimal_angular_separation = 7*pi/180;
odf_norm = 0.3*max(odf.w_bin{1}(:));

if FLAG_relaxation
    odf_s = struct;
    conn = struct;
    structure = struct;
    
    odf_s.n = odf.n;
    odf_s.x = odf.x;
    odf_s.y = odf.y;
    odf_s.z = odf.z;
    odf_s.c = abs([odf_s.x odf_s.y odf_s.z]);
    odf_s.tri = odf.tri;
    odf_s.w = squeeze(odf.w_bin{1}(vx,vy,vz,:))/odf_norm;
    odf_s.w_normalized = odf_s.w/sum(odf_s.w);
    odf_s.diso = squeeze(odf.diso_bin{1}(vx,vy,vz,:));
    odf_s.sqddelta = squeeze(odf.sqddelta_bin{1}(vx,vy,vz,:));
    if strcmp(method,'dtr2d')
        odf_s.r = squeeze(odf.r2_bin{1}(vx,vy,vz,:));
        odf_s.t = squeeze(odf.t2_bin{1}(vx,vy,vz,:));
    elseif strcmp(method,'dtr1d')
        odf_s.r = squeeze(odf.r1_bin{1}(vx,vy,vz,:));
        odf_s.t = squeeze(odf.t1_bin{1}(vx,vy,vz,:));
    end
    odf_s.points = [odf_s.x odf_s.y odf_s.z];
    odf_s.verts = repmat(odf_s.w,[1 3]).*[odf_s.x odf_s.y odf_s.z];
    odf_s.norm_verts = vecnorm(odf_s.verts')';
    % odf_s.norms = vertexNormal(triangulation(odf_s.tri,odf_s.verts),(1:odf_s.n)');
    
    % if nnz(isnan(odf_s.w)) > 0
    %    continue
    % end
    
    TR = triangulation(odf_s.tri, odf_s.verts);
    conn.basis = vertexAttachments(TR);
    indx = false(size(conn.basis,1),1);
    for i = 1:size(conn.basis,1)
        conn.tri = conn.basis{i};
        conn.verts = unique([odf_s.tri(conn.tri,1); odf_s.tri(conn.tri,2); odf_s.tri(conn.tri,3)]);
        if all(odf_s.norm_verts(i) >= odf_s.norm_verts(conn.verts))
            indx(i) = 1;
        end
    end
    
    % Filtering out low-probability peaks
    indw = odf_s.w_normalized/max(odf_s.w_normalized) >= threshold_w;
    
    odf_verts_peaks = odf_s.verts(indx & indw, :);
    diso_peaks = odf_s.diso(indx & indw, :);
    sqddelta_peaks = odf_s.sqddelta(indx & indw, :);
    r_peaks = odf_s.r(indx & indw, :);
    t_peaks = odf_s.t(indx & indw, :);
    w_peaks = odf_s.w_normalized(indx & indw, :);
    
    odf_verts_peaks_normalized = odf_verts_peaks;
    for i = 1:nnz(indx & indw)
        odf_verts_peaks_normalized(i,:) = odf_verts_peaks(i,:)/norm(odf_verts_peaks(i,:));
    end
    
    % Filtering out redundant antipodal points and taking up to n_peaks peaks with highest weights
    minimal_angular_separation = odf.minimal_angular_separation;
    relevant_ind = [];
    checked_ind = [];
    for i = 1:nnz(indx & indw)-1
        if ~ismember(i, checked_ind)
            for j = i+1:nnz(indx & indw)
                if abs(dot(odf_verts_peaks_normalized(i,:),odf_verts_peaks_normalized(j,:))) > cos(2*minimal_angular_separation)
                    if w_peaks(i) >= w_peaks(j)
                        relevant_ind = [relevant_ind i];
                    else
                        relevant_ind = [relevant_ind j];
                    end
                    checked_ind = [checked_ind j];
                end
            end
        end
    end
    
    odf_verts_peaks = odf_verts_peaks(relevant_ind, :);
    diso_peaks = diso_peaks(relevant_ind);
    sqddelta_peaks = sqddelta_peaks(relevant_ind);
    r_peaks = r_peaks(relevant_ind);
    t_peaks = t_peaks(relevant_ind);
    w_peaks = w_peaks(relevant_ind);
    
    norm_odf_verts_peaks = sqrt(odf_verts_peaks(:,1).^2+odf_verts_peaks(:,2).^2 + odf_verts_peaks(:,3).^2);
    [~, indx] = sort(norm_odf_verts_peaks,'descend');
    
    norm_peaks = double.empty;
    n_peaks = max_n_peaks;
    if numel(indx) < max_n_peaks
        odf_verts_peaks(1:numel(indx),:) = odf_verts_peaks(indx,:);
        norm_peaks(1:numel(indx),:) = norm_odf_verts_peaks(indx);
        diso_peaks = diso_peaks(indx);
        sqddelta_peaks = sqddelta_peaks(indx);
        r_peaks = r_peaks(indx);
        t_peaks = t_peaks(indx);
        w_peaks = w_peaks(indx);
        n_peaks = numel(indx);
    else
        odf_verts_peaks = odf_verts_peaks(indx(1:n_peaks),:);
        norm_peaks = norm_odf_verts_peaks(indx(1:n_peaks));
        diso_peaks = diso_peaks(indx(1:n_peaks));
        sqddelta_peaks = sqddelta_peaks(indx(1:n_peaks));
        r_peaks = r_peaks(indx(1:n_peaks));
        t_peaks = t_peaks(indx(1:n_peaks));
        w_peaks = w_peaks(indx(1:n_peaks));
    end
    
else
    odf_s = struct;
    conn = struct;
    structure = struct;
    
    odf_s.n = odf.n;
    odf_s.x = odf.x;
    odf_s.y = odf.y;
    odf_s.z = odf.z;
    odf_s.c = abs([odf_s.x odf_s.y odf_s.z]);
    odf_s.tri = odf.tri;
    odf_s.w = squeeze(odf.w_bin{1}(vx,vy,vz,:))/odf_norm;
    odf_s.w_normalized = odf_s.w/sum(odf_s.w);
    odf_s.diso = squeeze(odf.diso_bin{1}(vx,vy,vz,:));
    odf_s.sqddelta = squeeze(odf.sqddelta_bin{1}(vx,vy,vz,:));
    odf_s.points = [odf_s.x odf_s.y odf_s.z];
    odf_s.verts = repmat(odf_s.w,[1 3]).*[odf_s.x odf_s.y odf_s.z];
    odf_s.norm_verts = vecnorm(odf_s.verts')';
    % odf_s.norms = vertexNormal(triangulation(odf_s.tri,odf_s.verts),(1:odf_s.n)');
    
    % if nnz(isnan(odf_s.w)) > 0
    %    continue
    % end
    
    TR = triangulation(odf_s.tri, odf_s.verts);
    conn.basis = vertexAttachments(TR);
    indx = false(size(conn.basis,1),1);
    for i = 1:size(conn.basis,1)
        conn.tri = conn.basis{i};
        conn.verts = unique([odf_s.tri(conn.tri,1); odf_s.tri(conn.tri,2); odf_s.tri(conn.tri,3)]);
        if all(odf_s.norm_verts(i) >= odf_s.norm_verts(conn.verts))
            indx(i) = 1;
        end
    end
    
    % Filtering out low-probability peaks
    indw = odf_s.w_normalized/max(odf_s.w_normalized) >= threshold_w;
    odf_verts_peaks = odf_s.verts(indx & indw, :);
    diso_peaks = odf_s.diso(indx & indw, :);
    sqddelta_peaks = odf_s.sqddelta(indx & indw, :);
    w_peaks = odf_s.w_normalized(indx & indw, :);
    
    odf_verts_peaks_normalized = odf_verts_peaks;
    for i = 1:nnz(indx & indw)
        odf_verts_peaks_normalized(i,:) = odf_verts_peaks(i,:)/norm(odf_verts_peaks(i,:));
    end
    
    % Filtering out redundant antipodal points and taking up to n_peaks peaks with highest weights
    minimal_angular_separation = odf.minimal_angular_separation;
    relevant_ind = [];
    checked_ind = [];
    for i = 1:nnz(indx & indw)-1
        if ~ismember(i, checked_ind)
            for j = i+1:nnz(indx & indw)
                if abs(dot(odf_verts_peaks_normalized(i,:),odf_verts_peaks_normalized(j,:))) > cos(2*minimal_angular_separation)
                    if w_peaks(i) >= w_peaks(j)
                        relevant_ind = [relevant_ind i];
                    else
                        relevant_ind = [relevant_ind j];
                    end
                    checked_ind = [checked_ind j];
                end
            end
        end
    end
    
    odf_verts_peaks = odf_verts_peaks(relevant_ind, :);
    diso_peaks = diso_peaks(relevant_ind);
    sqddelta_peaks = sqddelta_peaks(relevant_ind);
    w_peaks = w_peaks(relevant_ind);
    
    norm_odf_verts_peaks = sqrt(odf_verts_peaks(:,1).^2+odf_verts_peaks(:,2).^2 + odf_verts_peaks(:,3).^2);
    [~, indx] = sort(norm_odf_verts_peaks,'descend');
    
    norm_peaks = double.empty;
    n_peaks = max_n_peaks;
    if numel(indx) < max_n_peaks
        odf_verts_peaks(1:numel(indx),:) = odf_verts_peaks(indx,:);
        norm_peaks(1:numel(indx),:) = norm_odf_verts_peaks(indx);
        diso_peaks = diso_peaks(indx);
        sqddelta_peaks = sqddelta_peaks(indx);
        w_peaks = w_peaks(indx);
        n_peaks = numel(indx);
    else
        odf_verts_peaks = odf_verts_peaks(indx(1:n_peaks),:);
        norm_peaks = norm_odf_verts_peaks(indx(1:n_peaks));
        diso_peaks = diso_peaks(indx(1:n_peaks));
        sqddelta_peaks = sqddelta_peaks(indx(1:n_peaks));
        w_peaks = w_peaks(indx(1:n_peaks));
    end
    
end

peaks = odf_verts_peaks;

color_peaks = peaks;
for i = 1:n_peaks
    color_peaks(i,:) = abs(color_peaks(i,:))/norm_peaks(i);
end

cind = (odf_s.diso-min(disomin))/(max(disomax)-min(disomin));
col = dist_cind2rgb_jet(cind);
odf_s.cdiso = [col.r col.g col.b];

cind = (odf_s.sqddelta-min(sqddeltamin))/(max(sqddeltamax)-min(sqddeltamin));
col = dist_cind2rgb_jet(cind);
odf_s.csqddelta = [col.r col.g col.b];

if FLAG_relaxation
    cind = (odf_s.r-min(rmin))/(max(rmax)-min(rmin));
    col = dist_cind2rgb_jet(cind);
    odf_s.cr = [col.r col.g col.b];
    
    cind = (odf_s.t-min(tmin))/(max(tmax)-min(tmin));
    col = dist_cind2rgb_jet(cind);
    odf_s.ct = [col.r col.g col.b];
end

peaks = length_peaks*peaks;

norm(peaks(1,:))

%% PLOT

set(0, 'defaultLegendInterpreter','latex');
lw = 2;
global_font_size = 16;
global_font_size_labels = 18;

width_ODF = 0.2;
height_ODF = 0.25;
width_peak_metrics = 0.15;
height_peak_metrics = 0.2;

inter_v_1 = 0.055;
inter_v_2 = 0.025;

offset_left = 0.037;
inter_h = 0.15;
inter_h_2 = 0.09;

y_bottom = (1 - inter_v_1 - inter_v_2 - 3*max(height_ODF, height_peak_metrics))/2;
y_top = y_bottom + inter_v_1 + inter_v_2 + 2*max(height_ODF, height_peak_metrics);
y1 = y_bottom;
y2 = y_bottom + inter_v_1 + max(height_ODF, height_peak_metrics);

x_top = (1 - width_ODF)/2;
x_left = offset_left + (1 - 2*width_ODF - 2*width_peak_metrics - inter_h - 2*inter_h_2)/2;
x_right = x_left + width_ODF + width_peak_metrics + inter_h;
x_left_2 = x_left + width_ODF + inter_h_2;
x_right_2 = x_right + width_ODF + inter_h_2;

azimuth = -50; %-50
elevation = 25; %25
    
f = figure('Position', [0,0,1500,1500], 'Units', 'pixels', 'visible', 'off');
f.PaperOrientation = 'landscape';
f.PaperUnits = 'normalized';
f.PaperPosition = [0 0 1 1];
set(gcf, 'Renderer', 'painters')

colormap jet

axh = axes('position',[x_top y_top width_ODF height_ODF]);
hold(axh, 'on')
axis tight, axis square, axis equal
view(azimuth,elevation)
p = patch('Faces',odf_s.tri,'Vertices',odf_s.verts);
set(p,'FaceColor','interp', 'FaceVertexCData', odf_s.c, 'EdgeColor', 'none', 'LineWidth',0.01)
tle = title(axh,'Local orientations', 'Interpreter', 'latex', 'Fontsize', global_font_size_labels);
set(tle, 'Units', 'Normalized', 'Position', [0.5, 0.85, 0]);
axis off

%% diso
axh = axes('position',[x_left y2 width_ODF height_ODF]);
hold(axh, 'on')
axis tight, axis square, axis equal
view(azimuth,elevation)
for i = 1:n_peaks
    plot3([peaks(i,1) -peaks(i,1)], [peaks(i,2) -peaks(i,2)], [peaks(i,3) -peaks(i,3)], 'LineStyle', '-', 'LineWidth', 2, 'Color', color_peaks(i,:))
end
p = patch('Faces',odf_s.tri,'Vertices',odf_s.verts);
set(p,'FaceColor','interp', 'FaceVertexCData', odf_s.cdiso, 'EdgeColor', 'none', 'LineWidth',0.01)

cb = colorbar(axh);
cb_lims = get(cb, 'Limits');
cb_pos = get(cb, 'Position');
set(cb, 'Ticks', cb_lims, 'TickLabels', {num2str(disomin*1e9), num2str(disomax*1e9)}, 'FontSize', global_font_size, 'TickLabelInterpreter', 'latex', 'Linewidth', 1)
set(cb, 'Position', [cb_pos(1) y2+(height_ODF-0.7*cb_pos(4))/2 cb_pos(3) 0.7*cb_pos(4)], 'Units', 'Normalized')
hold(axh, 'on')

axis off

%% sqddelta
axh = axes('position',[x_left y1 width_ODF height_ODF]);
hold(axh, 'on')
axis tight, axis square, axis equal
view(azimuth,elevation)
for i = 1:n_peaks
    plot3([peaks(i,1) -peaks(i,1)], [peaks(i,2) -peaks(i,2)], [peaks(i,3) -peaks(i,3)], 'LineStyle', '-', 'LineWidth', 2, 'Color', color_peaks(i,:))
end
p = patch('Faces',odf_s.tri,'Vertices',odf_s.verts);
set(p,'FaceColor','interp', 'FaceVertexCData', odf_s.csqddelta, 'EdgeColor', 'none', 'LineWidth',0.01)

cb = colorbar(axh);
cb_lims = get(cb, 'Limits');
cb_pos = get(cb, 'Position');
set(cb, 'Ticks', cb_lims, 'TickLabels', {num2str(sqddeltamin), num2str(sqddeltamax)}, 'FontSize', global_font_size, 'TickLabelInterpreter', 'latex', 'Linewidth', 1)
set(cb, 'Position', [cb_pos(1) y1+(height_ODF-0.7*cb_pos(4))/2 cb_pos(3) 0.7*cb_pos(4)], 'Units', 'Normalized')
hold(axh, 'on')

axis off

if FLAG_relaxation
    %% Relaxation
    axh = axes('position',[x_right y2 width_ODF height_ODF]);
    hold(axh, 'on')
    axis tight, axis square, axis equal
    view(azimuth,elevation)
    for i = 1:n_peaks
        plot3([peaks(i,1) -peaks(i,1)], [peaks(i,2) -peaks(i,2)], [peaks(i,3) -peaks(i,3)], 'LineStyle', '-', 'LineWidth', 2, 'Color', color_peaks(i,:))
    end
    p = patch('Faces',odf_s.tri,'Vertices',odf_s.verts);
    set(p,'FaceColor','interp', 'FaceVertexCData', odf_s.cr, 'EdgeColor', 'none', 'LineWidth',0.01)
    
    cb = colorbar(axh);
    cb_lims = get(cb, 'Limits');
    cb_pos = get(cb, 'Position');
    set(cb, 'Ticks', cb_lims, 'TickLabels', {num2str(rmin), num2str(rmax)}, 'FontSize', global_font_size, 'TickLabelInterpreter', 'latex', 'Linewidth', 1)
    set(cb, 'Position', [cb_pos(1) y2+(height_ODF-0.7*cb_pos(4))/2 cb_pos(3) 0.7*cb_pos(4)], 'Units', 'Normalized')
    hold(axh, 'on')
    
    axis off
    
    axh = axes('position',[x_right y1 width_ODF height_ODF]);
    hold(axh, 'on')
    axis tight, axis square, axis equal
    view(azimuth,elevation)
    for i = 1:n_peaks
        plot3([peaks(i,1) -peaks(i,1)], [peaks(i,2) -peaks(i,2)], [peaks(i,3) -peaks(i,3)], 'LineStyle', '-', 'LineWidth', 2, 'Color', color_peaks(i,:))
    end
    p = patch('Faces',odf_s.tri,'Vertices',odf_s.verts);
    set(p,'FaceColor','interp', 'FaceVertexCData', odf_s.ct, 'EdgeColor', 'none', 'LineWidth',0.01)
    
    cb = colorbar(axh);
    cb_lims = get(cb, 'Limits');
    cb_pos = get(cb, 'Position');
    set(cb, 'Ticks', cb_lims, 'TickLabels', {num2str(tmin), num2str(tmax)}, 'FontSize', global_font_size, 'TickLabelInterpreter', 'latex', 'Linewidth', 1)
    set(cb, 'Position', [cb_pos(1) y1+(height_ODF-0.7*cb_pos(4))/2 cb_pos(3) 0.7*cb_pos(4)], 'Units', 'Normalized')
    hold(axh, 'on')
    
    axis off
end

hold(axh, 'on')

if FLAG_relaxation
    param_names = {'diso'; 'sqddelta'; 'r'; 't'};
    if strcmp(method, 'dtr2d')
        param_names_plot = {'$\hat{\mathrm{E}}[D_\mathrm{iso}]\,(\mu\mathrm{m}^2/\mathrm{ms})$'; '$\hat{\mathrm{E}}[D_\Delta^2]$'; ['$\hat{\mathrm{E}}[' char_r_2 ']\,(\mathrm{s}^{-1})$']; ['$\hat{\mathrm{E}}[' char_t_2 ']\,(\mathrm{ms})$']};
    elseif strcmp(method, 'dtr1d')
        param_names_plot = {'$\hat{\mathrm{E}}[D_\mathrm{iso}]\,(\mu\mathrm{m}^2/\mathrm{ms})$'; '$\hat{\mathrm{E}}[D_\Delta^2]$'; ['$\hat{\mathrm{E}}[' char_r_2 ']\,(\mathrm{s}^{-1})$']; ['$\hat{\mathrm{E}}[' char_t_2 ']\,(\mathrm{s})$']};
    end
else
    param_names = {'diso'; 'sqddelta'};
    param_names_plot = {'$\hat{\mathrm{E}}[D_\mathrm{iso}]\,(\mu\mathrm{m}^2/\mathrm{ms})$'; '$\hat{\mathrm{E}}[D_\Delta^2]$'};    
end
Nparam = numel(param_names);

axh_diso = axes('position',[x_left_2 y2+(height_ODF - height_peak_metrics)/2 width_peak_metrics height_peak_metrics],'YAxisLocation','left');
axh_sqddelta = axes('position',[x_left_2 y1+(height_ODF - height_peak_metrics)/2 width_peak_metrics height_peak_metrics],'YAxisLocation','left');
ylim(axh_diso, ylim_diso)
ylim(axh_sqddelta, ylim_sqddelta)

if FLAG_relaxation
    axh_r = axes('position',[x_right_2 y2+(height_ODF - height_peak_metrics)/2 width_peak_metrics height_peak_metrics],'YAxisLocation','left');
    axh_t = axes('position',[x_right_2 y1+(height_ODF - height_peak_metrics)/2 width_peak_metrics height_peak_metrics],'YAxisLocation','left');
    ylim(axh_r, ylim_r)
    ylim(axh_t, ylim_t)
end

for nparam = 1:Nparam
    eval(['axh = axh_' param_names{nparam} ';'])
    hold(axh,'on')
    
    axh.XAxis.TickLabelInterpreter = 'latex';
    axh.YAxis.TickLabelInterpreter = 'latex';
    
    eval(['quantity = ' param_names{nparam} '_peaks;'])
    if strcmp(param_names{nparam},'diso')
        quantity = quantity*1e9;
    end
    
    for i = 1:n_peaks
        color = color_peaks(i,:);
        scatter(axh, i, quantity(i), 60, 'filled', 's', 'MarkerFaceColor', color, 'MarkerEdgeColor', 'k')
    end
    
    box(axh,'off')
    grid(axh,'on')
    set(axh,'XTick',0:n_peaks+1)
    set(axh,'XTickLabel',{' '})
    
    set(axh, 'Fontsize', global_font_size)
    xlim(axh, [0, n_peaks+1])
    
    %     tle = title(axh, param_names_plot{nparam}, 'Fontsize', global_font_size, 'Interpreter', 'latex');
    %     set(tle, 'Units', 'Normalized', 'Position', [0.5, 1.15, 0]);
    
    %     yyaxis right
    ylab = ylabel(axh, param_names_plot{nparam}, 'Fontsize', global_font_size_labels, 'Interpreter', 'latex');
    set(ylab, 'Units', 'Normalized', 'Position', [-1.75, 0.5, 0]);
    
%     if nparam == 2 || nparam == Nparam
    xlab = xlabel(axh, 'Peaks', 'Fontsize', global_font_size_labels, 'Interpreter', 'latex');
    set(xlab, 'Units', 'Normalized', 'Position', [0.5, -0.1, 0]);
%     end
    
    set(axh, 'LineWidth', lw)
    set(axh, 'TickDir','out');
end

saveas(gcf,fullfile(directory_output, strcat('ODF_metrics_nbMeshNodes', num2str(nb_MeshNodes), '_vx', num2str(vx), '_vy', num2str(vy), '_vz', num2str(vz), '.pdf')))
clf(f)


end