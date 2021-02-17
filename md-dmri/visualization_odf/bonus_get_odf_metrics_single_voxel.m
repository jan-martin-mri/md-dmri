function get_odf_metrics_single_voxel(vx, vy, vz, nb_MeshNodes, directory)

close all

% vx, vy, vz = 38, 56, 2 - CC+cingulum

n_peaks = 3;
threshold_w = 0.25;

disomin = 1e-10; disomax = 3.5e-9;
sqddeltamin = 0; sqddeltamax = 1;
r2min = 0; r2max = 30;
t2min = 0.03; t2max = 0.12;

%% Load ODFs
odf_bsmedian_fn = fullfile(directory, ['odf_bsmedian_' num2str(nb_MeshNodes)]);
load(odf_bsmedian_fn, 'odf_bsmedian'); 
odf = odf_bsmedian;
odf_norm = 0.3*max(odf.w_bin{1}(:));

odf_s.n = odf.n;
odf_s.x = odf.x;
odf_s.y = odf.y;
odf_s.z = odf.z;
odf_s.c = abs([odf_s.x odf_s.y odf_s.z]);
odf_s.tri = odf.tri;

% [nb_voxel_x, nb_voxel_y, nb_voxel_z, ~] = size(odf.w_bin{1});

% odf_s.w = squeeze(odf.w_bin{1}(vx,vy,vz,:));
odf_s.w = squeeze(odf.w_bin{1}(vx,vy,vz,:))/odf_norm;
odf_s.w_normalized = odf_s.w/sum(odf_s.w);
odf_s.diso = squeeze(odf.diso_bin{1}(vx,vy,vz,:));
odf_s.sqddelta = squeeze(odf.sqddelta_bin{1}(vx,vy,vz,:));
odf_s.r2 = squeeze(odf.r2_bin{1}(vx,vy,vz,:));
odf_s.t2 = squeeze(odf.t2_bin{1}(vx,vy,vz,:));
odf_s.points = [odf_s.x odf_s.y odf_s.z];
odf_s.verts = repmat(odf_s.w,[1 3]).*[odf_s.x odf_s.y odf_s.z];
odf_s.norm_verts = vecnorm(odf_s.verts')';
odf_s.norms = vertexNormal(triangulation(odf_s.tri,odf_s.verts),(1:odf_s.n)');

cind = (odf_s.diso-disomin)/(disomax-disomin);
col = dist_cind2rgb_jet(cind);
odf_s.cdiso = [col.r col.g col.b];
cind = (odf_s.sqddelta-sqddeltamin)/(sqddeltamax-sqddeltamin);
col = dist_cind2rgb_jet(cind);
odf_s.csqddelta = [col.r col.g col.b];
cind = (odf_s.r2-r2min)/(r2max-r2min);
col = dist_cind2rgb_jet(cind);
odf_s.cr2 = [col.r col.g col.b];
cind = (odf_s.t2-t2min)/(t2max-t2min);
col = dist_cind2rgb_jet(cind);
odf_s.ct2 = [col.r col.g col.b];

TR = triangulation(odf_s.tri, odf_s.verts);
conn.basis = vertexAttachments(TR);
% peaks = zeros(n_peaks,3);
indx = false(size(conn.basis,1),1);
for i = 1:size(conn.basis,1)
    conn.tri = conn.basis{i};
    conn.verts = unique([odf_s.tri(conn.tri,1); odf_s.tri(conn.tri,2); odf_s.tri(conn.tri,3)]);
    if all(odf_s.norm_verts(i) >= odf_s.norm_verts(conn.verts))
        indx(i) = 1;
    end
end

% Filtering out low-probability peaks
indw = odf_s.w_normalized >= threshold_w*max(odf_s.w_normalized);
odf_verts_peaks = odf_s.verts(indx & indw, :);
diso_peaks = odf_s.diso(indx & indw, :);
sqddelta_peaks = odf_s.sqddelta(indx & indw, :);
r2_peaks = odf_s.r2(indx & indw, :);
t2_peaks = odf_s.t2(indx & indw, :);

% Filtering out redundant antipodal points (keep z > 0) and taking up to n_peaks peaks with highest weights
ind_z_positive = odf_verts_peaks(:,3) > 0;
odf_verts_peaks = odf_verts_peaks(ind_z_positive, :);
diso_peaks = diso_peaks(ind_z_positive);
sqddelta_peaks = sqddelta_peaks(ind_z_positive);
r2_peaks = r2_peaks(ind_z_positive);
t2_peaks = t2_peaks(ind_z_positive);

norm_odf_verts_peaks = sqrt(odf_verts_peaks(:,1).^2+odf_verts_peaks(:,2).^2 + odf_verts_peaks(:,3).^2);
[~, indx] = sort(norm_odf_verts_peaks,'descend');
if numel(indx) < n_peaks
    peaks(1:numel(indx),:) = odf_verts_peaks(indx,:);
    norm_peaks(1:numel(indx),:) = norm_odf_verts_peaks(indx);
    diso_peaks = diso_peaks(indx);
    sqddelta_peaks = sqddelta_peaks(indx);
    r2_peaks = r2_peaks(indx);
    t2_peaks = t2_peaks(indx);
    n_peaks = numel(indx);
else
    peaks = odf_verts_peaks(indx(1:n_peaks),:);
    norm_peaks = norm_odf_verts_peaks(indx(1:n_peaks));
    diso_peaks = diso_peaks(indx(1:n_peaks));
    sqddelta_peaks = sqddelta_peaks(indx(1:n_peaks));
    r2_peaks = r2_peaks(indx(1:n_peaks));
    t2_peaks = t2_peaks(indx(1:n_peaks));
end

color_peaks = peaks;
for i = 1:n_peaks
    color_peaks(i,:) = abs(color_peaks(i,:))/norm_peaks(i);
end


% PLOT

set(0, 'defaultLegendInterpreter','latex');
lw = 2;
global_font_size = 16;
global_font_size_labels = 18;


width_left = 0.2;
height_left = 0.3;

width_middle = width_left;
height_middle = height_left;

width_right = 0.15;
height_right = 0.2;

inter_h_LM = 0.1;
inter_h_MR = 0.01;

inter_v = 0.10;
inter_v_middle = 0.01;

x_left = (1-inter_h_LM-inter_h_MR-width_left-width_middle-width_right)/2-0.05;
x_middle = x_left + width_left + inter_h_LM;  
x_right = x_middle + width_middle + inter_h_MR;

y_left = (1-height_left)/2;
y1_middle = (1-3*height_middle-2*inter_v_middle)/2;
y2_middle = y1_middle + inter_v_middle + height_middle;
y3_middle = y2_middle + inter_v_middle + height_middle;
y1_right = (1-3*height_right-2*inter_v)/2;
y2_right = y1_right + inter_v + height_right;
y3_right = y2_right + inter_v + height_right;

for k = 1:2
    f = figure('Position', [0,0,1500,1500], 'Units', 'pixels', 'visible', 'off');
    f.PaperOrientation = 'landscape';
    f.PaperUnits = 'normalized';
    f.PaperPosition = [0 0 1 1];
    set(gcf, 'Renderer', 'painters')
    
    axh = axes('position',[x_left y_left width_left height_left]);
    hold(axh, 'on')
    axis tight, axis square, axis equal
    view(-50,25)
    p = patch('Faces',odf_s.tri,'Vertices',odf_s.verts);
    set(p,'FaceColor','interp', 'FaceVertexCData', odf_s.c, 'EdgeColor', 'none', 'LineWidth',0.01)
    % for i = 1:n_peaks
    %     scatter3(peaks(i,1), peaks(i,2), peaks(i,3), 100, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', color_peaks(i,:))
    % end
    tle = title(axh,'Local orientations', 'Interpreter', 'latex', 'Fontsize', global_font_size_labels);
    set(tle, 'Units', 'Normalized', 'Position', [0.5, 0.85, 0]);
    axis off
    
    axh = axes('position',[x_middle y3_middle width_middle height_middle]);
    hold(axh, 'on')
    axis tight, axis square, axis equal
    view(-50,25)
    p = patch('Faces',odf_s.tri,'Vertices',odf_s.verts);
    set(p,'FaceColor','interp', 'FaceVertexCData', odf_s.cdiso, 'EdgeColor', 'none', 'LineWidth',0.01)
    for i = 1:n_peaks
        scatter3(peaks(i,1), peaks(i,2), peaks(i,3), 100, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', color_peaks(i,:))
    end
    axis off
    
    axh = axes('position',[x_middle y2_middle width_middle height_middle]);
    hold(axh, 'on')
    axis tight, axis square, axis equal
    view(-50,25)
    p = patch('Faces',odf_s.tri,'Vertices',odf_s.verts);
    set(p,'FaceColor','interp', 'FaceVertexCData', odf_s.csqddelta, 'EdgeColor', 'none', 'LineWidth',0.01)
    for i = 1:n_peaks
        scatter3(peaks(i,1), peaks(i,2), peaks(i,3), 100, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', color_peaks(i,:))
    end
    axis off
    
    if k == 1
        axh = axes('position',[x_middle y1_middle width_middle height_middle]);
        hold(axh, 'on')
        axis tight, axis square, axis equal
        view(-50,25)
        p = patch('Faces',odf_s.tri,'Vertices',odf_s.verts);
        set(p,'FaceColor','interp', 'FaceVertexCData', odf_s.cr2, 'EdgeColor', 'none', 'LineWidth',0.01)
        for i = 1:n_peaks
            scatter3(peaks(i,1), peaks(i,2), peaks(i,3), 100, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', color_peaks(i,:))
        end
        axis off
    elseif k == 2
        axh = axes('position',[x_middle y1_middle width_middle height_middle]);
        hold(axh, 'on')
        axis tight, axis square, axis equal
        view(-50,25)
        p = patch('Faces',odf_s.tri,'Vertices',odf_s.verts);
        set(p,'FaceColor','interp', 'FaceVertexCData', odf_s.ct2, 'EdgeColor', 'none', 'LineWidth',0.01)
        for i = 1:n_peaks
            scatter3(peaks(i,1), peaks(i,2), peaks(i,3), 100, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', color_peaks(i,:))
        end
        axis off
    end
    
    if k == 1
        param_names = {'diso'; 'sqddelta'; 'r2'};
        param_names_plot = {'$D_\mathrm{iso}\,(\mu\mathrm{m}^2/\mathrm{ms})$'; '$D_\Delta^2$'; '$R_2\,(\mathrm{s}^{-1})$'};
        Nparam = numel(param_names);
        
        axh_diso = axes('position',[x_right y3_right width_right height_right],'YAxisLocation','right');
        axh_sqddelta = axes('position',[x_right y2_right width_right height_right],'YAxisLocation','right');
        axh_r2 = axes('position',[x_right y1_right width_right height_right],'YAxisLocation','right');
    elseif k == 2
        param_names = {'diso'; 'sqddelta'; 't2'};
        param_names_plot = {'$D_\mathrm{iso}\,(\mu\mathrm{m}^2/\mathrm{ms})$'; '$D_\Delta^2$'; '$T_2\,(\mathrm{ms})$'};
        Nparam = numel(param_names);
        
        axh_diso = axes('position',[x_right y3_right width_right height_right],'YAxisLocation','right');
        axh_sqddelta = axes('position',[x_right y2_right width_right height_right],'YAxisLocation','right');
        axh_t2 = axes('position',[x_right y1_right width_right height_right],'YAxisLocation','right');
    end
    
    for nparam = 1:Nparam
        eval(['axh = axh_' param_names{nparam} ';'])
        hold(axh,'on')
        
        axh.XAxis.TickLabelInterpreter = 'latex';
        axh.YAxis.TickLabelInterpreter = 'latex';
        
        eval(['quantity = ' param_names{nparam} '_peaks;'])
        if strcmp(param_names{nparam},'diso')
            quantity = quantity*1e9;
        elseif strcmp(param_names{nparam},'t2')
            quantity = quantity*1e3;
        end
        
        for i = 1:n_peaks
            color = color_peaks(i,:);
            scatter(axh, i, quantity(i), 40, 'filled', 'MarkerFaceColor', color, 'MarkerEdgeColor', color)
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
        set(ylab, 'Units', 'Normalized', 'Position', [1.38, 0.5, 0]);
        
        set(axh, 'LineWidth', lw)
        set(axh, 'TickDir','out');
    end
    
    saveas(gcf,fullfile(directory, strcat('ODF_metrics',num2str(k),'_nbMeshNodes', num2str(nb_MeshNodes), '_vx', num2str(vx), '_vy', num2str(vy), '_vz', num2str(vz), '.pdf')))
    clf(f)
end

end