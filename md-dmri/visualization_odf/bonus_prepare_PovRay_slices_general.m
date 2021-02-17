%% Adapt this script for dtr2r1d (take a look at .r and .t variables).
%% Copy / Paste ODFplot_r and ODFplot_t and rename to _r1 and _r2 in PovRay

clear

data_path = '/Users/alexis_reymbaut/Dropbox/Research_Lund/Pipeline_general/Glioma_Spectrum/20200504083444';
povray_files_path = '/Users/alexis_reymbaut/Dropbox/Research_Lund/Pipeline_general/Pipeline/PovRay';
odf_path = fullfile(data_path, '3_odfs');
avg = 'median'; % median or mean
N_nodes = 1000;
method = 'dtd';

FLAG_relaxation = false;
if strcmp(method,'dtr2d') || strcmp(method,'dtr1d')
    FLAG_relaxation = true;
end

%% Load ODF file
odf_name = ['odf_bs' avg];
odf_fn = fullfile(odf_path, strcat('odf_bs', avg, '_', num2str(N_nodes), '.mat'));
load(odf_fn, odf_name);
eval(['odf = ' odf_name ';']);
sz = size(odf.w_bin{1});
norm_resc = 1.75; % axial - 1.5, 2.25; coronal - 1.75, 2.5, 2.75 re-scale odfs

%% Choose slices to be plotted

slice = 'axial';
dim_1 = 30:30; 
dim_2 = 1:sz(1); 
dim_3 = 1:sz(2); 

% slice = 'coronal';
% dim_1 = 1:sz(2); % 39
% dim_2 = 1:sz(1); 
% dim_3 = 1:sz(3); 

% slice = 'sagittal';
% dim_1 = 1:sz(1); % 38 and 44
% dim_2 = 1:sz(2);  
% dim_3 = 1:sz(3); 

%% Color limits
disomin = 1e-10; disomax = 3.5e-9;
sqddeltamin = 0.25; sqddeltamax = 1;

if strcmp(method,'dtr2d')
    rmin = 8; rmax = 30;
    tmin = .03; tmax = .12;
    char_r = 'r2';
    char_t = 't2';
elseif strcmp(method,'dtr1d')
    rmin = 0.1; rmax = 0.9;
    tmin = 1; tmax = 2;
    char_r = 'r1';
    char_t = 't1';
end

%% Create and save files for POVRAY
odf_norm = max(odf.w_bin{1}(:))/norm_resc;
path_pov = fullfile(odf_path, 'pov', slice);
msf_mkdir(path_pov);

if strcmp(slice,'axial') == 1
    for nk = dim_1
        nslice = nk;
        output_path = fullfile(path_pov, num2str(nslice));
        msf_mkdir(output_path);
        
        strings_to_replace = {'~~Ni~~'; '~~Nj~~'; '~~Nk~~'; '~~nslice~~'; '~~lookat1~~'; '~~lookat2~~'};
        replacement_strings = {num2str(sz(1)); num2str(sz(2)); num2str(sz(3)); num2str(nslice); num2str(floor(sz(1)/2)+0.5); num2str(floor(sz(2)/2)+0.5)};
        replace_string(fullfile(povray_files_path, slice, 'ODFplot_orientations.pov'), fullfile(output_path, 'ODFplot_orientations.pov'), strings_to_replace, replacement_strings);
        replace_string(fullfile(povray_files_path, slice, 'ODFplot_diso.pov'), fullfile(output_path, 'ODFplot_diso.pov'), strings_to_replace, replacement_strings);
        replace_string(fullfile(povray_files_path, slice, 'ODFplot_sqddelta.pov'), fullfile(output_path, 'ODFplot_sqddelta.pov'), strings_to_replace, replacement_strings);
        
        if FLAG_relaxation
            strings_to_replace = {'~~Ni~~'; '~~Nj~~'; '~~Nk~~'; '~~nslice~~'; '~~lookat1~~'; '~~lookat2~~'; '~~char_r~~'; '~~char_t~~'};
            replacement_strings = {num2str(sz(1)); num2str(sz(2)); num2str(sz(3)); num2str(nslice); num2str(floor(sz(1)/2)+0.5); num2str(floor(sz(2)/2)+0.5); char_r; char_t};
            replace_string(fullfile(povray_files_path, slice, 'ODFplot_r.pov'), fullfile(output_path, strcat('ODFplot_', char_r, '.pov')), strings_to_replace, replacement_strings);
            replace_string(fullfile(povray_files_path, slice, 'ODFplot_t.pov'), fullfile(output_path, strcat('ODFplot_', char_t, '.pov')), strings_to_replace, replacement_strings);
        end
        
        count = 0;
        clear odf_2
        odf_array.verts = [];
        odf_array.norms = [];
        odf_array.tri = [];
        odf_array.c = [];
        odf_array.diso = [];
        odf_array.sqddelta = [];
        
        if FLAG_relaxation
            odf_array.r = [];
            odf_array.t = [];
        end
        
        for ni = dim_2
            for nj = dim_3
                odf_2.n = odf.n;
                odf_2.x = odf.x;
                odf_2.y = odf.y;
                odf_2.z = odf.z;
                odf_2.c = abs([odf.x odf.y odf.z]);
                odf_2.c = odf_2.c./repmat(max(odf_2.c,[],2),[1 3]);
                odf_2.tri = odf.tri;
                
                odf_2.w = squeeze(odf.w_bin{1}(ni,nj,nk,:))/odf_norm;
                odf_2.diso = squeeze(odf.diso_bin{1}(ni,nj,nk,:));
                odf_2.sqddelta = squeeze(odf.sqddelta_bin{1}(ni,nj,nk,:));
                if strcmp(method, 'dtr2d')
                    odf_2.r = squeeze(odf.r2_bin{1}(ni,nj,nk,:));
                    odf_2.t = squeeze(odf.t2_bin{1}(ni,nj,nk,:));
                elseif strcmp(method, 'dtr1d')
                    odf_2.r = squeeze(odf.r1_bin{1}(ni,nj,nk,:));
                    odf_2.t = squeeze(odf.t1_bin{1}(ni,nj,nk,:));
                end
                
                if sum(odf_2.w)>0
                    odf_2.verts = repmat(odf_2.w,[1 3]).*[odf_2.x odf_2.y odf_2.z];
                    odf_2.verts = odf_2.verts + [ni*ones(odf.n,1) nj*ones(odf.n,1) nk*ones(odf.n,1)];
                    odf_2.norms = vertexNormal(triangulation(odf_2.tri,odf_2.verts),(1:odf_2.n)');
                    
                    odf_array.verts = cat(1,odf_array.verts,odf_2.verts);
                    odf_array.tri = cat(1,odf_array.tri,odf_2.tri+count*odf.n);
                    odf_array.c = cat(1,odf_array.c,odf_2.c);
                    odf_array.norms = cat(1,odf_array.norms,odf_2.norms);
                    odf_array.diso = cat(1,odf_array.diso,odf_2.diso);
                    odf_array.sqddelta = cat(1,odf_array.sqddelta,odf_2.sqddelta);
                    
                    if FLAG_relaxation
                        odf_array.r = cat(1,odf_array.r,odf_2.r);
                        odf_array.t = cat(1,odf_array.t,odf_2.t);
                    end
                    
                    count = count+1;
                end
            end
        end
        
        cind = (odf_array.diso-min(disomin))/(max(disomax)-min(disomin));
        col = dist_cind2rgb_jet(cind);
        odf_array.cdiso = [col.r col.g col.b];
        
        cind = (odf_array.sqddelta-min(sqddeltamin))/(max(sqddeltamax)-min(sqddeltamin));
        col = dist_cind2rgb_jet(cind);
        odf_array.csqddelta = [col.r col.g col.b];
        
        if FLAG_relaxation
            cind = (odf_array.r-min(rmin))/(max(rmax)-min(rmin));
            col = dist_cind2rgb_jet(cind);
            odf_array.cr = [col.r col.g col.b];
            
            cind = (odf_array.t-min(tmin))/(max(tmax)-min(tmin));
            col = dist_cind2rgb_jet(cind);
            odf_array.ct = [col.r col.g col.b];
        end
        
        fid = fopen(fullfile(output_path, strcat('ODF_verts_', num2str(nslice), '.txt')), 'w');
        N = size(odf_array.verts, 1);
        fprintf(fid, '%8.0i,\n', N);
        format = '%8.2f, %8.2f, %8.2f,\n';
        for n = 1:N
            fprintf(fid,format,odf_array.verts(n,:));
        end
        fclose(fid);
        
        fid = fopen(fullfile(output_path, strcat('ODF_tri_', num2str(nslice), '.txt')), 'w');
        N = size(odf_array.tri, 1);
        fprintf(fid, '%8.0i,\n', N);
        format = '%8.0f, %8.0f, %8.0f,\n';
        for n = 1:N
            fprintf(fid,format,odf_array.tri(n,:)-1);
        end
        fclose(fid);
        
        fid = fopen(fullfile(output_path, strcat('ODF_norms_', num2str(nslice), '.txt')), 'w');
        N = size(odf_array.norms, 1);
        fprintf(fid, '%8.0i,\n', N);
        format = '%8.2f, %8.2f, %8.2f,\n';
        for n = 1:N
            fprintf(fid,format,odf_array.norms(n,:));
        end
        fclose(fid);
        
        fid = fopen(fullfile(output_path, strcat('ODF_c_', num2str(nslice), '.txt')), 'w');
        N = size(odf_array.c, 1);
        fprintf(fid, '%8.0i,\n', N);
        format = '%8.3f, %8.3f, %8.3f,\n';
        for n = 1:N
            fprintf(fid,format,odf_array.c(n,:));
        end
        fclose(fid);
        
        fid = fopen(fullfile(output_path, strcat('ODF_cdiso_', num2str(nslice), '.txt')), 'w');
        N = size(odf_array.cdiso, 1);
        fprintf(fid, '%8.0i,\n', N);
        format = '%8.3f, %8.3f, %8.3f,\n';
        for n = 1:N
            fprintf(fid,format,odf_array.cdiso(n,:));
        end
        fclose(fid);
        
        fid = fopen(fullfile(output_path, strcat('ODF_csqddelta_', num2str(nslice), '.txt')), 'w');
        N = size(odf_array.csqddelta, 1);
        fprintf(fid, '%8.0i,\n', N);
        format = '%8.3f, %8.3f, %8.3f,\n';
        for n = 1:N
            fprintf(fid,format,odf_array.csqddelta(n,:));
        end
        fclose(fid);
        
        if FLAG_relaxation
            fid = fopen(fullfile(output_path, strcat('ODF_c', char_r, '_', num2str(nslice), '.txt')), 'w');
            N = size(odf_array.cr, 1);
            fprintf(fid, '%8.0i,\n', N);
            format = '%8.3f, %8.3f, %8.3f,\n';
            for n = 1:N
                fprintf(fid,format,odf_array.cr(n,:));
            end
            fclose(fid);
            
            fid = fopen(fullfile(output_path, strcat('ODF_c', char_t, '_', num2str(nslice), '.txt')), 'w');
            N = size(odf_array.ct, 1);
            fprintf(fid, '%8.0i,\n', N);
            format = '%8.3f, %8.3f, %8.3f,\n';
            for n = 1:N
                fprintf(fid,format,odf_array.ct(n,:));
            end
            fclose(fid);
        end     
    end
    
elseif strcmp(slice,'coronal')
    for nj = dim_1
        nslice = nj;
        output_path = fullfile(path_pov, num2str(nslice));
        msf_mkdir(output_path);
        
        strings_to_replace = {'~~Ni~~'; '~~Nj~~'; '~~Nk~~'; '~~nslice~~'; '~~lookat1~~'; '~~lookat2~~'};
        replacement_strings = {num2str(sz(1)); num2str(sz(2)); num2str(sz(3)); num2str(nslice); num2str(floor(sz(3)/2)+0.5); num2str(floor(sz(1)/2)+0.5)};
        replace_string(fullfile(povray_files_path, slice, 'ODFplot_orientations.pov'), fullfile(output_path, 'ODFplot_orientations.pov'), strings_to_replace, replacement_strings);
        replace_string(fullfile(povray_files_path, slice, 'ODFplot_diso.pov'), fullfile(output_path, 'ODFplot_diso.pov'), strings_to_replace, replacement_strings);
        replace_string(fullfile(povray_files_path, slice, 'ODFplot_sqddelta.pov'), fullfile(output_path, 'ODFplot_sqddelta.pov'), strings_to_replace, replacement_strings);
        
        if FLAG_relaxation
            strings_to_replace = {'~~Ni~~'; '~~Nj~~'; '~~Nk~~'; '~~nslice~~'; '~~lookat1~~'; '~~lookat2~~'; '~~char_r~~'; '~~char_t~~'};
            replacement_strings = {num2str(sz(1)); num2str(sz(2)); num2str(sz(3)); num2str(nslice); num2str(floor(sz(3)/2)+0.5); num2str(floor(sz(1)/2)+0.5); char_r; char_t};
            replace_string(fullfile(povray_files_path, slice, 'ODFplot_r.pov'), fullfile(output_path, strcat('ODFplot_', char_r, '.pov')), strings_to_replace, replacement_strings);
            replace_string(fullfile(povray_files_path, slice, 'ODFplot_t.pov'), fullfile(output_path, strcat('ODFplot_', char_t, '.pov')), strings_to_replace, replacement_strings);
        end
        
        count = 0;
        clear odf_2
        odf_array.verts = [];
        odf_array.norms = [];
        odf_array.tri = [];
        odf_array.c = [];
        odf_array.diso = [];
        odf_array.sqddelta = [];
        
        if FLAG_relaxation
            odf_array.r = [];
            odf_array.t = [];
        end
        
        for ni = dim_2
            for nk = dim_3      
                odf_2.n = odf.n;
                odf_2.x = odf.x;
                odf_2.y = odf.y;
                odf_2.z = odf.z;
                odf_2.c = abs([odf.x odf.y odf.z]);
                odf_2.c = odf_2.c./repmat(max(odf_2.c,[],2),[1 3]);
                odf_2.tri = odf.tri;
                
                odf_2.w = squeeze(odf.w_bin{1}(ni,nj,nk,:))/odf_norm;
                odf_2.diso = squeeze(odf.diso_bin{1}(ni,nj,nk,:));
                odf_2.sqddelta = squeeze(odf.sqddelta_bin{1}(ni,nj,nk,:));
                if strcmp(method, 'dtr2d')
                    odf_2.r = squeeze(odf.r2_bin{1}(ni,nj,nk,:));
                    odf_2.t = squeeze(odf.t2_bin{1}(ni,nj,nk,:));
                elseif strcmp(method, 'dtr1d')
                    odf_2.r = squeeze(odf.r1_bin{1}(ni,nj,nk,:));
                    odf_2.t = squeeze(odf.t1_bin{1}(ni,nj,nk,:));
                end
                
                if sum(odf_2.w)>0
                    odf_2.verts = repmat(odf_2.w,[1 3]).*[odf_2.x odf_2.y odf_2.z];
                    odf_2.verts = odf_2.verts + [ni*ones(odf.n,1) nj*ones(odf.n,1) nk*ones(odf.n,1)];
                    odf_2.norms = vertexNormal(triangulation(odf_2.tri,odf_2.verts),(1:odf_2.n)');
                    
                    odf_array.verts = cat(1,odf_array.verts,odf_2.verts);
                    odf_array.tri = cat(1,odf_array.tri,odf_2.tri+count*odf.n);
                    odf_array.c = cat(1,odf_array.c,odf_2.c);
                    odf_array.norms = cat(1,odf_array.norms,odf_2.norms);
                    odf_array.diso = cat(1,odf_array.diso,odf_2.diso);
                    odf_array.sqddelta = cat(1,odf_array.sqddelta,odf_2.sqddelta);
                    
                    if FLAG_relaxation
                        odf_array.r = cat(1,odf_array.r,odf_2.r);
                        odf_array.t = cat(1,odf_array.t,odf_2.t);
                    end
                    
                    count = count+1;
                end
            end
        end
        
        cind = (odf_array.diso-min(disomin))/(max(disomax)-min(disomin));
        col = dist_cind2rgb_jet(cind);
        odf_array.cdiso = [col.r col.g col.b];
        
        cind = (odf_array.sqddelta-min(sqddeltamin))/(max(sqddeltamax)-min(sqddeltamin));
        col = dist_cind2rgb_jet(cind);
        odf_array.csqddelta = [col.r col.g col.b];
        
        if FLAG_relaxation
            cind = (odf_array.r-min(rmin))/(max(rmax)-min(rmin));
            col = dist_cind2rgb_jet(cind);
            odf_array.cr = [col.r col.g col.b];
            
            cind = (odf_array.t-min(tmin))/(max(tmax)-min(tmin));
            col = dist_cind2rgb_jet(cind);
            odf_array.ct = [col.r col.g col.b];
        end
        
        fid = fopen(fullfile(output_path, strcat('ODF_verts_', num2str(nslice), '.txt')), 'w');
        N = size(odf_array.verts, 1);
        fprintf(fid, '%8.0i,\n', N);
        format = '%8.2f, %8.2f, %8.2f,\n';
        for n = 1:N
            fprintf(fid,format,odf_array.verts(n,:));
        end
        fclose(fid);
        
        fid = fopen(fullfile(output_path, strcat('ODF_tri_', num2str(nslice), '.txt')), 'w');
        N = size(odf_array.tri, 1);
        fprintf(fid, '%8.0i,\n', N);
        format = '%8.0f, %8.0f, %8.0f,\n';
        for n = 1:N
            fprintf(fid,format,odf_array.tri(n,:)-1);
        end
        fclose(fid);
        
        fid = fopen(fullfile(output_path, strcat('ODF_norms_', num2str(nslice), '.txt')), 'w');
        N = size(odf_array.norms, 1);
        fprintf(fid, '%8.0i,\n', N);
        format = '%8.2f, %8.2f, %8.2f,\n';
        for n = 1:N
            fprintf(fid,format,odf_array.norms(n,:));
        end
        fclose(fid);
        
        fid = fopen(fullfile(output_path, strcat('ODF_c_', num2str(nslice), '.txt')), 'w');
        N = size(odf_array.c, 1);
        fprintf(fid, '%8.0i,\n', N);
        format = '%8.3f, %8.3f, %8.3f,\n';
        for n = 1:N
            fprintf(fid,format,odf_array.c(n,:));
        end
        fclose(fid);
        
        fid = fopen(fullfile(output_path, strcat('ODF_cdiso_', num2str(nslice), '.txt')), 'w');
        N = size(odf_array.cdiso, 1);
        fprintf(fid, '%8.0i,\n', N);
        format = '%8.3f, %8.3f, %8.3f,\n';
        for n = 1:N
            fprintf(fid,format,odf_array.cdiso(n,:));
        end
        fclose(fid);
        
        fid = fopen(fullfile(output_path, strcat('ODF_csqddelta_', num2str(nslice), '.txt')), 'w');
        N = size(odf_array.csqddelta, 1);
        fprintf(fid, '%8.0i,\n', N);
        format = '%8.3f, %8.3f, %8.3f,\n';
        for n = 1:N
            fprintf(fid,format,odf_array.csqddelta(n,:));
        end
        fclose(fid);
        
        if FLAG_relaxation
            fid = fopen(fullfile(output_path, strcat('ODF_c', char_r, '_', num2str(nslice), '.txt')), 'w');
            N = size(odf_array.cr, 1);
            fprintf(fid, '%8.0i,\n', N);
            format = '%8.3f, %8.3f, %8.3f,\n';
            for n = 1:N
                fprintf(fid,format,odf_array.cr(n,:));
            end
            fclose(fid);
            
            fid = fopen(fullfile(output_path, strcat('ODF_c', char_t, '_', num2str(nslice), '.txt')), 'w');
            N = size(odf_array.ct, 1);
            fprintf(fid, '%8.0i,\n', N);
            format = '%8.3f, %8.3f, %8.3f,\n';
            for n = 1:N
                fprintf(fid,format,odf_array.ct(n,:));
            end
            fclose(fid);
        end  
    end
    
elseif strcmp(slice,'sagittal')
    for ni = dim_1
        nslice = ni;
        output_path = fullfile(path_pov, num2str(nslice));
        msf_mkdir(output_path);
        
        strings_to_replace = {'~~Ni~~'; '~~Nj~~'; '~~Nk~~'; '~~nslice~~'; '~~lookat1~~'; '~~lookat2~~'};
        replacement_strings = {num2str(sz(1)); num2str(sz(2)); num2str(sz(3)); num2str(nslice); num2str(floor(sz(2)/2)+0.5); num2str(floor(sz(3)/2)+0.5)};
        replace_string(fullfile(povray_files_path, slice, 'ODFplot_orientations.pov'), fullfile(output_path, 'ODFplot_orientations.pov'), strings_to_replace, replacement_strings);
        replace_string(fullfile(povray_files_path, slice, 'ODFplot_diso.pov'), fullfile(output_path, 'ODFplot_diso.pov'), strings_to_replace, replacement_strings);
        replace_string(fullfile(povray_files_path, slice, 'ODFplot_sqddelta.pov'), fullfile(output_path, 'ODFplot_sqddelta.pov'), strings_to_replace, replacement_strings);
        
        if FLAG_relaxation
            strings_to_replace = {'~~Ni~~'; '~~Nj~~'; '~~Nk~~'; '~~nslice~~'; '~~lookat1~~'; '~~lookat2~~'; '~~char_r~~'; '~~char_t~~'};
            replacement_strings = {num2str(sz(1)); num2str(sz(2)); num2str(sz(3)); num2str(nslice); num2str(floor(sz(2)/2)+0.5); num2str(floor(sz(3)/2)+0.5); char_r; char_t};
            replace_string(fullfile(povray_files_path, slice, 'ODFplot_r.pov'), fullfile(output_path, strcat('ODFplot_', char_r, '.pov')), strings_to_replace, replacement_strings);
            replace_string(fullfile(povray_files_path, slice, 'ODFplot_t.pov'), fullfile(output_path, strcat('ODFplot_', char_t, '.pov')), strings_to_replace, replacement_strings);
        end
        
        count = 0;
        clear odf_2
        odf_array.verts = [];
        odf_array.norms = [];
        odf_array.tri = [];
        odf_array.c = [];
        odf_array.diso = [];
        odf_array.sqddelta = [];
        
        if FLAG_relaxation
            odf_array.r = [];
            odf_array.t = [];
        end
        
        for nj = dim_2
            for nk = dim_3
                odf_2.n = odf.n;
                odf_2.x = odf.x;
                odf_2.y = odf.y;
                odf_2.z = odf.z;
                odf_2.c = abs([odf.x odf.y odf.z]);
                odf_2.c = odf_2.c./repmat(max(odf_2.c,[],2),[1 3]);
                odf_2.tri = odf.tri;
                
                odf_2.w = squeeze(odf.w_bin{1}(ni,nj,nk,:))/odf_norm;
                odf_2.diso = squeeze(odf.diso_bin{1}(ni,nj,nk,:));
                odf_2.sqddelta = squeeze(odf.sqddelta_bin{1}(ni,nj,nk,:));
                if strcmp(method, 'dtr2d')
                    odf_2.r = squeeze(odf.r2_bin{1}(ni,nj,nk,:));
                    odf_2.t = squeeze(odf.t2_bin{1}(ni,nj,nk,:));
                elseif strcmp(method, 'dtr1d')
                    odf_2.r = squeeze(odf.r1_bin{1}(ni,nj,nk,:));
                    odf_2.t = squeeze(odf.t1_bin{1}(ni,nj,nk,:));
                end
                
                if sum(odf_2.w)>0
                    
                    odf_2.verts = repmat(odf_2.w,[1 3]).*[odf_2.x odf_2.y odf_2.z];
                    odf_2.verts = odf_2.verts + [ni*ones(odf.n,1) nj*ones(odf.n,1) nk*ones(odf.n,1)];
                    odf_2.norms = vertexNormal(triangulation(odf_2.tri,odf_2.verts),(1:odf_2.n)');
                    
                    odf_array.verts = cat(1,odf_array.verts,odf_2.verts);
                    odf_array.tri = cat(1,odf_array.tri,odf_2.tri+count*odf.n);
                    odf_array.c = cat(1,odf_array.c,odf_2.c);
                    odf_array.norms = cat(1,odf_array.norms,odf_2.norms);
                    odf_array.diso = cat(1,odf_array.diso,odf_2.diso);
                    odf_array.sqddelta = cat(1,odf_array.sqddelta,odf_2.sqddelta);
                    
                    if FLAG_relaxation
                        odf_array.r = cat(1,odf_array.r,odf_2.r);
                        odf_array.t = cat(1,odf_array.t,odf_2.t);
                    end
                    
                    count = count+1;
                end
            end
        end
        
        cind = (odf_array.diso-min(disomin))/(max(disomax)-min(disomin));
        col = dist_cind2rgb_jet(cind);
        odf_array.cdiso = [col.r col.g col.b];
        
        cind = (odf_array.sqddelta-min(sqddeltamin))/(max(sqddeltamax)-min(sqddeltamin));
        col = dist_cind2rgb_jet(cind);
        odf_array.csqddelta = [col.r col.g col.b];
        
        if FLAG_relaxation
            cind = (odf_array.r-min(rmin))/(max(rmax)-min(rmin));
            col = dist_cind2rgb_jet(cind);
            odf_array.cr = [col.r col.g col.b];
            
            cind = (odf_array.t-min(tmin))/(max(tmax)-min(tmin));
            col = dist_cind2rgb_jet(cind);
            odf_array.ct = [col.r col.g col.b];
        end
        
        fid = fopen(fullfile(output_path, strcat('ODF_verts_', num2str(nslice), '.txt')), 'w');
        N = size(odf_array.verts, 1);
        fprintf(fid, '%8.0i,\n', N);
        format = '%8.2f, %8.2f, %8.2f,\n';
        for n = 1:N
            fprintf(fid,format,odf_array.verts(n,:));
        end
        fclose(fid);
        
        fid = fopen(fullfile(output_path, strcat('ODF_tri_', num2str(nslice), '.txt')), 'w');
        N = size(odf_array.tri, 1);
        fprintf(fid, '%8.0i,\n', N);
        format = '%8.0f, %8.0f, %8.0f,\n';
        for n = 1:N
            fprintf(fid,format,odf_array.tri(n,:)-1);
        end
        fclose(fid);
        
        fid = fopen(fullfile(output_path, strcat('ODF_norms_', num2str(nslice), '.txt')), 'w');
        N = size(odf_array.norms, 1);
        fprintf(fid, '%8.0i,\n', N);
        format = '%8.2f, %8.2f, %8.2f,\n';
        for n = 1:N
            fprintf(fid,format,odf_array.norms(n,:));
        end
        fclose(fid);
        
        fid = fopen(fullfile(output_path, strcat('ODF_c_', num2str(nslice), '.txt')), 'w');
        N = size(odf_array.c, 1);
        fprintf(fid, '%8.0i,\n', N);
        format = '%8.3f, %8.3f, %8.3f,\n';
        for n = 1:N
            fprintf(fid,format,odf_array.c(n,:));
        end
        fclose(fid);
        
        fid = fopen(fullfile(output_path, strcat('ODF_cdiso_', num2str(nslice), '.txt')), 'w');
        N = size(odf_array.cdiso, 1);
        fprintf(fid, '%8.0i,\n', N);
        format = '%8.3f, %8.3f, %8.3f,\n';
        for n = 1:N
            fprintf(fid,format,odf_array.cdiso(n,:));
        end
        fclose(fid);
        
        fid = fopen(fullfile(output_path, strcat('ODF_csqddelta_', num2str(nslice), '.txt')), 'w');
        N = size(odf_array.csqddelta, 1);
        fprintf(fid, '%8.0i,\n', N);
        format = '%8.3f, %8.3f, %8.3f,\n';
        for n = 1:N
            fprintf(fid,format,odf_array.csqddelta(n,:));
        end
        fclose(fid);
        
        if FLAG_relaxation
            fid = fopen(fullfile(output_path, strcat('ODF_c', char_r, '_', num2str(nslice), '.txt')), 'w');
            N = size(odf_array.cr, 1);
            fprintf(fid, '%8.0i,\n', N);
            format = '%8.3f, %8.3f, %8.3f,\n';
            for n = 1:N
                fprintf(fid,format,odf_array.cr(n,:));
            end
            fclose(fid);
            
            fid = fopen(fullfile(output_path, strcat('ODF_c', char_t, '_', num2str(nslice), '.txt')), 'w');
            N = size(odf_array.ct, 1);
            fprintf(fid, '%8.0i,\n', N);
            format = '%8.3f, %8.3f, %8.3f,\n';
            for n = 1:N
                fprintf(fid,format,odf_array.ct(n,:));
            end
            fclose(fid);
        end  
    end
else
    fprintf('Available choices for ''nslice'' variable: \n 1) axial\n 2) axial\n 3) sagittal\n')
    return
end

