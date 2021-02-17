function [NCLUST, ind_fibers_all_clusters, ind_fibers_all_halos, rho] = density_peaks_clustering_improved(NCLUST, theta_fibers_clustering, phi_fibers_clustering, w_fibers_clustering, struct_clustering)

max_nb_clusters = struct_clustering.max_nb_clusters;
min_weight_per_cluster = struct_clustering.min_weight_per_cluster;
w = w_fibers_clustering/sum(w_fibers_clustering);

%% Compute inter-point distance
N = length(w);
inter_point_distance = zeros(N);
for i = 1:N-1
    for j = i+1:N
        inter_point_distance(i,j) = acos(abs(cosine_angular_difference(theta_fibers_clustering(i),phi_fibers_clustering(i),theta_fibers_clustering(j),phi_fibers_clustering(j))));
    end
end
inter_point_distance = inter_point_distance + inter_point_distance';

%% Find optimal sigma for the gaussian potential
% sigma_samples = linspace(1e-11,max_sigma,500);
sigma_samples = logspace(-11,-0.5,300);
entropy_list = make_entropy_gaussian(sigma_samples, inter_point_distance, w);
[~, ind_min_entropy] = min(entropy_list);
optimal_sigma = sigma_samples(ind_min_entropy);
d_cutoff = 3*optimal_sigma;

% radius = 2/sqrt(N); % Average radius occupied on the sphere
% d_cutoff = radius; 

%% Compute density using a Gaussian kernel
rho = zeros(N,1);
for i=1:N
    rho(i) = sum(w .* exp(-(inter_point_distance(i,:)/d_cutoff).^2/2));
end
[~, ord_rho] = sort(rho,'descend');

%% Compute delta-distance and nearest neighbors with higher density
delta = zeros(N,1);
nn = zeros(N,1);
for i = 2:N
    t_index = ord_rho(1:i-1);
    [delta(ord_rho(i)), min_idx] = min(inter_point_distance(ord_rho(i),ord_rho(1:i-1)));
    nn(ord_rho(i)) = t_index(min_idx);
end
delta(ord_rho(1)) = max(delta(2:end));

%% Number of clusters (given by the number of ODF peaks)
if NCLUST > max_nb_clusters
    NCLUST = max_nb_clusters;
end

% Using delta (isolate points of largest delta-distance and sufficiently high density)
% threshold_rho = struct_clustering.threshold_rho;
% ind_rho_relevant = (rho >= threshold_rho*max(rho));
% [~, ord_delta] = sort(delta, 'descend'); 
% ord_delta_relevant = ord_delta(ind_rho_relevant(ord_delta));

% Using gamma
gamma = rho/max(rho).*delta/max(delta);
[~, ord_gamma] = sort(gamma, 'descend'); 

%% Clustering
repeat_clustering = 1;
while repeat_clustering && NCLUST > 0
    
    % ind_cluster_centers = ord_delta_relevant(1:NCLUST);
    ind_cluster_centers = ord_gamma(1:NCLUST);
    
    % If highest-density point is not a cluster centroid
    if ~ismember(ord_rho(1), ind_cluster_centers)
        [~, min_idx] = min(inter_point_distance(ord_rho(1), ind_cluster_centers));
        nn(ord_rho(1)) = ind_cluster_centers(min_idx);
    end
    
    %% Cluster centers (centroids)
    cluster_result = zeros(N,1);
    null_index = 1:N;
    centroid_index = null_index(ind_cluster_centers);
    for i = 1:NCLUST
        cluster_result(centroid_index(i)) = i;
    end
    
    %% Assign nearest neighbors with higher density to their cluster center
    for i = 1:N
        if cluster_result(ord_rho(i)) == 0
            cluster_result(ord_rho(i)) = cluster_result(nn(ord_rho(i)));
        end
    end
    
    %% Define cluster-border density for future halo
    rho_border = zeros(NCLUST,1);
    for i = 1:N-1
        for j = i+1:N
            if (cluster_result(i) ~= cluster_result(j)) && (inter_point_distance(i,j) < d_cutoff)
                rho_average = (w(i)*rho(i) + w(j)*rho(j))/(w(i) + w(j));
                rho_border(cluster_result(i)) = max(rho_border(cluster_result(i)), rho_average);
                rho_border(cluster_result(j)) = max(rho_border(cluster_result(j)), rho_average);
            end
        end
    end
    
    %% Define halo
    halo = false(N,1);
    for i = 1:N
        if rho(i) < rho_border(cluster_result(i))
            halo(i) = true;
        end
    end
    
    %% Final step
    ind_fibers_all_halos = null_index(halo); % Fibers in clusters' halos
    ind_fibers_not_in_halos = null_index(~halo);
    cluster_assignation = cluster_result(ind_fibers_not_in_halos);
    ind_fibers_all_clusters = cell([1,NCLUST]); % Fibers in clusters' cores
    total_cluster_weights = zeros([1,NCLUST]);
    sum_weights = 0;
    for i = 1:NCLUST
        ind_fibers_all_clusters{i} = ind_fibers_not_in_halos(cluster_assignation == i);
        total_cluster_weights(i) = sum(w(ind_fibers_all_clusters{i}));
        sum_weights = sum_weights + total_cluster_weights(i);
    end
    
    for i = 1:NCLUST
        if (total_cluster_weights(i)/sum_weights < min_weight_per_cluster) || (length(ind_fibers_all_clusters{i}) < 3) % You need at least three points to later compute medians
            repeat_clustering = 1;
            NCLUST = NCLUST - 1;
            break
        else
            repeat_clustering = 0;
        end
    end
end

if NCLUST == 0 || NCLUST == -1
    ind_fibers_all_clusters = [];
    ind_fibers_all_halos = [];
end
   
end