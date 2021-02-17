function list_colors = cluster_orientational_colors(median_x,median_y,median_z)

nb_clusters = length(median_x);

list_colors = cell([1, nb_clusters]);
for n = 1:nb_clusters
    median_r = [median_x(n),median_y(n),median_z(n)];
    list_colors{n} = abs(median_r)/max(abs(median_r));
end

end

% if optimal_K == 1
%     list_colors = red;
% elseif optimal_K == 2
%     list_colors = [red,green];
% elseif optimal_K == 3
%     list_colors = [red,green, blue];
% elseif optimal_K == 4
%     list_colors = [red,green, blue, orange];
% elseif optimal_K == 5
%     list_colors = [red,green, blue, orange, yellow];
% end