function entropy_list = make_entropy_gaussian(sigma_samples, inter_point_distance, w_fibers_clustering)

nb_fibers = length(w_fibers_clustering);

entropy_list = zeros(1,length(sigma_samples));
for n = 1:length(sigma_samples)
    sigma_test = sigma_samples(n);
    potential_list = zeros(1,nb_fibers);
    for i = 1:nb_fibers 
        potential_list(i) = sum(w_fibers_clustering .* exp(-(inter_point_distance(i,:)/sigma_test).^2/2)); 
    end
    Z = sum(potential_list);
    entropy_list(n) = -sum(potential_list/Z.*log(potential_list/Z));
end

end