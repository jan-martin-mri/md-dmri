function weight_quant = weighted_quantile(list, w, q_list)

    ind_not_nan = logical(~isnan(list).*~isnan(w));
    list = list(ind_not_nan);
    w = w(ind_not_nan);
    
    [list_sorted, ind_sorted] = sort(list);
    w = w/sum(w);
    w = w(ind_sorted);
    c = cumsum(w);
    
    weight_quant = zeros([1 length(q_list)]);
    for m = 1:length(q_list)
        q = q_list(m);
        i = find(c >= q, 1, 'first');
        if i == 1
           weight_quant(m) = list_sorted(1); 
        else
            if abs(c(i) - q) <= abs(c(i-1) - q)
                weight_quant(m) = list_sorted(i);
            else
                weight_quant(m) = list_sorted(i-1);
            end
        end
    end
end