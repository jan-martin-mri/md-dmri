function weight_median = weighted_median(list, w)
    [list_sorted, ind_sorted] = sort(list);
    w = w/sum(w);
    w = w(ind_sorted);
    c = cumsum(w);
    i = find(c >= 0.5, 1, 'first');
    if abs(c(i) - 0.5) <= abs(c(i-1) - 0.5)
        weight_median = list_sorted(i);
    else
        weight_median = list_sorted(i-1);
    end
end