function dtr2d_weights = dtr2d_dist2nodes(dtr2d)

[n,~,~,~,~,~,w] = dtr2d_dist2par(dtr2d);

dtr2d_weights = [w'];
dtr2d_weights = [n; dtr2d_weights(:)];
