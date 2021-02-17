function cosine_beta = cosine_angular_difference(theta_1, phi_1, theta_2, phi_2)

% Spherical law of cosines
cosine_beta = cos(theta_1).*cos(theta_2) + sin(theta_1).*sin(theta_2).*cos(phi_1-phi_2);

% A numerical error could make cosine_beta slightly overshoot its natural bounds 
cosine_beta(cosine_beta > 1) = 1;
cosine_beta(cosine_beta < -1) = -1;
