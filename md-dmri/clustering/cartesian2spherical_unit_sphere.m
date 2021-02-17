function [theta, phi] = cartesian2spherical_unit_sphere(x,y,z)

[azimuth,elevation,~] = cart2sph(x,y,z);

theta = pi/2 - elevation;

phi = zeros(1,length(azimuth));
for i = 1:length(azimuth)
    angle = azimuth(i);
    if angle < 0
        phi(i) = angle+2*pi;
    else
        phi(i) = angle;
    end
end
