function [dt] = find_dt_please(rho, c, k, h, dx, dy)
% function to find the most restrictive time step from the 6 unique
% criteria

% first, calculate the 5 unique time steps
dt1 = 1 ./ (2./rho./c .* (h./dy + h./dx + k./dy.^2 + k./dx.^2)); % dt criterion for top corners
dt2 = 1 ./ (2./rho./c .* (h./dx + k./dy.^2 + k./dx.^2)); % dt criterion for sides
dt3 = 1 ./ (2./rho./c .* (h./dy + k./dx.^2 + k./dy.^2)); % dt criterion for top
dt4 = 1 ./ (2./rho./c .* (k./dy.^2 + k./dx.^2 + h./dx)); % dt criterion for bottom corners
dt5 = 1 ./ (k./rho./c .* (2./dy.^2 + 3./dx.^2)); % dt criterion for bottom
dt6 = 1 ./ (2.*k./rho./c .* (1./dy.^2 + 1./dx.^2)); % dt criterion for interior

% put them all in an array
dt_all = [dt1; dt2; dt3; dt4; dt5; dt6];

% find the minimum
dt = zeros(length(rho),1);
for k = 1:length(rho)
    dt(k) = min(dt_all(k,:));
end





end