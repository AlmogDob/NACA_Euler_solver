clc; clear; close all;

db_file = 'NACA.db';
db = sqlread(sqlite(db_file),'NACA_data');
%%
wanted_NACA = '0010';
wanted_Mach = '0.6';
wanted_alpha = '1';

results = fetch(sqlite(db_file),sprintf('SELECT * FROM NACA_data where NACA = %s and Mach_inf = %s and delta_t = 1e-5 and angle_of_attack_deg = %s;', wanted_NACA, wanted_Mach, wanted_alpha));
results = sortrows(results, 'angle_of_attack_deg');
ID = results.ID;

[ni, nj, num_points_on_airfoil, gamma, rho_inf, Mach_inf, p_inf, a_inf, u_inf, alpha_deg, NACA, x, y, rho, u, v, e, p, a, M, p0, CL, CD] = read_matrixes_from_DB(ID, results);

figure('Position',[300,100,1200,750])
hold all
contourf(x, y, M, 300, 'LineStyle','none');
colorbar;
colormap('turbo');
axis equal
title(sprintf('NACA: %s,  Mach: %.2f,  $\\alpha=%.2f^\\circ$', NACA, Mach_inf, alpha_deg),'Interpreter','latex')

i_LE  = (ni-1) / 2;
i_TEL = i_LE - num_points_on_airfoil / 2;
i_TEU = i_LE + num_points_on_airfoil / 2;
num_points_half = num_points_on_airfoil / 2 - 1;

points_on_airfoil   = [x(1,i_TEL+2:i_TEU)', y(1,i_TEL+2:i_TEU)'];
p_on_airfoil_points = p(1,i_TEL+2:i_TEU)';

parallel_vecs = [];
normal_vecs = [];
force_vecs = [];
p_on_airfoil_panels = [];

for i = (1:num_points_half)+1
    parallel_vecs(i,1) = points_on_airfoil(i-1,1) - points_on_airfoil(i,1);
    parallel_vecs(i,2) = points_on_airfoil(i-1,2) - points_on_airfoil(i,2);
end
for i = num_points_on_airfoil:-1:(num_points_on_airfoil - num_points_half + 1)
    parallel_vecs(i,1) = points_on_airfoil(i-1,1) - points_on_airfoil(i,1);
    parallel_vecs(i,2) = points_on_airfoil(i-1,2) - points_on_airfoil(i,2);
end
for i = 1:num_points_on_airfoil
    normal_vecs(i, 1) = - parallel_vecs(i, 2);
    normal_vecs(i, 2) = + parallel_vecs(i, 1);
end
for i = 2:num_points_on_airfoil
    p_on_airfoil_panels(i, 1) = (p_on_airfoil_points(i, 1) + p_on_airfoil_points(i-1, 1)) / 2;
end
for i = 1:num_points_on_airfoil
    force_vecs(i, :) = normal_vecs(i, :) * p_on_airfoil_panels(i, 1);
end
force_vecs_srink = force_vecs(2:end,:);
tot_force = sum(force_vecs);
tot_coeff = tot_force / (0.5 * rho_inf * u_inf^2);
CL = - tot_coeff(1) * sind(alpha_deg) + tot_coeff(2) * cosd(alpha_deg);
CD = + tot_coeff(1) * cosd(alpha_deg) + tot_coeff(2) * sind(alpha_deg);

mid_points = [(points_on_airfoil(1:end-1,1)+points_on_airfoil(2:end,1))/2,(points_on_airfoil(1:end-1,2)+points_on_airfoil(2:end,2))/2];

force_magnitudes = sqrt(force_vecs_srink(:,1).^2+force_vecs_srink(:,2).^2);
center_of_force = sum(mid_points .* force_magnitudes, 1) ./ sum(force_magnitudes);

mom_func = @(x_cp, y_cp) sum( ...
    (mid_points(:,1) - x_cp) .* force_vecs_srink(:,2) - ...
    (mid_points(:,2) - y_cp) .* force_vecs_srink(:,1) );

% % --- Compute the line of action (moment/force vector relationship) ---
% % Any point [x_CoP, y_CoP] on the line satisfies:
% % moment_origin = x_CoP*F_total(2) - y_CoP*F_total(1)
% % Solve the least-squares zero-moment condition in 2D
% moment_origin = mom_func(0,0);
% A = [ tot_force(2), -tot_force(1) ];
% b = moment_origin;
% CoP = (A' * b) / (A*A');

% Split the airfoil coordinates
xU = x(1, i_LE+2:i_TEU);   yU = y(1, i_LE+2:i_TEU);
xL = x(1, i_TEL+2:i_LE);   yL = y(1, i_TEL+2:i_LE);
% Ensure they are in same x-direction
if xU(1) > xU(end)
    xU = flip(xU);  yU = flip(yU);
end
if xL(1) > xL(end)
    xL = flip(xL);  yL = flip(yL);
end
% % Interpolate lower y's to match upper x's
% yL_interp = interp1(xL, yL, xU, 'linear', 'extrap');
% Camber line coordinates
x_camber = 0.5 * (xU + xL);
y_camber = 0.5 * (yU + yL);
% Interpolant for camber line shape
y_camber_fun = @(xq) interp1(x_camber, y_camber, xq, 'linear', 'extrap');
% Anonymous function for moment about camber point at given x
moment_on_camber = @(x_cp) mom_func(x_cp, y_camber_fun(x_cp));
% Bracket the interval across the chord (e.g. between LE and TE)
x_interval = [min(x_camber), max(x_camber)];
% Solve M=0 along camber line
x_CoP_camber = fzero(moment_on_camber, mean(x_interval));
y_CoP_camber = y_camber_fun(x_CoP_camber);
CoP = [x_CoP_camber, y_CoP_camber];

fprintf('CL = %g, CD = %g\n', CL, CD)
mom_func(CoP(1), CoP(2))

% figure
hold all
quiver(mid_points(1:end,1), mid_points(1:end,2), force_vecs_srink(:,1), force_vecs_srink(:,2), 'Color', 'k')
plot(points_on_airfoil(1:end,1), points_on_airfoil(1:end,2),'k-+')
axis equal
% plot([CoP(1),0.25], [CoP(2),y_camber_fun(0.25)], '--*m', 'LineWidth',1.5)
% text(0.25, y_camber_fun(0.25), sprintf('\n0.25C'), 'Color', 'm', 'FontSize', 16)
plot(CoP(1), CoP(2), 'ro', 'MarkerFaceColor', 'r')
text(CoP(1), CoP(2), sprintf('CoP\n'), 'Color', 'r', 'FontSize', 16)
plot(x_camber, y_camber)





% figure
% quiver(points_on_airfoil(1:end,1), points_on_airfoil(1:end,2), parallel_vecs(:,1), parallel_vecs(:,2), 'off')
% hold on
% quiver(points_on_airfoil(1:end,1), points_on_airfoil(1:end,2), force_vecs(:,1), force_vecs(:,2))
% hold on
% plot(points_on_airfoil(1:end,1), points_on_airfoil(1:end,2),'k')
% axis equal





