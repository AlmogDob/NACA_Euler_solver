clc; clear; close all;

db_file = 'NACA.db';
db = sqlread(sqlite(db_file),'NACA_data');

results = fetch(sqlite(db_file),'SELECT * FROM NACA_data where NACA = 0012 and Mach_inf = 0.6 and delta_t = 1e-5 and angle_of_attack_deg = 7;');
results = sortrows(results, 'angle_of_attack_deg');
ID = results.ID;

[ni, nj, num_points_on_airfoil, gamma, rho_inf, Mach_inf, p_inf, a_inf, u_inf, alpha_deg, NACA, x, y, rho, u, v, e, p, a, M, p0, CL, CD] = read_matrixes_from_DB(ID, results);

figure
hold all
contourf(x, y, M, 300, 'LineStyle','none');
colorbar;
colormap('turbo');
axis equal
title(sprintf('NACA: %s,  Mach: %.2f,  $\\alpha=%.2f^\\circ$', NACA, Mach_inf, alpha_deg),'Interpreter','latex')
% 
% factor = 0.0001;
% quiver(x, y, u*factor, v*factor,"off", "Color", "#000000");


% plot(x(:,end), y(:,end),'-k')%,'Color',"#7E2F8E")
% plot(x(:,1), y(:,1),'-k')
% plot(x(end,:), y(end,:),'-k')
% plot(x(1,:), y(1,:),'-k')
% 
% for i = 2:(length(x(1,:))-1)
%     plot(x(:,i), y(:,i),'-k')%,'Color',"#0072BD")
% end
% for i = 2:(length(x(:,1))-1)
%     plot(x(i,:), y(i,:),'-k')%,'Color',"#0072BD")
% end

%%

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
tot_force = sum(force_vecs);
tot_coeff = tot_force / (0.5 * rho_inf * u_inf^2)
CL = - tot_coeff(1) * sind(alpha_deg) + tot_coeff(2) * cosd(alpha_deg)
CD = + tot_coeff(1) * cosd(alpha_deg) + tot_coeff(2) * sind(alpha_deg)

figure
quiver(points_on_airfoil(1:end,1), points_on_airfoil(1:end,2), normal_vecs(:,1), normal_vecs(:,2), 'off')
hold on
quiver(points_on_airfoil(1:end,1), points_on_airfoil(1:end,2), parallel_vecs(:,1), parallel_vecs(:,2), 'off')
hold on
plot(points_on_airfoil(1:end,1), points_on_airfoil(1:end,2),'k')
axis equal






