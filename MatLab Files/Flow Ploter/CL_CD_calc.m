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
title(sprintf('NACA: %04d,  Mach: %.2f,  $\\alpha=%.2f^\\circ$', NACA, Mach_inf, alpha_deg),'Interpreter','latex')
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





function [ni, nj, num_points_on_airfoil, gamma, rho_inf, Mach_inf, p_inf, a_inf, u_inf, alpha_deg, NACA, x, y, rho, u, v, e, p, a, M, p0, CL, CD] = read_matrixes_from_DB(ID, db_table)    
    index = find(db_table.ID == ID);
    x_array = typecast(db_table.x_2Dmat{index},'double');
    y_array = typecast(db_table.y_2Dmat{index},'double');
    rhos    = typecast(db_table.rho_2Dmat{index},'double');
    us      = typecast(db_table.u_2Dmat{index},'double');
    vs      = typecast(db_table.v_2Dmat{index},'double');
    es      = typecast(db_table.e_2Dmat{index},'double');
    
    ni = db_table.ni(index);
    nj = db_table.nj(index);
    num_points_on_airfoil = db_table.num_points_on_airfoil(index);
    gamma     = db_table.Gamma(index);
    rho_inf   = db_table.delta_y(index);
    Mach_inf  = db_table.Mach_inf(index);
    p_inf     = db_table.environment_pressure(index);
    a_inf     = sqrt(gamma * p_inf / rho_inf);
    u_inf     = Mach_inf * a_inf;
    alpha_deg = db_table.angle_of_attack_deg(index);
    NACA      = db_table.NACA(index);
    CL        = db_table.CL(index);
    CD        = db_table.CD(index);
    
    for i = 1:nj
        rho(i,:) = rhos(1+(i-1)*ni:ni+(i-1)*ni);
    end
    rho;
    
    for i = 1:nj
        u(i,:) = us(1+(i-1)*ni:ni+(i-1)*ni);
    end
    u;
    
    for i = 1:nj
        v(i,:) = vs(1+(i-1)*ni:ni+(i-1)*ni);
    end
    v;
    
    for i = 1:nj
        e(i,:) = es(1+(i-1)*ni:ni+(i-1)*ni);
    end
    e;
    
    for i = 1:nj
        x(i,:) = x_array(1+(i-1)*ni:ni+(i-1)*ni);
    end
    x;
    
    for i = 1:nj
        y(i,:) = y_array(1+(i-1)*ni:ni+(i-1)*ni);
    end
    y;
    
    p  = (gamma - 1).*(e - 0.5.*rho.*(u.^2 + v.^2));
    a = sqrt(gamma.*p./rho);
    M = sqrt(u.^2 + v.^2)./a;
    p0 = (p) .* (1 + (gamma-1) / 2 * M.^2).^(gamma / (gamma-1));
end
