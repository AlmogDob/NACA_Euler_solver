clc; clear; close all;

db_file = 'NACA.db';
db_table = sqlread(sqlite(db_file),'NACA_data');
unique(fetch(sqlite(db_file),'SELECT NACA FROM NACA_data;'))

NACA = '54112';
alpha_deg = 0;

sql_temp = sprintf('SELECT * FROM NACA_data where NACA = %s and angle_of_attack_deg = %f and delta_t <= 1e-5;', NACA, alpha_deg)
results = fetch(sqlite(db_file),sql_temp);
results = sortrows(results, 'Mach_inf');
IDs = results.ID;

fig1 = figure('Name','1','Position',[300,150,900,600]);
font_size = 16;
what_to_show = 'M';
for ID = IDs'
% for ID = 18
    ID;
    [ni, nj, num_points_on_airfoil, gamma, rho_inf, Mach_inf, p_inf, a_inf, u_inf, alpha_deg, NACA, x, y, rho, u, v, e, p, a, M, p0, CL, CD] = read_matrixes_from_DB(ID, db_table);
    clf(fig1);

    if what_to_show == 'M'
        contourf(x, y, M, 200, 'LineStyle','none');
        hold on
        cbar = colorbar;
        cbar.Label.String = 'M';
        colormap('turbo');
        axis equal
        title(sprintf('NACA: %s,  Mach: %.2f,  $\\alpha=%.2f^\\circ$', NACA, Mach_inf, alpha_deg), 'FontSize',font_size,'Interpreter','latex')
    elseif what_to_show == 'p0'
        contourf(x, y, p0, 200, 'LineStyle','none');
        hold on
        cbar = colorbar;
        cbar.Label.String = 'p0';
        colormap('turbo');
        axis equal
        title(sprintf('NACA: %s,  Mach: %.2f,  $\\alpha=%.2f^\\circ$', NACA, Mach_inf, alpha_deg), 'FontSize',font_size,'Interpreter','latex')
    
    end
    
    
    
    % quiver(x, y, u, v, "Color", "#000000");

    drawnow
    % input("press");

    pause(0.5);
end






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
