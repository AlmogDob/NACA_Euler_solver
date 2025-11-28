clc; clear; close all;

db_file = 'NACA.db';
db_table = sqlread(sqlite(db_file),'NACA_data');

%%

naca_list = unique(fetch(sqlite(db_file),'SELECT NACA FROM NACA_data;'));
naca_list = naca_list.NACA(:)

NACA = '9410';
alpha_deg = 4;

sql_temp = sprintf('SELECT * FROM NACA_data where NACA = %s and angle_of_attack_deg = %f and delta_t <= 1e-5;', NACA, alpha_deg)
results = fetch(sqlite(db_file),sql_temp);
results = sortrows(results, 'Mach_inf');
IDs = results.ID;

fig1 = figure('Name','1','Position',[300,150,900,600]);
font_size = 16;
what_to_show = 'M';
for ID = IDs'
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
        title(sprintf('NACA: %s\nMach: %.2f,  $\\alpha=%.2f^\\circ$', NACA, Mach_inf, alpha_deg), 'FontSize',font_size,'Interpreter','latex')
    elseif what_to_show == 'p0'
        contourf(x, y, p0, 200, 'LineStyle','none');
        hold on
        cbar = colorbar;
        cbar.Label.String = 'p0';
        colormap('turbo');
        axis equal
        title(sprintf('NACA: %s,  Mach: %.2f,  $\\alpha=%.2f^\\circ$', NACA, Mach_inf, alpha_deg), 'FontSize',font_size,'Interpreter','latex')
    elseif what_to_show == 'p'
        contourf(x, y, p, 200, 'LineStyle','none');
        hold on
        cbar = colorbar;
        cbar.Label.String = 'p';
        colormap('turbo');
        axis equal
        title(sprintf('NACA: %s,  Mach: %.2f,  $\\alpha=%.2f^\\circ$', NACA, Mach_inf, alpha_deg), 'FontSize',font_size,'Interpreter','latex')
    end
    
    
    
    % quiver(x, y, u, v, "Color", "#000000");

    drawnow
    % input("press");

    pause(0.5);
end







