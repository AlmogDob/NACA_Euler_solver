clc; clear; close all;

db_file = 'NACA.db';
db_table = sqlread(sqlite(db_file),'NACA_data');
NACA_list = unique(fetch(sqlite(db_file),'SELECT NACA FROM NACA_data;'));

NACA = '44012';
Mach = 0.9;

sql_temp = sprintf('SELECT * FROM NACA_data where NACA = %s and Mach_inf = %f and delta_t <= 1e-5;', NACA, Mach)
results = fetch(sqlite(db_file),sql_temp);
% results = fetch(sqlite(db_file),sprintf('SELECT * FROM NACA_data where NACA = %s and Mach_inf = %f and delta_t <= 1e-5;', NACA, Mach));
% results = fetch(sqlite(db_file),sprintf('SELECT * FROM NACA_data where NACA = %s and Mach_inf = %f and delta_t = 1e-5 and angle_of_attack_deg = 0;', NACA, Mach));
results = sortrows(results, 'angle_of_attack_deg');
IDs = results.ID;

fig1 = figure('Name','1','Position',[300,150,900,600]);
for ID = IDs'
% for ID = 18
    ID;
    [ni, nj, num_points_on_airfoil, gamma, rho_inf, Mach_inf, p_inf, a_inf, u_inf, alpha_deg, NACA, x, y, rho, u, v, e, p, a, M, p0, CL, CD] = read_matrixes_from_DB(ID, db_table);
    clf(fig1);

    hold all
    % M(M>=1)=1;
    contourf(x, y, M, 200, 'LineStyle','none');
    colorbar;
    colormap('turbo');
    axis equal
    title(sprintf('NACA: %s,  Mach: %.2f,  $\\alpha=%.2f^\\circ$', NACA, Mach_inf, alpha_deg),'Interpreter','latex')
    % quiver(x, y, u, v, "Color", "#000000");

    drawnow
    % input("press");

    % pause(0.5);
end





