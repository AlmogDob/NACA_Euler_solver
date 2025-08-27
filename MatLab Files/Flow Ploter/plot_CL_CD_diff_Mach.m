clc; clear; close all;

db_file = 'NACA.db';
db_table = sqlread(sqlite(db_file),'NACA_data');

NACA = '0020';

Mach_inf_list = fetch(sqlite(db_file),sprintf('SELECT Mach_inf FROM NACA_data where NACA = %s and delta_t = 1e-5;', NACA));
Mach_inf_list = unique(Mach_inf_list.Mach_inf(:));

fig1 = figure('Name','1 CL_CD','Position',[300,100,1200,750]);
font_size = 16;
colors = hsv(length(Mach_inf_list))*0.9;
subplot(2,2,1)
    hold all
    lg = {};
    for Mach_index = 1:length(Mach_inf_list)
        Mach_inf = Mach_inf_list(Mach_index);
        results = fetch(sqlite(db_file), sprintf('SELECT * FROM NACA_data where NACA = %s and Mach_inf = %f and delta_t = 1e-5;', NACA, Mach_inf));
    
        results = sortrows(results, 'angle_of_attack_deg');
        IDs = results.ID;
        alphas = results.angle_of_attack_deg;
        CLs = results.CL;
        CDs = results.CD;
        title(sprintf('NACA: %s\n', results.NACA(1)))
        
        plot(alphas, CLs, 'LineWidth',1, 'Color', colors(Mach_index,:))
        % plot(alphas, CDs, 'LineWidth',1, 'Color', colors(Mach_index,:))
        % plot(CDs, CLs, 'LineWidth',1, 'Color', colors(Mach_index,:))
        % plot(alphas, CLs./CDs, 'LineWidth',1, 'Color', colors(Mach_index,:))
        lg{end+1} = sprintf('\n%.1f  $%.1f\\le\\alpha\\le %.1f$\n', Mach_inf, alphas(1), alphas(end));
    end
    xlabel('$\alpha$','FontSize',font_size,'Interpreter','latex')
    ylabel('$CL$','FontSize',font_size,'Interpreter','latex')
    legend(lg, 'Location','eastoutside','Interpreter','latex')
    grid on
    grid minor
    box on
subplot(2,2,2)
    hold all
    lg = {};
    for Mach_index = 1:length(Mach_inf_list)
        Mach_inf = Mach_inf_list(Mach_index);
        results = fetch(sqlite(db_file), sprintf('SELECT * FROM NACA_data where NACA = %s and Mach_inf = %f and delta_t = 1e-5;', NACA, Mach_inf));
    
        results = sortrows(results, 'angle_of_attack_deg');
        IDs = results.ID;
        alphas = results.angle_of_attack_deg;
        CLs = results.CL;
        CDs = results.CD;
        title(sprintf('NACA: %s\n', results.NACA(1)))
        
        % plot(alphas, CLs, 'LineWidth',1, 'Color', colors(Mach_index,:))
        plot(alphas, CDs, 'LineWidth',1, 'Color', colors(Mach_index,:))
        % plot(CDs, CLs, 'LineWidth',1, 'Color', colors(Mach_index,:))
        % plot(alphas, CLs./CDs, 'LineWidth',1, 'Color', colors(Mach_index,:))
        lg{end+1} = sprintf('\n%.1f  $%.1f\\le\\alpha\\le %.1f$\n', Mach_inf, alphas(1), alphas(end));
    end
    xlabel('$\alpha$','FontSize',font_size,'Interpreter','latex')
    ylabel('$CD$','FontSize',font_size,'Interpreter','latex')
    % legend(lg, 'Location','eastoutside','Interpreter','latex')
    grid on
    grid minor
    box on
subplot(2,2,3)
    hold all
    lg = {};
    for Mach_index = 1:length(Mach_inf_list)
        Mach_inf = Mach_inf_list(Mach_index);
        results = fetch(sqlite(db_file), sprintf('SELECT * FROM NACA_data where NACA = %s and Mach_inf = %f and delta_t = 1e-5;', NACA, Mach_inf));
    
        results = sortrows(results, 'angle_of_attack_deg');
        IDs = results.ID;
        alphas = results.angle_of_attack_deg;
        CLs = results.CL;
        CDs = results.CD;
        title(sprintf('NACA: %s\n', results.NACA(1)))
        
        % plot(alphas, CLs, 'LineWidth',1, 'Color', colors(Mach_index,:))
        % plot(alphas, CDs, 'LineWidth',1, 'Color', colors(Mach_index,:))
        plot(CDs, CLs, 'LineWidth',1, 'Color', colors(Mach_index,:))
        % plot(alphas, CLs./CDs, 'LineWidth',1, 'Color', colors(Mach_index,:))
        lg{end+1} = sprintf('\n%.1f  $%.1f\\le\\alpha\\le %.1f$\n', Mach_inf, alphas(1), alphas(end));
    end
    xlabel('$CD$','FontSize',font_size,'Interpreter','latex')
    ylabel('$CL$','FontSize',font_size,'Interpreter','latex')
    % legend(lg, 'Location','eastoutside','Interpreter','latex')
    grid on
    grid minor
    box on
subplot(2,2,4)
    hold all
    lg = {};
    for Mach_index = 1:length(Mach_inf_list)
        Mach_inf = Mach_inf_list(Mach_index);
        results = fetch(sqlite(db_file), sprintf('SELECT * FROM NACA_data where NACA = %s and Mach_inf = %f and delta_t = 1e-5;', NACA, Mach_inf));
    
        results = sortrows(results, 'angle_of_attack_deg');
        IDs = results.ID;
        alphas = results.angle_of_attack_deg;
        CLs = results.CL;
        CDs = results.CD;
        title(sprintf('NACA: %s\n', results.NACA(1)))
        
        % plot(alphas, CLs, 'LineWidth',1, 'Color', colors(Mach_index,:))
        % plot(alphas, CDs, 'LineWidth',1, 'Color', colors(Mach_index,:))
        % plot(CDs, CLs, 'LineWidth',1, 'Color', colors(Mach_index,:))
        plot(alphas, CLs./CDs, 'LineWidth',1, 'Color', colors(Mach_index,:))
        lg{end+1} = sprintf('\n%.1f  $%.1f\\le\\alpha\\le %.1f$\n', Mach_inf, alphas(1), alphas(end));
    end
    xlabel('$\alpha$','FontSize',font_size,'Interpreter','latex')
    ylabel('$\frac{CL}{CD}$','FontSize',font_size,'Interpreter','latex')
    % legend(lg, 'Location','eastoutside','Interpreter','latex')
    grid on
    grid minor
    box on
;







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
