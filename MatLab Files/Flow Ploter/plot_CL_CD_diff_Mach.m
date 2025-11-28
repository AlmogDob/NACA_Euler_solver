clc; clear; close all;

db_file = 'NACA.db';
db_table = sqlread(sqlite(db_file),'NACA_data');
%%
NACA = '9410';

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







