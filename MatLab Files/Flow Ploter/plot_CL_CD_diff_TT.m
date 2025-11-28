clc; clear; close all;

% db_file = 'NACA_int.db';
db_file = 'NACA.db';
db_table = sqlread(sqlite(db_file),'NACA_data');
%%

NACA_start = '94';
NACA_end = '';
Mach = 0.9;
NACA_list = unique(fetch(sqlite(db_file),sprintf('SELECT NACA FROM NACA_data where NACA GLOB "%s??" and Mach_inf = %f;', NACA_start, Mach)));
NACA_list = NACA_list.NACA(:);
% % 
% NACA_start = '55';
% NACA_end = '12';
% Mach = 0.9;
% NACA_list = unique(fetch(sqlite(db_file),sprintf('SELECT NACA FROM NACA_data where NACA GLOB "%s*%s" and Mach_inf = %f;', NACA_start, NACA_end, Mach)));
% NACA_list = NACA_list.NACA(:);

% NACA_start = '';
% NACA_end = '412';
% Mach = 0.8;
% NACA_list = unique(fetch(sqlite(db_file),sprintf('SELECT NACA FROM NACA_data where NACA GLOB "?%s" and Mach_inf = %f;', NACA_end, Mach)));
% NACA_list = NACA_list.NACA(:);

fig1 = figure('Name','1 CL_CD','Position',[200,100,1400,750]);
font_size = 16;
sgtitle(sprintf('Mach: %.1f\nNACA: %s * %s',Mach, NACA_start, NACA_end),'FontSize',font_size,'Interpreter','latex');
colors = hsv(length(NACA_list))*0.9;
subplot(2,2,1)
    hold all
    lg = {};
    % Mach_inf_list = 0.6;
    for NACA_index = 1:length(NACA_list)
        NACA = NACA_list(NACA_index);
        results = fetch(sqlite(db_file), sprintf('SELECT * FROM NACA_data where NACA = %s and Mach_inf = %f and delta_t <= 1e-5;', NACA, Mach));
    
        results = sortrows(results, 'angle_of_attack_deg');
        IDs = results.ID;
        alphas = results.angle_of_attack_deg;
        CLs = results.CL;
        CDs = results.CD;
        
        plot(alphas, CLs, 'LineWidth',1, 'Color', colors(NACA_index,:))
        % plot(alphas, CDs, 'LineWidth',1, 'Color', colors(Mach_index,:))
        % plot(CDs, CLs, 'LineWidth',1, 'Color', colors(Mach_index,:))
        % plot(alphas, CLs./CDs, 'LineWidth',1, 'Color', colors(Mach_index,:))
        lg{end+1} = sprintf('\n%s  $%.1f\\le\\alpha\\le %.1f$\n', NACA, alphas(1), alphas(end));
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
    for NACA_index = 1:length(NACA_list)
        NACA = NACA_list(NACA_index);
        results = fetch(sqlite(db_file), sprintf('SELECT * FROM NACA_data where NACA = %s and Mach_inf = %f and delta_t <= 1e-5;', NACA, Mach));
    
        results = sortrows(results, 'angle_of_attack_deg');
        IDs = results.ID;
        alphas = results.angle_of_attack_deg;
        CLs = results.CL;
        CDs = results.CD;
        
        % plot(alphas, CLs, 'LineWidth',1, 'Color', colors(Mach_index,:))
        plot(alphas, CDs, 'LineWidth',1, 'Color', colors(NACA_index,:))
        % plot(CDs, CLs, 'LineWidth',1, 'Color', colors(Mach_index,:))
        % plot(alphas, CLs./CDs, 'LineWidth',1, 'Color', colors(Mach_index,:))
        lg{end+1} = sprintf('\n%s  $%.1f\\le\\alpha\\le %.1f$\n', NACA, alphas(1), alphas(end));
    end
    xlabel('$\alpha$','FontSize',font_size,'Interpreter','latex')
    ylabel('$CD$','FontSize',font_size,'Interpreter','latex')
    legend(lg, 'Location','eastoutside','Interpreter','latex')
    grid on
    grid minor
    box on
subplot(2,2,3)
    hold all
    lg = {};
    for NACA_index = 1:length(NACA_list)
        NACA = NACA_list(NACA_index);
        results = fetch(sqlite(db_file), sprintf('SELECT * FROM NACA_data where NACA = %s and Mach_inf = %f and delta_t <= 1e-5;', NACA, Mach));
    
        results = sortrows(results, 'angle_of_attack_deg');
        IDs = results.ID;
        alphas = results.angle_of_attack_deg;
        CLs = results.CL;
        CDs = results.CD;
        
        % plot(alphas, CLs, 'LineWidth',1, 'Color', colors(Mach_index,:))
        % plot(alphas, CDs, 'LineWidth',1, 'Color', colors(Mach_index,:))
        plot(CDs, CLs, 'LineWidth',1, 'Color', colors(NACA_index,:))
        % plot(alphas, CLs./CDs, 'LineWidth',1, 'Color', colors(Mach_index,:))
        lg{end+1} = sprintf('\n%s  $%.1f\\le\\alpha\\le %.1f$\n', NACA, alphas(1), alphas(end));
    end
    xlabel('$CD$','FontSize',font_size,'Interpreter','latex')
    ylabel('$CL$','FontSize',font_size,'Interpreter','latex')
    legend(lg, 'Location','eastoutside','Interpreter','latex')
    grid on
    grid minor
    box on
subplot(2,2,4)
    hold all
    lg = {};
    for NACA_index = 1:length(NACA_list)
        NACA = NACA_list(NACA_index);
        results = fetch(sqlite(db_file), sprintf('SELECT * FROM NACA_data where NACA = %s and Mach_inf = %f and delta_t <= 1e-5;', NACA, Mach));
    
        results = sortrows(results, 'angle_of_attack_deg');
        IDs = results.ID;
        alphas = results.angle_of_attack_deg;
        CLs = results.CL;
        CDs = results.CD;
        
        % plot(alphas, CLs, 'LineWidth',1, 'Color', colors(Mach_index,:))
        % plot(alphas, CDs, 'LineWidth',1, 'Color', colors(Mach_index,:))
        % plot(CDs, CLs, 'LineWidth',1, 'Color', colors(Mach_index,:))
        plot(alphas, CLs./CDs, 'LineWidth',1, 'Color', colors(NACA_index,:))
        lg{end+1} = sprintf('\n%s  $%.1f\\le\\alpha\\le %.1f$\n', NACA, alphas(1), alphas(end));
    end
    xlabel('$\alpha$','FontSize',font_size,'Interpreter','latex')
    ylabel('$\frac{CL}{CD}$','FontSize',font_size,'Interpreter','latex')
    legend(lg, 'Location','eastoutside','Interpreter','latex')
    grid on
    grid minor
    box on
;








