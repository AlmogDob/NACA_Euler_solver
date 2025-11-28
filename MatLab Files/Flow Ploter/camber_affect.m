clc; clear; close all;

db_file = 'NACA.db';
db_table = sqlread(sqlite(db_file),'NACA_data');

NACA_second_digit = '4';
NACA_thickness    = '12';
Mach = 0.9;

NACA_4_digit = strcat(NACA_second_digit, NACA_thickness)
NACA_list_4_digit = unique(fetch(sqlite(db_file),sprintf('SELECT NACA FROM NACA_data where NACA GLOB "?%s" and Mach_inf = %f;', NACA_4_digit, Mach)));
NACA_list_4_digit = NACA_list_4_digit.NACA(:);

NACA_5_digit_not_reflexed = strcat(NACA_second_digit, '0', NACA_thickness)
NACA_list_5_digit_not_reflexed = unique(fetch(sqlite(db_file),sprintf('SELECT NACA FROM NACA_data where NACA GLOB "?%s" and Mach_inf = %f;', NACA_5_digit_not_reflexed, Mach)));
NACA_list_5_digit_not_reflexed = NACA_list_5_digit_not_reflexed.NACA(:);

NACA_5_digit_reflexed = strcat(NACA_second_digit, '1', NACA_thickness)
NACA_list_5_digit_reflexed = unique(fetch(sqlite(db_file),sprintf('SELECT NACA FROM NACA_data where NACA GLOB "?%s" and Mach_inf = %f;', NACA_5_digit_reflexed, Mach)));
NACA_list_5_digit_reflexed = NACA_list_5_digit_reflexed.NACA(:);

% NACA_list = [0010;1410;2410;3410;4410;5410;6410;7410;8410;9410]

fig1 = figure('Name','1 CL_CD','Position',[200,100,1400,750]);
font_size = 16;
sgtitle(sprintf('Mach: %.1f\nNACA: *%s $|$ *%s $|$ *%s',Mach, NACA_4_digit, NACA_5_digit_not_reflexed, NACA_5_digit_reflexed),'FontSize',font_size,'Interpreter','latex');
colors = hsv(length(NACA_list_4_digit))*0.9;
subplot(3,5,1)
    hold all
    lg = {};
    % Mach_inf_list = 0.6;
    for NACA_index = 1:length(NACA_list_4_digit)
        NACA = NACA_list_4_digit(NACA_index);
        results = fetch(sqlite(db_file), sprintf('SELECT * FROM NACA_data where NACA = %s and Mach_inf = %f and delta_t <= 1e-5;', NACA, Mach));
    
        results = sortrows(results, 'angle_of_attack_deg');
        IDs = results.ID;
        alphas = results.angle_of_attack_deg;
        CLs = results.CL;
        CDs = results.CD;
        
        plot(alphas, CLs, 'LineWidth',1, 'Color', colors(NACA_index,:))
        lg{end+1} = sprintf('\n%s $%.1f\\le\\alpha\\le %.1f$\n', NACA, alphas(1), alphas(end));
    end
    xlabel('$\alpha$','FontSize',font_size,'Interpreter','latex')
    ylabel('$CL$','FontSize',font_size,'Interpreter','latex')
    % legend(lg, 'Location','eastoutside','FontSize',font_size-6,'Interpreter','latex')
    grid on
    grid minor
    box on
subplot(3,5,2)
    hold all
    lg = {};
    for NACA_index = 1:length(NACA_list_4_digit)
        NACA = NACA_list_4_digit(NACA_index);
        results = fetch(sqlite(db_file), sprintf('SELECT * FROM NACA_data where NACA = %s and Mach_inf = %f and delta_t <= 1e-5;', NACA, Mach));
    
        results = sortrows(results, 'angle_of_attack_deg');
        IDs = results.ID;
        alphas = results.angle_of_attack_deg;
        CLs = results.CL;
        CDs = results.CD;
        
        plot(alphas, CDs, 'LineWidth',1, 'Color', colors(NACA_index,:))
        lg{end+1} = sprintf('\n%s $%.1f\\le\\alpha\\le %.1f$\n', NACA, alphas(1), alphas(end));
    end
    xlabel('$\alpha$','FontSize',font_size,'Interpreter','latex')
    ylabel('$CD$','FontSize',font_size,'Interpreter','latex')
    % legend(lg, 'Location','eastoutside','FontSize',font_size-6,'Interpreter','latex')
    grid on
    grid minor
    box on
subplot(3,5,3)
    hold all
    lg = {};
    for NACA_index = 1:length(NACA_list_4_digit)
        NACA = NACA_list_4_digit(NACA_index);
        results = fetch(sqlite(db_file), sprintf('SELECT * FROM NACA_data where NACA = %s and Mach_inf = %f and delta_t <= 1e-5;', NACA, Mach));
    
        results = sortrows(results, 'angle_of_attack_deg');
        IDs = results.ID;
        alphas = results.angle_of_attack_deg;
        CLs = results.CL;
        CDs = results.CD;
        
        plot(alphas, CLs./CDs, 'LineWidth',1, 'Color', colors(NACA_index,:))
        lg{end+1} = sprintf('\n%s $%.1f\\le\\alpha\\le %.1f$\n', NACA, alphas(1), alphas(end));
    end
    xlabel('$\alpha$','FontSize',font_size,'Interpreter','latex')
    ylabel('$\frac{CL}{CD}$','FontSize',font_size,'Interpreter','latex')
    % legend(lg, 'Location','eastoutside','FontSize',font_size-6,'Interpreter','latex')
    grid on
    grid minor
    box on
subplot(3,5,4)
    hold all
    lg = {};
    for NACA_index = 1:length(NACA_list_4_digit)
        NACA = NACA_list_4_digit(NACA_index);
        results = fetch(sqlite(db_file), sprintf('SELECT * FROM NACA_data where NACA = %s and Mach_inf = %f and delta_t <= 1e-5;', NACA, Mach));
    
        results = sortrows(results, 'angle_of_attack_deg');
        IDs = results.ID;
        alphas = results.angle_of_attack_deg;
        CLs = results.CL;
        CDs = results.CD;
        
        plot(CDs, CLs, 'LineWidth',1, 'Color', colors(NACA_index,:))
        lg{end+1} = sprintf('\n%s $%.1f\\le\\alpha\\le %.1f$\n', NACA, alphas(1), alphas(end));
    end
    xlabel('$CD$','FontSize',font_size,'Interpreter','latex')
    ylabel('$CL$','FontSize',font_size,'Interpreter','latex')
    % legend(lg, 'Location','eastoutside','FontSize',font_size-6,'Interpreter','latex')
    grid on
    grid minor
    box on
subplot(3,5,5)
    hold all
    lg = {};
    for NACA_index = 1:length(NACA_list_4_digit)
        NACA = NACA_list_4_digit(NACA_index);
        results = fetch(sqlite(db_file), sprintf('SELECT * FROM NACA_data where NACA = %s and Mach_inf = %f and delta_t <= 1e-5;', NACA, Mach));
    
        results = sortrows(results, 'angle_of_attack_deg');
        IDs = results.ID;
        alphas = results.angle_of_attack_deg;
        CLs = results.CL;
        CDs = results.CD;
        
        plot(NaN, NaN, 'LineWidth',1, 'Color', colors(NACA_index,:))
        lg{end+1} = sprintf('\n%s $%.1f\\le\\alpha\\le %.1f$\n', NACA, alphas(1), alphas(end));
    end
    % xlabel('$CD$','FontSize',font_size,'Interpreter','latex')
    % ylabel('$CL$','FontSize',font_size,'Interpreter','latex')
    legend(lg, 'Location','eastoutside','FontSize',font_size-6,'Interpreter','latex')
    grid on
    grid minor
    box on
    % ---------------------------------------------------------------------
colors = hsv(length(NACA_list_5_digit_not_reflexed))*0.9;
subplot(3,5,6)
    hold all
    lg = {};
    % Mach_inf_list = 0.6;
    for NACA_index = 1:length(NACA_list_5_digit_not_reflexed)
        NACA = NACA_list_5_digit_not_reflexed(NACA_index);
        results = fetch(sqlite(db_file), sprintf('SELECT * FROM NACA_data where NACA = %s and Mach_inf = %f and delta_t <= 1e-5;', NACA, Mach));
    
        results = sortrows(results, 'angle_of_attack_deg');
        IDs = results.ID;
        alphas = results.angle_of_attack_deg;
        CLs = results.CL;
        CDs = results.CD;
        
        plot(alphas, CLs, 'LineWidth',1, 'Color', colors(NACA_index,:))
        lg{end+1} = sprintf('\n%s $%.1f\\le\\alpha\\le %.1f$\n', NACA, alphas(1), alphas(end));
    end
    xlabel('$\alpha$','FontSize',font_size,'Interpreter','latex')
    ylabel('$CL$','FontSize',font_size,'Interpreter','latex')
    % legend(lg, 'Location','eastoutside','FontSize',font_size-6,'Interpreter','latex')
    grid on
    grid minor
    box on
subplot(3,5,7)
    hold all
    lg = {};
    for NACA_index = 1:length(NACA_list_5_digit_not_reflexed)
        NACA = NACA_list_5_digit_not_reflexed(NACA_index);
        results = fetch(sqlite(db_file), sprintf('SELECT * FROM NACA_data where NACA = %s and Mach_inf = %f and delta_t <= 1e-5;', NACA, Mach));
    
        results = sortrows(results, 'angle_of_attack_deg');
        IDs = results.ID;
        alphas = results.angle_of_attack_deg;
        CLs = results.CL;
        CDs = results.CD;
        
        plot(alphas, CDs, 'LineWidth',1, 'Color', colors(NACA_index,:))
        lg{end+1} = sprintf('\n%s $%.1f\\le\\alpha\\le %.1f$\n', NACA, alphas(1), alphas(end));
    end
    xlabel('$\alpha$','FontSize',font_size,'Interpreter','latex')
    ylabel('$CD$','FontSize',font_size,'Interpreter','latex')
    % legend(lg, 'Location','eastoutside','FontSize',font_size-6,'Interpreter','latex')
    grid on
    grid minor
    box on
subplot(3,5,8)
    hold all
    lg = {};
    for NACA_index = 1:length(NACA_list_5_digit_not_reflexed)
        NACA = NACA_list_5_digit_not_reflexed(NACA_index);
        results = fetch(sqlite(db_file), sprintf('SELECT * FROM NACA_data where NACA = %s and Mach_inf = %f and delta_t <= 1e-5;', NACA, Mach));
    
        results = sortrows(results, 'angle_of_attack_deg');
        IDs = results.ID;
        alphas = results.angle_of_attack_deg;
        CLs = results.CL;
        CDs = results.CD;
        
        plot(alphas, CLs./CDs, 'LineWidth',1, 'Color', colors(NACA_index,:))
        lg{end+1} = sprintf('\n%s $%.1f\\le\\alpha\\le %.1f$\n', NACA, alphas(1), alphas(end));
    end
    xlabel('$\alpha$','FontSize',font_size,'Interpreter','latex')
    ylabel('$\frac{CL}{CD}$','FontSize',font_size,'Interpreter','latex')
    % legend(lg, 'Location','eastoutside','FontSize',font_size-6,'Interpreter','latex')
    grid on
    grid minor
    box on
subplot(3,5,9)
    hold all
    lg = {};
    for NACA_index = 1:length(NACA_list_5_digit_not_reflexed)
        NACA = NACA_list_5_digit_not_reflexed(NACA_index);
        results = fetch(sqlite(db_file), sprintf('SELECT * FROM NACA_data where NACA = %s and Mach_inf = %f and delta_t <= 1e-5;', NACA, Mach));
    
        results = sortrows(results, 'angle_of_attack_deg');
        IDs = results.ID;
        alphas = results.angle_of_attack_deg;
        CLs = results.CL;
        CDs = results.CD;
        
        plot(CDs, CLs, 'LineWidth',1, 'Color', colors(NACA_index,:))
        lg{end+1} = sprintf('\n%s $%.1f\\le\\alpha\\le %.1f$\n', NACA, alphas(1), alphas(end));
    end
    xlabel('$CD$','FontSize',font_size,'Interpreter','latex')
    ylabel('$CL$','FontSize',font_size,'Interpreter','latex')
    % legend(lg, 'Location','eastoutside','FontSize',font_size-6,'Interpreter','latex')
    grid on
    grid minor
    box on
subplot(3,5,10)
    hold all
    lg = {};
    for NACA_index = 1:length(NACA_list_5_digit_not_reflexed)
        NACA = NACA_list_5_digit_not_reflexed(NACA_index);
        results = fetch(sqlite(db_file), sprintf('SELECT * FROM NACA_data where NACA = %s and Mach_inf = %f and delta_t <= 1e-5;', NACA, Mach));
    
        results = sortrows(results, 'angle_of_attack_deg');
        IDs = results.ID;
        alphas = results.angle_of_attack_deg;
        CLs = results.CL;
        CDs = results.CD;
        
        plot(NaN, NaN, 'LineWidth',1, 'Color', colors(NACA_index,:))
        lg{end+1} = sprintf('\n%s $%.1f\\le\\alpha\\le %.1f$\n', NACA, alphas(1), alphas(end));
    end
    % xlabel('$CD$','FontSize',font_size,'Interpreter','latex')
    % ylabel('$CL$','FontSize',font_size,'Interpreter','latex')
    legend(lg, 'Location','eastoutside','FontSize',font_size-6,'Interpreter','latex')
    grid on
    grid minor
    box on
    % ---------------------------------------------------------------------
colors = hsv(length(NACA_list_5_digit_reflexed))*0.9;
subplot(3,5,11)
    hold all
    lg = {};
    % Mach_inf_list = 0.6;
    for NACA_index = 1:length(NACA_list_5_digit_reflexed)
        NACA = NACA_list_5_digit_reflexed(NACA_index);
        results = fetch(sqlite(db_file), sprintf('SELECT * FROM NACA_data where NACA = %s and Mach_inf = %f and delta_t <= 1e-5;', NACA, Mach));
    
        results = sortrows(results, 'angle_of_attack_deg');
        IDs = results.ID;
        alphas = results.angle_of_attack_deg;
        CLs = results.CL;
        CDs = results.CD;
        
        plot(alphas, CLs, 'LineWidth',1, 'Color', colors(NACA_index,:))
        lg{end+1} = sprintf('\n%s $%.1f\\le\\alpha\\le %.1f$\n', NACA, alphas(1), alphas(end));
    end
    xlabel('$\alpha$','FontSize',font_size,'Interpreter','latex')
    ylabel('$CL$','FontSize',font_size,'Interpreter','latex')
    % legend(lg, 'Location','eastoutside','FontSize',font_size-6,'Interpreter','latex')
    grid on
    grid minor
    box on
subplot(3,5,12)
    hold all
    lg = {};
    for NACA_index = 1:length(NACA_list_5_digit_reflexed)
        NACA = NACA_list_5_digit_reflexed(NACA_index);
        results = fetch(sqlite(db_file), sprintf('SELECT * FROM NACA_data where NACA = %s and Mach_inf = %f and delta_t <= 1e-5;', NACA, Mach));
    
        results = sortrows(results, 'angle_of_attack_deg');
        IDs = results.ID;
        alphas = results.angle_of_attack_deg;
        CLs = results.CL;
        CDs = results.CD;
        
        plot(alphas, CDs, 'LineWidth',1, 'Color', colors(NACA_index,:))
        lg{end+1} = sprintf('\n%s $%.1f\\le\\alpha\\le %.1f$\n', NACA, alphas(1), alphas(end));
    end
    xlabel('$\alpha$','FontSize',font_size,'Interpreter','latex')
    ylabel('$CD$','FontSize',font_size,'Interpreter','latex')
    % legend(lg, 'Location','eastoutside','FontSize',font_size-6,'Interpreter','latex')
    grid on
    grid minor
    box on
subplot(3,5,13)
    hold all
    lg = {};
    for NACA_index = 1:length(NACA_list_5_digit_reflexed)
        NACA = NACA_list_5_digit_reflexed(NACA_index);
        results = fetch(sqlite(db_file), sprintf('SELECT * FROM NACA_data where NACA = %s and Mach_inf = %f and delta_t <= 1e-5;', NACA, Mach));
    
        results = sortrows(results, 'angle_of_attack_deg');
        IDs = results.ID;
        alphas = results.angle_of_attack_deg;
        CLs = results.CL;
        CDs = results.CD;
        
        plot(alphas, CLs./CDs, 'LineWidth',1, 'Color', colors(NACA_index,:))
        lg{end+1} = sprintf('\n%s $%.1f\\le\\alpha\\le %.1f$\n', NACA, alphas(1), alphas(end));
    end
    xlabel('$\alpha$','FontSize',font_size,'Interpreter','latex')
    ylabel('$\frac{CL}{CD}$','FontSize',font_size,'Interpreter','latex')
    % legend(lg, 'Location','eastoutside','FontSize',font_size-6,'Interpreter','latex')
    grid on
    grid minor
    box on
subplot(3,5,14)
    hold all
    lg = {};
    for NACA_index = 1:length(NACA_list_5_digit_reflexed)
        NACA = NACA_list_5_digit_reflexed(NACA_index);
        results = fetch(sqlite(db_file), sprintf('SELECT * FROM NACA_data where NACA = %s and Mach_inf = %f and delta_t <= 1e-5;', NACA, Mach));
    
        results = sortrows(results, 'angle_of_attack_deg');
        IDs = results.ID;
        alphas = results.angle_of_attack_deg;
        CLs = results.CL;
        CDs = results.CD;
        
        plot(CDs, CLs, 'LineWidth',1, 'Color', colors(NACA_index,:))
        lg{end+1} = sprintf('\n%s $%.1f\\le\\alpha\\le %.1f$\n', NACA, alphas(1), alphas(end));
    end
    xlabel('$CD$','FontSize',font_size,'Interpreter','latex')
    ylabel('$CL$','FontSize',font_size,'Interpreter','latex')
    % legend(lg, 'Location','eastoutside','FontSize',font_size-6,'Interpreter','latex')
    grid on
    grid minor
    box on
subplot(3,5,15)
    hold all
    lg = {};
    for NACA_index = 1:length(NACA_list_5_digit_reflexed)
        NACA = NACA_list_5_digit_reflexed(NACA_index);
        results = fetch(sqlite(db_file), sprintf('SELECT * FROM NACA_data where NACA = %s and Mach_inf = %f and delta_t <= 1e-5;', NACA, Mach));
    
        results = sortrows(results, 'angle_of_attack_deg');
        IDs = results.ID;
        alphas = results.angle_of_attack_deg;
        CLs = results.CL;
        CDs = results.CD;
        
        plot(NaN, NaN, 'LineWidth',1, 'Color', colors(NACA_index,:))
        lg{end+1} = sprintf('\n%s $%.1f\\le\\alpha\\le %.1f$\n', NACA, alphas(1), alphas(end));
    end
    % xlabel('$CD$','FontSize',font_size,'Interpreter','latex')
    % ylabel('$CL$','FontSize',font_size,'Interpreter','latex')
    legend(lg, 'Location','eastoutside','FontSize',font_size-6,'Interpreter','latex')
    grid on
    grid minor
    box on
;









