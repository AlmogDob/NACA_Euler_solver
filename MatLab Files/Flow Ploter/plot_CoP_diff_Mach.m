clc; clear; close all;

db_file = 'NACA.db';
db = sqlread(sqlite(db_file),'NACA_data');
%%
clc; close all;

NACA = '0012';

Mach_inf_list = fetch(sqlite(db_file),sprintf('SELECT Mach_inf FROM NACA_data where NACA = %s and delta_t = 1e-5;', NACA));
Mach_inf_list = unique(Mach_inf_list.Mach_inf(:));


fig1 = figure('Name','1','Position',[300,100,1200,750]);
font_size = 16;
colors = hsv(floor(length(Mach_inf_list)*1.5))*0.9;
subplot(2,1,1)
    hold all
    lg = {};
    for Mach_index = 1:length(Mach_inf_list)
        Mach_inf = Mach_inf_list(Mach_index);
        results = fetch(sqlite(db_file), sprintf('SELECT * FROM NACA_data where NACA = %s and Mach_inf = %f and delta_t = 1e-5;', NACA, Mach_inf));
        results = sortrows(results, 'angle_of_attack_deg');
        
        IDs = results.ID(:);
        alphas = results.angle_of_attack_deg;
        CoPs = [length(IDs),2];
        for i = 1:length(IDs)
            CoPs(i,:) = get_CoP(IDs(i), results);
        end
        
        title(sprintf('NACA: %s\n', results.NACA(1)))
        
        plot(alphas, CoPs(:,1), 'LineWidth',1, 'Color', colors(Mach_index,:))

        lg{end+1} = sprintf('\n%.1f $|$ $%.1f\\le\\alpha\\le %.1f$\n', Mach_inf, alphas(1), alphas(end));
    end
    xlabel('$\alpha$','FontSize',font_size,'Interpreter','latex')
    ylabel('$CoP_x$','FontSize',font_size,'Interpreter','latex')
    legend(lg, 'Location','eastoutside','Interpreter','latex')
    grid on
    grid minor
    box on
subplot(2,1,2)
    hold all
    lg = {};
    for Mach_index = 1:length(Mach_inf_list)
        Mach_inf = Mach_inf_list(Mach_index);
        results = fetch(sqlite(db_file), sprintf('SELECT * FROM NACA_data where NACA = %s and Mach_inf = %f and delta_t = 1e-5;', NACA, Mach_inf));
        results = sortrows(results, 'angle_of_attack_deg');
        
        IDs = results.ID(:);
        alphas = results.angle_of_attack_deg;
        CoPs = [length(IDs),2];
        for i = 1:length(IDs)
            CoPs(i,:) = get_CoP(IDs(i), results);
        end
        
        title(sprintf('NACA: %s\n', results.NACA(1)))
        
        plot(alphas, CoPs(:,2), 'LineWidth',1, 'Color', colors(Mach_index,:))
        
        lg{end+1} = sprintf('\n%.1f $|$ $%.1f\\le\\alpha\\le %.1f$\n', Mach_inf, alphas(1), alphas(end));
    end
    xlabel('$\alpha$','FontSize',font_size,'Interpreter','latex')
    ylabel('$CoP_y$','FontSize',font_size,'Interpreter','latex')
    legend(lg, 'Location','eastoutside','Interpreter','latex')
    grid on
    grid minor
    box on
;






function CoP = get_CoP(ID, results)
    [ni, ~, num_points_on_airfoil, ~, ~, ~, ~, ~, ~, ~, ~, x, y, ~, ~, ~, ~, p, ~, ~, ~, ~, ~] = read_matrixes_from_DB(ID, results);
    
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

    mid_points = [(points_on_airfoil(1:end-1,1)+points_on_airfoil(2:end,1))/2,(points_on_airfoil(1:end-1,2)+points_on_airfoil(2:end,2))/2];
    
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
    mom_func(x_CoP_camber, y_CoP_camber);

end

