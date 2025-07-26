clc; clear; close all;

db_file = 'NACA.db';
db = sqlread(sqlite(db_file),'NACA_data');

results = fetch(sqlite(db_file),'SELECT * FROM NACA_data where NACA = 5410 and Mach_inf = 0.7 and delta_t = 1e-5 and angle_of_attack_deg = 10;');
results = sortrows(results, 'angle_of_attack_deg');
ID = results.ID;

[ni, nj, num_points_on_airfoil, gamma, rho_inf, Mach_inf, p_inf, a_inf, u_inf, alpha_deg, NACA, x, y, rho, u, v, e, p, a, M, p0, CL, CD] = read_matrixes_from_DB(ID, results);

figure
hold all
contourf(x, y, M, 300, "LineStyle","none");
colorbar;
colormap("turbo");

factor = 0.0001;
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
% figure
% for i = 1:(15)
%     i_index_from_free_end = i;
%     plot([flip(p0(:,i_index_from_free_end));p0(:,ni-i_index_from_free_end)])
%     hold on
% end
% ylabel('p_0')

figure
% for j = 1:(nj-6)
for j = 25
    j_index_from_free_end = j;
    plot(x(j_index_from_free_end,:), p0(j_index_from_free_end,:))
    % plot(p0(j_index_from_free_end,:))
    hold on
end
ylabel('p_0')

figure
% for j = 1:(nj-6)
for j = 1
    j_index_from_free_end = j;
    plot(x(j_index_from_free_end,:), first_order_first_derive(x(j_index_from_free_end,:), p0(j_index_from_free_end,:)))
    % plot(first_order_first_derive(x(j_index_from_free_end,:), p0(j_index_from_free_end,:)))
    hold on
end
ylabel('$\displaystyle\frac{\partial p_0}{\partial x}$','FontSize',15, 'Interpreter','latex')
% xlim([0.1,1])





% FUNCTIONS ===============================================================
function dy_dx = first_order_first_derive(x, y)
    if length(x) ~= length(y)   
        dy_dx = nan;
        return
    end
    dy_dx = zeros(length(x),1);
    for i = 1:length(x)-1
        dy_dx(i) = (y(i+1) - y(i)) / (x(i+1) - x(i));
    end
    dy_dx(end) = (y(end) - y(end-1)) / (x(end) - x(end-1));

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
