clc; clear; close all;

db_file = 'NACA.db';
db = sqlread(sqlite(db_file),'NACA_data');

results = fetch(sqlite(db_file),'SELECT * FROM NACA_data where NACA = 4410 and Mach_inf = 0.7 and delta_t = 1e-5 and angle_of_attack_deg = 5;');
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

