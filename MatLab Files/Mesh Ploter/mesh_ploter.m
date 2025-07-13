clc;
clear;
close all;

NACA = 2360;
result = get_result(NACA);

x = result.x_mat_init;
y = result.y_mat_init;

hold all
axis equal
plot(x(:,end), y(:,end),'-','Color',"#7E2F8E")
plot(x(:,1), y(:,1),'-m')
plot(x(end,:), y(end,:),'-r')
plot(x(1,:), y(1,:),'-k')


for i = 2:(length(x(1,:))-1)
    plot(x(:,i), y(:,i),'-','Color',"#0072BD")
end
for i = 2:(length(x(:,1))-1)
    plot(x(i,:), y(i,:),'-','Color',"#0072BD")
end

grid on
grid minor
title("The Initial Mesh");




fig1 = figure ('Name', '1','Position',[100 300 900 500]);

x = result.x_mat;
y = result.y_mat;

hold all
axis equal

plot(x(:,end), y(:,end),'-','Color',"#7E2F8E")
plot(x(:,1), y(:,1),'-m')
plot(x(end,:), y(end,:),'-r')
plot(x(1,:), y(1,:),'-k')


for i = 2:(length(x(1,:))-1)
    plot(x(:,i), y(:,i),'-','Color',"#0072BD")
end
for i = 2:(length(x(:,1))-1)
    plot(x(i,:), y(i,:),'-','Color',"#0072BD")
end

axis equal
grid on
grid minor
title("The Mesh");
ylabel("Y [-]")
xlabel("X [-]")
% legend({'The Mesh'},'FontSize',11 ,'Location','southeast')
%exportgraphics(fig1, 'grap1.png','Resolution',1200);
