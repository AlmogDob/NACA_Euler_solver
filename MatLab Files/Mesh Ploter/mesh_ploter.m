clc;
clear;
close all;

ni = 100+1;
nj = 47+1;
formatSpec = '%f';
fileID = fopen("matrices/x_mat_init.txt", "r");
xs = fscanf(fileID,formatSpec);
fclose(fileID);
for i = 1:nj
    x(i,:) = xs(1+(i-1)*ni:ni+(i-1)*ni);
end
x;

fileID = fopen("matrices/y_mat_init.txt", "r");
ys = fscanf(fileID,formatSpec);
fclose(fileID);
for i = 1:nj
    y(i,:) = ys(1+(i-1)*ni:ni+(i-1)*ni);
end
y;

hold all
axis equal
plot(x(:,end), y(:,end),'-*','Color',"#7E2F8E")
plot(x(:,1), y(:,1),'-*m')
plot(x(end,:), y(end,:),'-*r')
plot(x(1,:), y(1,:),'-*k')


for i = 2:(length(x(1,:))-1)
    plot(x(:,i), y(:,i),'-','Color',"#0072BD")
end
for i = 2:(length(x(:,1))-1)
    plot(x(i,:), y(i,:),'-','Color',"#0072BD")
end

title("init")

fig1 = figure ("Name","The Mesh for psi=BC phi=BC | r = 0.01, omega = 1",'Position',[100 300 900 500]);

formatSpec = '%f';
fileID = fopen("matrices/x_mat.txt", "r");
xs = fscanf(fileID,formatSpec);
fclose(fileID);
for i = 1:nj
    x(i,:) = xs(1+(i-1)*ni:ni+(i-1)*ni);
end
x;

fileID = fopen("matrices/y_mat.txt", "r");
ys = fscanf(fileID,formatSpec);
fclose(fileID);
for i = 1:nj
    y(i,:) = ys(1+(i-1)*ni:ni+(i-1)*ni);
end
y;

hold all
axis equal

plot(x(:,end), y(:,end),'-*','Color',"#7E2F8E")
plot(x(:,1), y(:,1),'-*m')
plot(x(end,:), y(end,:),'-*r')
plot(x(1,:), y(1,:),'-*k')


for i = 2:(length(x(1,:))-1)
    plot(x(:,i), y(:,i),'-','Color',"#0072BD")
end
for i = 2:(length(x(:,1))-1)
    plot(x(i,:), y(i,:),'-','Color',"#0072BD")
end

axis equal
% xlim([-5 5])
% ylim([-4.5 4.5])
% xlim([-0.3 1.3])
% ylim([-0.4 0.4])
title("The Mesh for psi=BC phi=BC | r = 0.01, omega = 1");
ylabel("Y [-]")
xlabel("X [-]")
subtitle("Almog Dobrescu 214254252")
% legend({'The Mesh'},'FontSize',11 ,'Location','southeast')
%exportgraphics(fig1, 'grap1.png','Resolution',1200);
