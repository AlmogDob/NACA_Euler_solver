clc;
clear;
close all;

formatSpec = '%f';
fileID = fopen("matrices/x_mat_init.txt", "r");
xs = fscanf(fileID,formatSpec);
fclose(fileID);
for i = 1:26
    x(i,:) = xs(1+(i-1)*51:51+(i-1)*51);
end
x;

fileID = fopen("matrices/y_mat_init.txt", "r");
ys = fscanf(fileID,formatSpec);
fclose(fileID);
for i = 1:26
    y(i,:) = ys(1+(i-1)*51:51+(i-1)*51);
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

%%
fig1 = figure ("Name","The Mesh for psi=BC phi=BC | r = 0.01, omega = 1",'Position',[100 300 900 500]);

formatSpec = '%f';
for index = 2410
    index
    clf
    dirX = sprintf("matrices/x_mat_%d.txt", index);
    dirY = sprintf("matrices/y_mat_%d.txt", index);
    fileID = fopen(dirX, "r");
    xs = fscanf(fileID,formatSpec);
    fclose(fileID);
    for i = 1:26
        x(i,:) = xs(1+(i-1)*51:51+(i-1)*51);
    end
    x;
    
    fileID = fopen(dirY, "r");
    ys = fscanf(fileID,formatSpec);
    fclose(fileID);
    for i = 1:26
        y(i,:) = ys(1+(i-1)*51:51+(i-1)*51);
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
    % drawnow
    % input("press");
end
% axis equal
xlim([-5 5])
ylim([-4.5 4.5])
% xlim([-0.3 1.3])
% ylim([-0.4 0.4])
title("The Mesh for psi=BC phi=BC | r = 0.01, omega = 1");
ylabel("Y [-]")
xlabel("X [-]")
% legend({'The Mesh'},'FontSize',11 ,'Location','southeast')
%exportgraphics(fig1, 'grap1.png','Resolution',1200);

%% 
Ls = readtable("matrices/Ls_valuse");
Lx = Ls.Var1;
Ly = Ls.Var2;

fig2 = figure ("Name","Residual Convergens for psi=BC phi=BC | r = 0.01, omega = 1",'Position',[100 300 900 500]);

semilogy(Lx,'LineWidth',2)
hold on
semilogy(Ly,'LineWidth',2)

grid on 
grid minor
title("Residual Convergens for phi=BC psi=BC | r = 0.01, omega = 1");
ylabel("L [-]")
xlabel("Number of Iterations [-]")
legend({'Lx', 'Ly'},'FontSize',11 ,'Location','northeast')
% exportgraphics(fig2, 'grap1.2.png','Resolution',1200);

%%

r = logspace(-3, 0,20);
omega = linspace(1e-2, 2-1e-3, 20);

omega_fixed = readtable("omeg 1 different r.txt");
omega_fixed = omega_fixed.Var1;

r_fixed = readtable("r 0.01 different omeg.txt");
r_fixed = r_fixed.Var1;

fig3 = figure ("Name","Number of Iterations for fixed r and different omega| psi=BC phi=BC",'Position',[250 300 900 500]);

loglog(omega, r_fixed,'LineWidth',2)

grid on 
grid minor
title("Number of Iterations for fixed r and different omega| psi=BC phi=BC");
ylabel("Number of Iterations [-]")
xlabel("omega [-]")
% legend({''},'FontSize',11 ,'Location','northeast')
% exportgraphics(fig3, 'grap3.1.png','Resolution',1200);

