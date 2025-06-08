%% Convergence history mach 1.5
clc; clear; close all

data = readmatrix("results\\iterations.txt");

fig3 = figure ("Name","Convergence History for Mach = 1.5, dt = 3.8e-4",'Position',[100 300 900 500]);

loglog(data(:,2))

title("Convergence History for Mach = 1.5, dt = 3.8e-4");
ylabel("L2Norm [-]")
xlabel("num of iterations [-]")
subtitle("Almog Dobrescu 214254252")
% exportgraphics(fig3, 'grap8.1.png','Resolution',1200);

%% Mach distribution mach 1.5
clc; clear; close all;

gamma = 1.4;

fig4 = figure ("Name","Mach Distribution On the Airfoil for Mach = 1.5, dt = 3.8e-4",'Position',[250 300 900 500]);

param = readmatrix("results\\ni_nj.txt");
Q0 = readmatrix("results\\Q0_mat.txt");
Q1 = readmatrix("results\\Q1_mat.txt");
Q2 = readmatrix("results\\Q2_mat.txt");
Q3 = readmatrix("results\\Q3_mat.txt");

u = Q1./Q0;
v = Q2./Q0;
p = (gamma - 1).*(Q3 - 0.5.*Q0.*(u.^2 + v.^2));
a = sqrt(gamma.*p./Q0);
M = sqrt(u.^2 + v.^2)./a;
i_TEL = param(3);
i_LE = param(4);
i_TEU = param(5);

plot(i_TEL:1:i_TEU, M(1,i_TEL+1:1:i_TEU+1))

title("Mach Distribution On the Airfoil for Mach 1.5, dt = 3.8e-4");
ylabel("Mach number [-]")
xlabel("xi index [-]")
subtitle("Almog Dobrescu 214254252")
% exportgraphics(fig4, 'grap9.png','Resolution',1200);

%% Pressure distribution mach 1.5
clc; clear;

gamma = 1.4;

fig5 = figure ("Name","Pressure Distribution On the Airfoil for Mach = 1.5, dt = 3.8e-4",'Position',[250 300 900 500]);

param = readmatrix("results\\ni_nj.txt");
Q0 = readmatrix("results\\Q0_mat.txt");
Q1 = readmatrix("results\\Q1_mat.txt");
Q2 = readmatrix("results\\Q2_mat.txt");
Q3 = readmatrix("results\\Q3_mat.txt");

u = Q1./Q0;
v = Q2./Q0;
p = (gamma - 1).*(Q3 - 0.5.*Q0.*(u.^2 + v.^2));

i_TEL = param(3);
i_LE = param(4);
i_TEU = param(5);

plot(i_TEL:1:i_TEU, p(1,i_TEL+1:1:i_TEU+1))

title("Pressure Distribution On the Airfoil for Mach = 1.5, dt = 3.8e-4");
ylabel("Pressure [Pa]")
xlabel("xi index [-]")
subtitle("Almog Dobrescu 214254252")
% exportgraphics(fig5, 'grap10.png','Resolution',1200);

%% Effect of Time Step on Convergence History for Mach = 0.9
clc; clear; close all

data = cell(1);
num_of_test = 20;
times = logspace(-5, -2.5, 20);
colors_temp = parula(num_of_test+5);
% colors = colors_temp(1:num_of_test,:);
colors = colors_temp;

fig6 = figure ("Name","Effect of Time Step on Convergence History for Mach = 0.9",'Position',[400 300 900 500]);

for i = 1:num_of_test
    data{i} = readmatrix("results\" + sprintf("results%d\\iterations.txt", i));
end
leg = cell(1,length(data));
figure(1)
for i = length(data):-1:1
    % loglog(data{i}(:,2), "LineWidth",1.5, 'Color', [(255*i/length(data))/255,(100*i/length(data)+55)/255,1])
    loglog(data{i}(:,2), "LineWidth",1.5, 'Color', colors(i,:))
    hold on
    leg{length(data)-i+1} = sprintf("%d",times(i));
end
legend(leg,'FontSize',13 ,'Location','southwest')
title("Effect of Time Step on Convergence History for Mach = 0.9");
ylabel("L2Norm [-]")
xlabel("num of iterations [-]")
subtitle("Almog Dobrescu 214254252");
% exportgraphics(fig6, 'grap11.png','Resolution',1200);

%% Effect of Time Step on Mach Distribution for Mach = 0.9
clc; clear; close all

data = cell(1);
num_of_test = 20;
times = logspace(-5, -2.5, 20);
gamma = 1.4;
line_types = {"-", "--"};

fig7 = figure ("Name","Effect of Time Step on Mach Distribution for Mach = 0.9",'Position',[250 300 900 500]);

title("Effect of Time Step on Mach Distribution for Mach = 0.9");
ylabel("Mach number [-]")
xlabel("xi index [-]")
subtitle("Almog Dobrescu 214254252")
% exportgraphics(fig7, 'grap13.png','Resolution',1200);

for i = 1:num_of_test
    Q0 = readmatrix("results\" + sprintf("results%d\\Q0_mat.txt", i));
    Q1 = readmatrix("results\" + sprintf("results%d\\Q1_mat.txt", i));
    Q2 = readmatrix("results\" + sprintf("results%d\\Q2_mat.txt", i));
    Q3 = readmatrix("results\" + sprintf("results%d\\Q3_mat.txt", i));
    u = Q1./Q0;
    v = Q2./Q0;
    p = (gamma - 1).*(Q3 - 0.5.*Q0.*(u.^2 + v.^2));
    a = sqrt(gamma.*p./Q0);
    M = sqrt(u.^2 + v.^2)./a;
    Machs{i} = M;
end
hold all
j = 1;
for i = num_of_test:-19:1
    param = readmatrix("results\\ni_nj.txt");
    i_TEL = param(3);
    i_LE = param(4);
    i_TEU = param(5);
    plot(i_TEL:1:i_TEU, Machs{i}(1,i_TEL+1:1:i_TEU+1), line_types{j}, "LineWidth",2, 'Color', [(255*i/num_of_test)/255,(100*i/num_of_test+55)/255,1])
    % leg{num_of_test-i+1} = sprintf("%d",times(i));
    j = j + 1;
end
% legend(leg,'FontSize',13 ,'Location','southwest')
legend({sprintf("%d",times(20)), sprintf("%d",times(1))},'FontSize',13 ,'Location','southwest')

%% Effect of Time Step on Pressure Distribution for Mach = 0.9
clc; clear; close all

data = cell(1);
num_of_test = 20;
times = logspace(-5, -2.5, 20);
gamma = 1.4;
line_types = {"-", "--"};

fig8 = figure ("Name","Effect of Time Step on Pressure Distribution for Mach = 0.9",'Position',[250 300 900 500]);

title("Effect of Time Step on Pressure Distribution for Mach = 0.9");
ylabel("Pressure [Pa]")
xlabel("xi index [-]")
subtitle("Almog Dobrescu 214254252")
% exportgraphics(fig8, 'grap14.png','Resolution',1200);

for i = 1:num_of_test
    Q0 = readmatrix("results\" + sprintf("results%d\\Q0_mat.txt", i));
    Q1 = readmatrix("results\" + sprintf("results%d\\Q1_mat.txt", i));
    Q2 = readmatrix("results\" + sprintf("results%d\\Q2_mat.txt", i));
    Q3 = readmatrix("results\" + sprintf("results%d\\Q3_mat.txt", i));
    u = Q1./Q0;
    v = Q2./Q0;
    p = (gamma - 1).*(Q3 - 0.5.*Q0.*(u.^2 + v.^2));
    Ps{i} = p;
end
hold all
j = 1;
for i = num_of_test:-19:1
    param = readmatrix("results\\ni_nj.txt");
    i_TEL = param(3);
    i_LE = param(4);
    i_TEU = param(5);
    plot(i_TEL:1:i_TEU, Ps{i}(1,i_TEL+1:1:i_TEU+1), line_types{j}, "LineWidth",2, 'Color', [(255*i/num_of_test)/255,(100*i/num_of_test+55)/255,1])
    % leg{num_of_test-i+1} = sprintf("%d",times(i));
    j = j + 1;
end
% legend(leg,'FontSize',13 ,'Location','southwest')
legend({sprintf("%d",times(20)), sprintf("%d",times(1))},'FontSize',13 ,'Location','northwest')