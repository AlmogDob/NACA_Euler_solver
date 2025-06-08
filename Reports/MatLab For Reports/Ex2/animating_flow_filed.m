clc; clear; close all

data = cell(1);
num_of_tests = 24;
gamma = 1.4;

figure ("Name","animating flow field: Mach = [0.7...3] ",'Position',[100 100 900 700]);

for i = 1:num_of_tests
    clf
    param = readmatrix("results\" + sprintf("results%d\\ni_nj.txt", i));
    Q0 = readmatrix("results\" + sprintf("results%d\\Q0_mat.txt", i));
    Q1 = readmatrix("results\" + sprintf("results%d\\Q1_mat.txt", i));
    Q2 = readmatrix("results\" + sprintf("results%d\\Q2_mat.txt", i));
    Q3 = readmatrix("results\" + sprintf("results%d\\Q3_mat.txt", i));
    x = readmatrix("results\" + sprintf("results%d\\x_mat.txt", i));
    y = readmatrix("results\" + sprintf("results%d\\y_mat.txt", i));
    u = Q1./Q0;
    v = Q2./Q0;
    p = (gamma - 1).*(Q3 - 0.5.*Q0.*(u.^2 + v.^2));
    a = sqrt(gamma.*p./Q0);
    M = sqrt(u.^2 + v.^2)./a;
    Machs{i} = M;
    
    ni = param(1);
    nj = param(2);
    i_TEL = param(3);
    i_LE = param(4);
    i_TEU = param(5);

    axis equal
    hold all
    contourf(x, y, M, 3000, "LineStyle","none");
    colorbar;
    colormap("turbo");
    
    scale = 0.1/10000;
    q = quiver(x, y, u*scale, v*scale, "off", "Color", "#FFFFFF");
    xlim([-1.75, 2.5])
    ylim([-2, 2])
    
    plot(x(:,end), y(:,end),'-*','Color',"#7E2F8E")
    plot(x(:,1), y(:,1),'-*m')
    plot(x(end,:), y(end,:),'-*r')
    plot(x(1,:), y(1,:),'-k')

    title(sprintf("%d: Mach = %f", i, M(((ni-1)/2), (nj))))

    drawnow
    % input("press");
end


