clc; clear; close all;
ni = 51;
nj = 26;
gamma = 1.4;

formatSpec = '%f';

fileID = fopen("results/Q0_mat.txt", "r");
Q0s = fscanf(fileID,formatSpec);
fclose(fileID);
for i = 1:nj
    Q0(i,:) = Q0s(1+(i-1)*ni:ni+(i-1)*ni);
end
Q0;

fileID = fopen("results/Q1_mat.txt", "r");
Q1s = fscanf(fileID,formatSpec);
fclose(fileID);
for i = 1:nj
    Q1(i,:) = Q1s(1+(i-1)*ni:ni+(i-1)*ni);
end
Q1;

fileID = fopen("results/Q2_mat.txt", "r");
Q2s = fscanf(fileID,formatSpec);
fclose(fileID);
for i = 1:nj
    Q2(i,:) = Q2s(1+(i-1)*ni:ni+(i-1)*ni);
end
Q2;

fileID = fopen("results/Q3_mat.txt", "r");
Q3s = fscanf(fileID,formatSpec);
fclose(fileID);
for i = 1:nj
    Q3(i,:) = Q3s(1+(i-1)*ni:ni+(i-1)*ni);
end
Q3;

fileID = fopen("results/x_mat.txt", "r");
xs = fscanf(fileID,formatSpec);
fclose(fileID);
for i = 1:nj
   x(i,:) = xs(1+(i-1)*ni:ni+(i-1)*ni);
end
x;

fileID = fopen("results/y_mat.txt", "r");
ys = fscanf(fileID,formatSpec);
fclose(fileID);
for i = 1:nj
    y(i,:) = ys(1+(i-1)*ni:ni+(i-1)*ni);
end
y;

u = Q1./Q0;
v = Q2./Q0;
p = (gamma - 1).*(Q3 - 0.5.*Q0.*(u.^2 + v.^2));
a = sqrt(gamma.*p./Q0);
M = sqrt(u.^2 + v.^2)./a;

hold all
contourf(x, y, M, 2000, "LineStyle","none");
colorbar;
colormap("jet");

% plot(x(:,end), y(:,end),'-*','Color',"#7E2F8E")
% plot(x(:,1), y(:,1),'-*m')
% plot(x(end,:), y(end,:),'-*r')
% plot(x(1,:), y(1,:),'-k')
% 
% for i = 2:(length(x(1,:))-1)
%     plot(x(:,i), y(:,i),'-','Color',"#0072BD")
% end
% for i = 2:(length(x(:,1))-1)
%     plot(x(i,:), y(i,:),'-','Color',"#0072BD")
% end