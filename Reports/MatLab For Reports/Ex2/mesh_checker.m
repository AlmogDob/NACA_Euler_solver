clc; clear; close all;
formatSpec = '%f';

fileID = fopen("results\ni_nj.txt", "r");
nd = fscanf(fileID,formatSpec);
ni = nd(1);
nj = nd(2);

fileID = fopen("results\x_mat.txt", "r");
xs = fscanf(fileID,formatSpec);
fclose(fileID);
for i = 1:nj
   x(i,:) = xs(1+(i-1)*ni:ni+(i-1)*ni);
end
x;

fileID = fopen("results\y_mat.txt", "r");
ys = fscanf(fileID,formatSpec);
fclose(fileID);
for i = 1:nj
    y(i,:) = ys(1+(i-1)*ni:ni+(i-1)*ni);
end
y;

fig1 = figure ("Name","1",'Position',[500 100 900 500]);
hold all
axis equal

plot(x(:,end), y(:,end),'-*','Color',"#7E2F8E")
plot(x(:,1), y(:,1),'-*m')
plot(x(end,:), y(end,:),'-*r')
plot(x(1,:), y(1,:),'-k')

for i = 2:(length(x(1,:))-1)
    plot(x(:,i), y(:,i),'-','Color',"#0072BD")
end
for i = 2:(length(x(:,1))-1)
    plot(x(i,:), y(i,:),'-','Color',"#0072BD")
end