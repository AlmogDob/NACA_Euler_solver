clc; clear;

db_file = 'NACA.db';
db = sqlread(sqlite(db_file),'NACA_data');

wanted_NACA = 0012;
ID = find(db.NACA == wanted_NACA);
% ID = 5;
if isempty(ID)
    error('No NACA=%d in Data Base', wanted_NACA);
end

x_array = typecast(db.x_2Dmat{ID},'double');
y_array = typecast(db.y_2Dmat{ID},'double');
ni = db.ni(ID);
nj = db.nj(ID);

for i = 1:nj
    x(i,:) = x_array(1+(i-1)*ni:ni+(i-1)*ni);
end
x;

for i = 1:nj
    y(i,:) = y_array(1+(i-1)*ni:ni+(i-1)*ni);
end
y;

figure
hold all
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

% axis equal
grid on
grid minor
title("The Mesh");
ylabel("Y [-]")
xlabel("X [-]")