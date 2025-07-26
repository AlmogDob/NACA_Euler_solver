clc; clear; close all;

db_file = 'NACA.db';
db = sqlread(sqlite(db_file),'NACA_data');

% wanted_NACA = 0012;
% ID = find(db.NACA == wanted_NACA);
% ID = ID(end);
ID = 3;
if isempty(ID)
    error('No NACA=%d in Data Base', wanted_NACA);
end

x_array = typecast(db.x_2Dmat{ID},'double');
y_array = typecast(db.y_2Dmat{ID},'double');
rhos    = typecast(db.rho_2Dmat{ID},'double');
us      = typecast(db.u_2Dmat{ID},'double');
vs      = typecast(db.v_2Dmat{ID},'double');
es      = typecast(db.e_2Dmat{ID},'double');

ni = db.ni(ID);
nj = db.nj(ID);
num_points_on_airfoil = db.num_points_on_airfoil(ID);
gamma     = db.Gamma(ID);
rho_inf   = db.delta_y(ID);
Mach_inf  = db.Mach_inf(ID);
p_inf     = db.environment_pressure(ID);
a_inf     = sqrt(gamma * p_inf / rho_inf);
u_inf     = Mach_inf * a_inf;
alpha_deg = db.angle_of_attack_deg(ID);
NACA      = db.NACA(ID);

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

figure
hold all
contourf(x, y, M, 300, 'LineStyle','none');
colorbar;
colormap('turbo');
axis equal
title(sprintf('NACA: %04d,  Mach: %.2f,  $\\alpha=%.2f^\\circ$', NACA, Mach_inf, alpha_deg),'Interpreter','latex')
% 
% factor = 0.0001;
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