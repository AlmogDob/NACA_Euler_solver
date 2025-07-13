clc; clear; close all

NACA = '4710';
num_of_points = 100;

% ========================

m = str2num(NACA(1))/100;
p = str2num(NACA(2))/10;
t = str2num(NACA(3:4))/100;

delta_x = 1 / (num_of_points-1);
x        = zeros(num_of_points,1);
dy_c__dx = zeros(num_of_points,1);
y_t      = zeros(num_of_points,1);
theta    = zeros(num_of_points,1);
y_c      = zeros(num_of_points,1);
x_L      = zeros(num_of_points,1);
x_U      = zeros(num_of_points,1);
y_L      = zeros(num_of_points,1);
y_U      = zeros(num_of_points,1);
for i = 0:num_of_points-1
    x(i+1)= delta_x * i;
    y_t(i+1) = 5 * t * (0.2969 * sqrt(x(i+1)) - 0.1260 * x(i+1) - 0.3516 * x(i+1)^2 + 0.2843 * x(i+1)^3 - 0.1036 * x(i+1)^4);
    
    if x(i+1) <= p
        y_c(i+1)      = m / p^2 * (2 * p * x(i+1) - x(i+1)^2);
        dy_c__dx(i+1) = m / p^2 * (p - x(i+1));
    else
        y_c(i+1)      = m / (1 - p)^2 * ((1 - 2 * p) + 2 * p * x(i+1) - x(i+1)^2);
        dy_c__dx(i+1) = 2 * m / (1 - p)^2 * (p - x(i+1));
    end
    theta(i+1) = atan(dy_c__dx(i+1));

    x_U(i+1) = x(i+1) - y_t(i+1) * sin(theta(i+1));
    x_L(i+1) = x(i+1) + y_t(i+1) * sin(theta(i+1));
    y_U(i+1) = y_c(i+1) + y_t(i+1) * cos(theta(i+1));
    y_L(i+1) = y_c(i+1) - y_t(i+1) * cos(theta(i+1));
end

fig1 = figure(1);
hold all
plot(x_U, y_U)
plot(x, y_c, '--k')
plot(x_L, y_L)
axis equal
grid on
grid minor
box on
