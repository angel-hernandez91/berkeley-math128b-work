clear;
clc;
N = 20;
th = zeros(1, N);
X = 2;
dy = 0.2;
v0 = 0;
TOL = 1e-2;
max_iterations = 100;
th_approx = q1_newton_system(th, X, dy, v0, TOL, N, max_iterations);


%plotting
x = [0, cumsum(dy*tan(th_approx))];
y = -dy*(0:N)';

plot(x, y, '.-') 
axis equal
drawnow