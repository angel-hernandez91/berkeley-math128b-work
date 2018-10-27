clear;
clc;
N = 50;
th = zeros(1, N);
mu = 0.2;
X = 2;
dy = 0.05;
v0 = 0.2;
TOL = 1e-8;
max_iterations = 100;

%T is the table of errors ||theta_k - theta_{k + 1}||_inf
[th_approx, T] = newton_system(th, mu, X, dy, v0, TOL, N, max_iterations);

%plotting
x = [0, cumsum(dy*tan(th_approx))];
y = -dy*(0:N)';

plot(x, y, '.-') 
axis equal
drawnow