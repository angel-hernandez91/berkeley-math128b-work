clear;
clc;
N = 50;
th = zeros(1, N);
mu = 0.2;
X = 2;
dy = 0.02;
v0 = 0.2;
TOL = 1e-8;
max_iterations = 100;
Nc = 2; %converges for Nc = 2, Nc = 4. Fails to converge for Nc = 1;
rk_th_approx = rk_continuation(th, mu, X, dy, v0, N, Nc);

%T is the table of errors ||theta_k - theta_{k + 1}||_inf
[th_approx, T] = newton_system(rk_th_approx, mu, X, dy, v0, TOL, N, max_iterations);

%plotting
x = [0, cumsum(dy*tan(th_approx))];
y = -dy*(0:N)';

plot(x, y, '.-') 
axis equal
drawnow