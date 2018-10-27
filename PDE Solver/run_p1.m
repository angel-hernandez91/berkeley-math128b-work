n = 40;
h = 1/n;
x = 0: h: 1;
y = 0: h: 1;
z = zeros(n+1, n+1);

TOL = 1e-8;

%boundary conditions
for i = 1: n
    z(i, 1) = 3*y(i)*(1-y(i));
    z(i, n+ 1) = 3*y(i)*(1-y(i));
end

[z, error] = pde_solver(z, h, n, TOL);

%plotting
[X, Y] = ndgrid(x, y);
surf(X, Y, z')
axis equal
shading interp
cameramenu
lighting phong
camlight right