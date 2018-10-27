n = [5, 10, 20, 40, 80];

TOL = 1e-8;

%compute "exact" solution.
    h = 1/n(end);
    x = 0: h: 1;
    y = 0: h: 1;
    z = zeros(n(end)+1, n(end)+1);



    %boundary conditions
    for i = 1: n(end)
        z(i, 1) = 3*y(i)*(1-y(i));
        z(i, n(end)+ 1) = 3*y(i)*(1-y(i));
    end

    z_exact = pde_solver(z, h, n(end), TOL);
    
%compute z for n = [5, 10, 20, 40]
for k = 1:(length(n) -1)
    h = 1/n(k);
    x = 0: h: 1;
    y = 0: h: 1;
    z = zeros(n(k)+1, n(k)+1);
    
    for i = 1: n(k)
        z(i, 1) = 3*y(i)*(1-y(i));
        z(i, n(k) + 1) = 3*y(i)*(1-y(i));
    end
    z = pde_solver(z, h, n(k), TOL);
    e{k} = z; %cell array of z for each grid size
end
%compute max norm error at overlapping gridpoints
e5 = e{1} - z_exact(1:16:end, 1:16:end);
e10 = e{2} - z_exact(1:8:end, 1:8:end);
e20 = e{3} - z_exact(1:4:end, 1:4:end);
e40 = e{4} - z_exact(1:2:end, 1:2:end);

max_norm = [norm(e5(:), inf), norm(e10(:), inf), norm(e20(:), inf), norm(e40(:), inf)]; 

%compute approximate slope of data points
h = 1./n;
h(end) = [];
p = polyfit(h, max_norm, 1);
slope = exp(p(1));
disp('approximate slope of curve: ');
disp(slope)

%plotting error vs h
loglog(max_norm, h, '.-')


