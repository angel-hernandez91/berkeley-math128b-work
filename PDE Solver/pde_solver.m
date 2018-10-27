function [z, error] = pde_solver(z, h, n, TOL)
    for k = 1:100000
       z_old = z; z_new = z;
        for i = 2:n
            for j = 2:n
                z_x = (z_old(i+1, j) - z_old(i-1, j))/(2*h);
                z_y = (z_old(i, j+1) - z_old(i, j-1))/(2*h);
                z_xy = (z_new(i+1,j+1) -z_new(i-1,j+1) -z_new(i+1,j-1)+z_new(i-1,j-1))/(4*h^2);
                z_xx = z_new(i -1, j) + z_new(i + 1, j);
                z_yy = z_new(i, j -1) + z_new(i, j + 1);

                z_new(i, j) = (2*h^2*(z_x)*(z_y)*(z_xy) - (1 + z_x^2)*(z_yy) - (1 + z_y^2)*(z_xx))/(-2*(1 + z_x^2) - 2*(1 + z_y^2));


            end
        end
        if norm(z_new - z_old, inf) < TOL
            break
        else
            error(k) = norm(z_new - z_old, inf);
            z_old = z; z = z_new;
        end
    end
end

