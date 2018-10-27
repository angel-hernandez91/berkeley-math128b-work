function [j_iters, N] = Jacobi_fn(n)
    N = n^2;
    A = delsq(numgrid('S', n + 2)) * (n + 1)^2;
    b = ones(n^2, 1);
    j_iters = 0;
    
    x = zeros(n^2, 1);
    x_exact = A\b;
    init_error = norm((x - x_exact), Inf);
    Ej = (norm((x - x_exact), Inf));
    L = -tril(A, -1);
    U = -triu(A, 1);
    D = diag(diag(A));
   
    
    while Ej > (1e-6)*init_error
        x = D\(L + U)*x + D\b;
        Ej = (norm((x - x_exact), Inf));
        j_iters = j_iters + 1;
        if j_iters >= 1000
            break;
        end
    end

end