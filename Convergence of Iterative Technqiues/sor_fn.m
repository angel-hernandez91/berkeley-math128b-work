function [sor_iters, N] = sor_fn(n)
    N = n^2;
    A = delsq(numgrid('S', n + 2)) * (n + 1)^2;
    b = ones(n^2, 1);
    sor_iters = 0;

    x = zeros(n^2, 1);
    x_exact = A\b;
    init_error = norm((x - x_exact), Inf);
    Ew = (norm((x - x_exact), Inf));
    L = -tril(A, -1);
    U = -triu(A, 1);
    D = diag(diag(A));
    w = 2/(1 + sqrt(8)/(n + 1));
    

    Tw = (D-w*L)\((1-w)*D + w*U);
    cw = w*((D-w*L)\b);
    while Ew > (1e-6)*init_error
        x = Tw*x + cw;
        Ew = (norm((x - x_exact), Inf));
        sor_iters = sor_iters + 1;
        if sor_iters >= 1000
            break;
        end
    end

end