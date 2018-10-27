function [gs_iters, N] = gauss_seidel_fn(n)
    N = n^2;
    A = delsq(numgrid('S', n + 2)) * (n + 1)^2;
    b = ones(n^2, 1);
    gs_iters = 0;

    x = zeros(n^2, 1);
    x_exact = A\b;
    init_error = norm((x - x_exact), Inf);
    Eg = (norm((x - x_exact), Inf));
    L = -tril(A, -1);
    U = -triu(A, 1);
    D = diag(diag(A));

    Tg = (D-L)\U;
    cg = (D-L)\b;
    while Eg > (1e-6)*init_error
        x = Tg*x + cg;
        Eg = (norm((x - x_exact), Inf));
        gs_iters = gs_iters + 1;
        if gs_iters >= 1000
            break;
        end
    end

end