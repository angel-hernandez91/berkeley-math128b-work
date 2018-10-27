function  th_approx = rk_continuation(th, mu, X, dy, v0, N, Nc)
    h = 1/Nc;
    b = -h*p13F(th, mu, X, dy, v0, N);
    for i = 1:N
        A = p13J(th, mu, X, dy, v0, N);
        k1 = A\b;
        
        A = p13J(th + 0.5*k1', mu, X, dy, v0, N);
        k2 = A\b;
        
        A = p13J(th + 0.5*k2', mu, X, dy, v0, N);
        k3 = A\b;
        
        A = p13J(th + k3', mu, X, dy, v0, N);
        k4 = A\b;
        
        th = th + (k1' + 2*k2' + 2*k3' + k4')/6;
    end
    th_approx = th;
end

