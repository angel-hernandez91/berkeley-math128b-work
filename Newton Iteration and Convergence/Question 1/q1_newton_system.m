function th_approx = q1_newton_system(th, X, dy, v0, TOL, N, iter)
    k = 1;
    while k <= iter
       y = q1J(th, X, dy, v0, N)\(-q1F(th, X, dy, v0, N));
       th = th + y';
       if norm(y, inf) < TOL
           break
       else           
           k = k + 1;
       end
    end
    th_approx = th;
end 