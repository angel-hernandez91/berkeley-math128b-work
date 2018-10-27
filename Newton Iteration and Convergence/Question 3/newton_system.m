function [th_approx, T] = newton_system(th, mu, X, dy, v0, TOL, N, iter)
    k = 1;
    while k <= iter
       F = p13F(th, mu, X, dy, v0, N);
       J = p13J(th, mu, X, dy, v0, N);
       y = J\(-F);
       th_old = th;
       th = th + y';
       if norm(y, inf) < TOL
           break
       else
           error(k) = norm(th_old - th, inf);
           iters(k) = k;
           k = k + 1;
       end
    end
    th_approx = th;
    T(:, 1) = iters';
    T(:, 2) = error';
end 