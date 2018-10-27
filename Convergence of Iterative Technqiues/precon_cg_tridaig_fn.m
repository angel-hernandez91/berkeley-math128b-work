function [pct_iters, N] = precon_cg_tridaig_fn(n)
    N = n^2;
    A = delsq(numgrid('S', n + 2)) * (n + 1)^2;
    b = ones(n^2, 1);
    pct_iters = 0;
    
    

    x=0*b;
    x_exact = A\b;
    init_error = norm((x - x_exact), Inf);
    Epc = norm((x - x_exact), Inf);
    M = gallery('tridiag', diag(A,-1), diag(A), diag(A,1));
    x=0*b;
    r=b;
    p=M\b;
    z = p;


    while Epc > (1e-6)*init_error
      Ap=A*p;
      rdotr = r'*z;

      alpha = rdotr/(p'*Ap);

      x = x + alpha*p;
      r = r - alpha*Ap;
      z = M\r;
      beta = (r'*z)/rdotr;
      p = z + beta*p;

      Epc = norm((x - x_exact), Inf);
      pct_iters = pct_iters + 1;
    end
end

