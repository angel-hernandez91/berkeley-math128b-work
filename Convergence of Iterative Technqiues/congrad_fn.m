function [cg_iters, N] = congrad_fn(n)
    N = n^2;
    A = delsq(numgrid('S', n + 2)) * (n + 1)^2;
    b = ones(n^2, 1);
    cg_iters = 0;
    
    

    x=0*b;
    x_exact = A\b;
    init_error = norm((x - x_exact), Inf);
    Ec = norm((x - x_exact), Inf);
    r=b;
    p=b;


    while Ec > (1e-6)*init_error
      Ap=A*p;
      rdotr=r'*r;

      alpha=rdotr/(p'*Ap);
      x=x+alpha*p;
      r=r-alpha*Ap;
      beta=(r'*r)/rdotr;
      p=r+beta*p;

      Ec = norm((x - x_exact), Inf);
      cg_iters = cg_iters + 1;
    end
end