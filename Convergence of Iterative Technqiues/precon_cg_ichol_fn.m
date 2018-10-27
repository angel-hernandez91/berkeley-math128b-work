function [pci_iters, N] = precon_cg_ichol_fn(n)
    N = n^2;
    A = delsq(numgrid('S', n + 2)) * (n + 1)^2;
    b = ones(n^2, 1);
    pci_iters = 0;
    
    x=0*b;
    x_exact = A\b;
    Epci = norm((x - x_exact), Inf);
    init_error = norm((x - x_exact), Inf);
    R = ichol(A);
    
    
    r=b;
    p=R\(R'\b);
    z = p;


    while Epci > (1e-6)*init_error  
      Ap=A*p;
      rdotr = r'*z;

      alpha = rdotr/(p'*Ap);

      x = x + alpha*p;
      r = r - alpha*Ap;
      z = R\(R'\r);
      beta = (r'*z)/rdotr;
      p = z + beta*p;

      Epci = norm((x - x_exact), Inf);
      pci_iters = pci_iters + 1;
    end

end