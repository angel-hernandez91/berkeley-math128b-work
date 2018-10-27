n = 50;
A = delsq(numgrid('S', n + 2)) * (n + 1)^2;
b = ones(n^2, 1);

x_exact = A\b;
Epci = zeros(1, 1000);
R = ichol(A);
M = R'*R;
x=0*b;
r=b;
p=R\(R'\b);
z = p;


for i = 1:1000
  Ap=A*p;
  rdotr = r'*z;

  alpha = rdotr/(p'*Ap);
  
  x = x + alpha*p;
  r = r - alpha*Ap;
  z = R\(R'\r);
  beta = (r'*z)/rdotr;
  p = z + beta*p;
  
  Epci(i) = norm((x - x_exact), Inf);
end