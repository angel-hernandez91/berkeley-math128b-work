n = 50;
A = delsq(numgrid('S', n + 2)) * (n + 1)^2;
b = ones(n^2, 1);

x_exact = A\b;
Ec = zeros(1, 1000);

x=0*b;
r=b;
p=b;


for i = 1:1000
  Ap=A*p;
  rdotr=r'*r;

  alpha=rdotr/(p'*Ap);
  x=x+alpha*p;
  r=r-alpha*Ap;
  beta=(r'*r)/rdotr;
  p=r+beta*p;
  
  Ec(i) = norm((x - x_exact), Inf);
end