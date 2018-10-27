n = 50;
A = delsq(numgrid('S', n + 2)) * (n + 1)^2;
b = ones(n^2, 1);

x_exact = A\b;
Epc = zeros(1, 1000);
M = gallery('tridiag', diag(A,-1), diag(A), diag(A,1));
x=0*b;
r=b;
p=M\b;
z = p;


for i = 1:1000
  Ap=A*p;
  rdotr = r'*z;

  alpha = rdotr/(p'*Ap);
  
  x = x + alpha*p;
  r = r - alpha*Ap;
  z = M\r;
  beta = (r'*z)/rdotr;
  p = z + beta*p;
  
  Epc(i) = norm((x - x_exact), Inf);
end