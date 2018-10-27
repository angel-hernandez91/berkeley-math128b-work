n = 50;
A = delsq(numgrid('S', n + 2)) * (n + 1)^2;
b = ones(n^2, 1);

x_exact = A\b;
L = -tril(A, -1);
U = -triu(A, 1);
D = diag(diag(A));
w = 2/(1 + sqrt(8)/(n + 1));
x = zeros(n^2, 1);
Ew = zeros(1, 1000);

Tw = (D-w*L)\((1-w)*D + w*U);
cw = w*((D-w*L)\b);
for i = 1:1000
    x = Tw*x + cw;
    Ew(i) = (norm((x - x_exact), Inf));
end
