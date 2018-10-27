n = 50;
A = delsq(numgrid('S', n + 2)) * (n + 1)^2;
b = ones(n^2, 1);

x_exact = A\b;
L = -tril(A, -1);
U = -triu(A, 1);
D = diag(diag(A));
x = zeros(n^2, 1);
Ej = zeros(1, 1000);
for i = 1:1000
    x = D\(L + U)*x + D\b;
    Ej(i) = (norm((x - x_exact), Inf));
    
end

