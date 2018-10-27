function J = p13J(th, mu, X, dy, v0, N)
    J = zeros(N, N);
    I = eye(N);
    delta = 1e-6;
    for j = 1:N
        J(:, j) = (1/(2*delta))*(p13F(th' + delta*I(:,j), mu, X, dy, v0, N) - p13F(th' - delta*I(:,j), mu, X, dy, v0, N));
    end
end