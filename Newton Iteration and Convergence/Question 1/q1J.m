function J = q1J(th, X, dy, v0, N)
    J = zeros(N, N);
    I = eye(N);
    delta = 1e-6;
    for j = 1:N
        J(:, j) = (q1F(th' + delta*I(:,j), X, dy, v0, N) - q1F(th' - delta*I(:,j), X, dy, v0, N))/(2*delta);
    end
end