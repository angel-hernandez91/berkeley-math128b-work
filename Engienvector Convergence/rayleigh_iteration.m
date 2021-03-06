clear;
clc;
format long
%mu = 0.0281; %below e-value
%mu = 0.1013; %above e-value
%mu = 0.052;  %between e-values 0.048 and 0.058
%mu = 0.083;  %close to e-value with multiplicity > 1
mu_vector = [0.0281, 0.1013, 0.052, 0.084];
K = assemble(20);
N = size(K, 1);

for j = 1:4
    lambda_old = mu_vector(j);
    v = randn(1, N);
    v0 = v/(norm(v));
    lambda = 0;
    lambda_store = 0;
    x = 0;
    for i = 1:20
        w = (K - lambda_old*speye(N, N))\v0';
        v_new = w/norm(w);
        v0 = v_new';
        lambda = v_new'*K*v_new;
        lambda_store(i) = lambda; 
        x(i) = i;
        if abs(lambda - lambda_old) < 1e-14
            error = zeros(1, length(x));
            for k = 1: length(x)
                error(k) = abs(lambda_store(k) - lambda_store(end));
            end
            xplot{j} = x;
            error_matrix{j} = error;
            break
        end
        lambda_old = lambda_store(i); 
    end   
end
semilogy(xplot{1},error_matrix{1}, 'b');
hold on
semilogy(xplot{2},error_matrix{2}, 'r');
hold on
semilogy(xplot{3},error_matrix{3}, 'g');
hold on
semilogy(xplot{4},error_matrix{4}, 'k');
title('Part 3 Plots')
xlabel('Iterations')
ylabel('Error: lambda_k - lambda')
legend('below', 'above', 'between', 'multiplicity > 1')
