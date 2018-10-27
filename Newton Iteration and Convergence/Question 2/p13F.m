function F = p13F(th, mu, X, dy, v0, N)
    %initializing vectors
    F = zeros(1, N);
    tan_vector = zeros(1, N);
    v = zeros(1, N);
    w = zeros(1, N);
    cosine_sum = zeros(1, N);
    v_cos_sum = zeros(1, N);
    g = 32.17;
    
    %compute v
    for i = 1:N
        for j = 1:i
            cosine_sum(i) = 1/cos(th(i));
        end
        v(i) = sqrt(v0^2 + 2*g*i*dy - 2*mu*dy*sum(cosine_sum));
        
    end
    
    %compute 1/v^2*cos(theta)
    for i = 1:N
        v_cos_sum(i) = 1/((v(i)^3)*(cos(th(i))));
    end
    
    %compute w and tan(theta)
    for i = 1:N
        w(i) = -dy*v(i)*sum(v_cos_sum);
        tan_vector(i) = tan(th(i));
    end
    
    %compute F
    for i = 1: N-1
        F(i) = (sin(th(i+1))/v(i+1))*(1 - mu*w(i + 1)) - (sin(th(i))/v(i))*(1 - mu*w(i));
    end
    
    F(N) = dy*sum(tan_vector) - X; %last element in F
    F = F';
end