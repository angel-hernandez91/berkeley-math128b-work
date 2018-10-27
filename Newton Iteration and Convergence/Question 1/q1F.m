function F = q1F(th, X, dy, v0, N)
    F = zeros(1, N);
    v = zeros(1, N);
    tan_vector = zeros(1, N);
    g = 32.17; %gravity constant in ft/s^2
    for i = 1: N
        v(i) = sqrt(v0^2 + 2*g*i*dy);
        tan_vector(i) = tan(th(i));
    end
    for i = 1: N-1
        F(i) = sin(th(i + 1))/v(i +1) - sin(th(i))/v(i);
        
    end
    F(N) = dy*sum(tan_vector) - X;
    F = F';
end