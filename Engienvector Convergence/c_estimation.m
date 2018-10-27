%Answer to part 2 regarding convergence:
%Ran after running inverse_iteration.m

    format long
    for i=1:4
        p = polyfit(xplot{i},error_matrix{i}, 1); %linear approximation for the corresponding plot of 
                                                  %each choosen value of mu
        p_matrix{i} = p;
        slope_vector(i) = exp(p_matrix{i}(1)); %approximate slope, i.e. C
    end

    disp(slope_vector)

    %the values in slope_vector are approximately equal to 1. Hence, we can say
    %C ~ 1 (that is, C is approximately 1). So the convergence is linear. 

    %Since K is a symmetric matrix exact convergence is
    %O(|lambda_2/lambda_1|^2m) where m is the max number of iterations. Since
    %the eigenvalues are exact m = 1. So the expected convergence rate is
    %exp(|lambda_2/lambda_1|^2) = exp((0.4675)^2) = exp(0.2186)
    