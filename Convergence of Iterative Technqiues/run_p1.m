Jacobi;
Gauss_Siedel;
SOR;
ConGrad;
preconditioned_cg_tridiag;
preconditioned_cg_ichol;

semilogy(1:1000,Ej, 1:1000,Eg, 1:1000,Ew, 1:1000,Ec, 1:1000,Epc, 1:1000,Epci);
legend('Jacobi', 'Gauss-Seidel', 'SOR', 'CG', 'Tridiag PreConCG', 'Ichol PreconCG');
title('Convergence of Iterative Techniques')
xlabel('Iteration');
ylabel('Error Norm');