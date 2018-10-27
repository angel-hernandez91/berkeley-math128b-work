system_size = [5, 10, 20, 50, 100];
for i = 1:length(system_size)
    [j_iters, N] = Jacobi_fn(system_size(i));
    [gs_iters, N] = gauss_seidel_fn(system_size(i));
    [sor_iters, N] = sor_fn(system_size(i));
    [cg_iters, N] = congrad_fn(system_size(i));
    [pct_iters, N] = precon_cg_tridaig_fn(system_size(i));
    [pci_iters, N] = precon_cg_ichol_fn(system_size(i));
      j_iters(i) = j_iters;
      gs_iters(i) = gs_iters;
      sor_iters(i) = sor_iters;
      cg_iters(i) = cg_iters;
      pct_iters(i) = pct_iters;
      pci_iters(i) = pci_iters;
      loglog(j_iters(i), N, '-o', gs_iters(i), N, '-o', sor_iters(i), N, '-o', cg_iters, N, '-o', pct_iters, N, '-o', pci_iters, N, '-o');
      hold on
      legend('Jacobi', 'Gauss Seidel', 'SOR', 'CG', 'PreCon CG Tridiag', 'PreCon CG Ichol');
      title('Iteration Count for a 10^-6 Error Tolerence vs System Size');
      xlabel('System Size (N)');
      ylabel('Number of Iterations');
end
