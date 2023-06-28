## methylation and transcriptomic genotypes in poplar ----------------------

# Data
epit_sparse <- readRDS(paste0(path_data, "epit-spca-select.rds"))

#Selection of sparsity penalty parameter
mb_out      <-
  neighbor_select(
    data = epit_sparse,
    lambda_min = 1e-3,
    nlambda = 50,
    nresamples = 50,
    verbose = TRUE,
    estim_var = 0.05
  )

lambda1_genot     <- mb_out$lambda_opt

pen_params_genot <-
  seq_l2_l1_fixed(
    dt = epit_sparse,
    l1 = lambda1_genot,
    nl2 = 20,
    l2_max = 30.94
  )

mgl_epit_sparse_geno <- lapply(pen_params_genot,
                               FUN = mglasso_pair_param,
                               X_ = epit_sparse,
                               type = "initial")

saveRDS(paste0(path_data, "mgl_epit_sparse_geno.rds"))