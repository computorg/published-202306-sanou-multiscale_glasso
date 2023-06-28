get_optim_algo_infos <- function(results_method, prec, l1, l2) {

  res <- results_method
  
  dtf <- data.frame(
    method = res[["method"]],
    time_cumu = res[["time_cumu"]],
    niter = res[["niter"]],
    costf = unlist(res[["fvalue"]]),
    prec = prec,
    l1 = l1,
    l2 = l2,
    converged = res[["converged"]]
  )
  
  return(dtf)
}

#pen_params: list in which each element is a vector with l1 and l2 seq_l2_l1
# add pair param in one config
one_simu_convergence <- function(list_ii_rho, verbose = FALSE, 
                                 max_iter_if_no_conv = 1e4, pen_params=c(0.5, 10),
                                 pi, alpha, path_conv) {
  
  ii    = list_ii_rho$nsimu
  n     = list_ii_rho$nn
  rho   = list_ii_rho$rho_sbm
  p = list_ii_rho$pp
  model = as.character(list_ii_rho$model)
  prec = list_ii_rho$eps
  
  KK        <- length(alpha)
  config    <- one_config(n, p, pi, alpha, rho)
  
  data      <- sim_data(p, n/p, model, prob_mat = pi, alpha = alpha, rho = rho)
  
  rho.admm = ifelse(pen_params[2] == 0, 1, pen_params[2]*2*1e-3)
  
  res_conesta <- conesta_algo(scale(data$X), lam1 = pen_params[1], lam2 = pen_params[2], 
                              max_iter = max_iter_if_no_conv, prec = prec)
  res_sgd <- sgd_algo(X=scale(data$X), lam1=pen_params[1], lam2=pen_params[2], max_iter=max_iter_if_no_conv, prec=prec)
  res_admm <- admm_algo(X=scale(data$X), lam1 = pen_params[1], lam2 = pen_params[2], 
                        max_iter = max_iter_if_no_conv, prec = prec, rho = rho.admm)

  res_raw <- list(res_sgd = res_sgd, res_admm = res_admm, res_conesta = res_conesta, pen_params = pen_params, 
                  config = config, data = data)
  
  res_sgd <- get_optim_algo_infos(res_sgd, prec, pen_params[1], pen_params[2])
  res_sgd$simu_label <- ii
  res_sgd$n <- n
  res_sgd$p <- p
  res_sgd$model <- model
  
  res_admm <- get_optim_algo_infos(res_admm, prec, pen_params[1], pen_params[2])
  res_admm$simu_label <- ii
  res_admm$n <- n
  res_admm$p <- p
  res_admm$model <- model
  
  res_conesta <- get_optim_algo_infos(res_conesta, prec, pen_params[1], pen_params[2])
  res_conesta$simu_label <- ii
  res_conesta$n <- n
  res_conesta$p <- p
  res_conesta$model <- model
  res_processed <- rbind(res_sgd, res_admm, res_conesta)
  
  file_name <- paste0("conv_simu_", ii, "_n", n, "_p", p, "_eps", prec, model, ".rds")
  path <- paste0(path_conv, file_name)
  
  assign(file_name, res_raw)
  do.call(saveRDS, c(lapply(file_name, as.name), file=path))
  
  return(res_processed)
}

#' One simulation configuration
one_config <- function(n, p, pi = NULL, alpha = NULL, rho = NULL){
  return(list(n = n, p = p, pi = pi, alpha = alpha, rho = rho))
}
