## microbial associations in gut data -------------------------------------
#Load  
load(paste0(path_data, "amgut1.filt.rda"))
dta <- amgut1.filt

#Transform data using `clr` approach or non paranormal method
dta <- t(clr(dta + 1 , 1))

#Stars selection  
mb_out      <- neighbor_select(data = dta, lambda_min = 1e-2, nlambda = 20, 
                               nresamples = 50, verbose = TRUE, estim_var = 0.05)
lambda1     <- mb_out$lambda_opt

## Sequence 1 of lambda2 values 
pen_params <- seq_l2_l1_fixed(dt = dta, l1 = lambda1, nl2 = 20, l2_max = TRUE) 
system.time(mgl_amgut_rev <- lapply(pen_params, 
                                    FUN = mglasso_pair_param, 
                                    X_ = dta, 
                                    type = "initial"))
save(mgl_amgut_rev, 
     file = paste0(path_real_data, "mgl_amgut_rev_l2_seq0to1_20val.RData"))

## Sequence 2 of lambda2 values 
pen_params <- seq_l2_l1_fixed(dt = dta, l1 = lambda1, nl2 = 20, l2_max = 20) 
pen_params[[1]] <- NULL # remove l2 = 0, already computed in prev list
system.time(mgl_amgut_rev_set2 <- mclapply(pen_params, 
                                           FUN = mglasso_pair_param, 
                                           X_ = dta, 
                                           type = "initial",
                                           mc.cores = 19))
save(mgl_amgut_rev_set2, 
     file = paste0(path_real_data, "mgl_amgut_rev_l2_seq1to20_20val.RData"))

## Sequence 3 of lambda2 values 
pen_params <- seq_l2_l1_fixed(dt = dta, l1 = lambda1, nl2 = 20, l2_max = 4.3) 
pen_params[[1]] <- NULL # remove l2 = 0, already computed in prev list
system.time(mgl_amgut_rev_set3 <- mclapply(pen_params, 
                                           FUN = mglasso_pair_param, 
                                           X_ = dta, 
                                           type = "initial",
                                           mc.cores = 19))
save(mgl_amgut_rev_set3, 
     file = paste0(path_real_data, "mgl_amgut_rev_l2_seq0to4_20val.RData"))
