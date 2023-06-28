## Scale free Performances calculation -----------------------------------------
# Launched on a cluster using 72 cores

### Settings ----------------------------------------------------------------
### Model -------------------------------------------------------------------
# NA values for some parameter mean they are not relevant
p         <- 40
seq_n     <- c(20, 40, 80)
seq_rho   <- 0.95
seq_dnsty <- NA 
type      <- NA
alpha     <- rep(1/5, 5)
ngroup    <- length(alpha)
pi        <- diag(0.75, ngroup)

### Simulation --------------------------------------------------------------
n_simu      <- 50
list_ii_rho <- configs_simu(n_simu, seq_rho, seq_dnsty, seq_n, type)
no_cores    <- min(72, length(list_ii_rho))


### Scale-Free --------------------------------------------------------------
runtime_roc_config_p40_scalefree <- system.time(
  roc_config_p40_scalefree <- mclapply(
    list_ii_rho, 
    FUN = one_simu_ROC, 
    model = "scale_free",
    p=p,
    pi=pi,
    alpha=alpha,
    path_roc=path_roc,
    mc.cores = no_cores)
)

save(roc_config_p40_scalefree, 
     file = paste0(path_data, "roc_config_p40_scalefree.RData"))
save(runtime_roc_config_p40_scalefree, 
     file = paste0(path_data, "runtime_roc_config_p40_scalefree.RData"))

load(paste0(path_data, "roc_config_p40_scalefree.RData")) 
dt_full <- roc_config_p40_scalefree

### Merge in one graph
index <- seq(1, 150, by = 3)
roc_dt20 <- dt_full[index]

index <- seq(2, 150, by = 3)
roc_dt40 <- dt_full[index]

index <- seq(3, 150, by = 3)
roc_dt80 <- dt_full[index]

roc_dt20 <- get_mean_ROC_stat(roc_dt20)
roc_dt40 <- get_mean_ROC_stat(roc_dt40)
roc_dt80 <- get_mean_ROC_stat(roc_dt80)

roc_dt20 <- reformat_roc_res_for_ggplot(roc_dt20)
roc_dt20$np <- 0.5
roc_dt40 <- reformat_roc_res_for_ggplot(roc_dt40)
roc_dt40$np <- 1
roc_dt80 <- reformat_roc_res_for_ggplot(roc_dt80)
roc_dt80$np <- 2

roc_dtf_sfree <- rbind(roc_dt20, roc_dt40, roc_dt80)

### Save
save(roc_dtf_sfree,
     file = paste0(path_data, "roc_dtf_sfree.RData"))
