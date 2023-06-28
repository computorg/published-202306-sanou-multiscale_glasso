## Erdos Performances calculation ------------------------------------------------
## Launched on a cluster using 72 cores

### Settings ----------------------------------------------------------------
### Model -------------------------------------------------------------------
## NA values for some parameter mean they are not relevant
p         <- 40
seq_n     <- c(20, 40, 80)
seq_rho   <- 0.95
seq_dnsty <- NA 
type      <- NA
alpha     <- rep(1/5, 5)
ngroup    <- length(alpha)
pi        <- diag(0.75, ngroup)

config    <- one_config(n, p, pi, alpha, rho)

### Simulation --------------------------------------------------------------
n_simu      <- 50
list_ii_rho <- configs_simu(n_simu, seq_rho, seq_dnsty, seq_n, type)
no_cores    <- min(detectCores() - 2, length(list_ii_rho)) ## update according to the 
## number of available cores 

### Erdos -------------------------------------------------------------------
path_roc <- path_data # path for saving results 
runtime_roc_config_p40_erdos01 <- system.time(
  roc_config_p40_erdos01 <- mclapply(
    list_ii_rho, 
    FUN = one_simu_ROC, 
    model = "erdos",
    p=p,
    pi=pi,
    alpha=alpha,
    path_roc=path_roc,
    mc.cores = no_cores)
)

save(roc_config_p40_erdos01, 
     file = paste0(path_data, "roc_config_p40_erdos01.RData"))
save(runtime_roc_config_p40_erdos01, 
     file = paste0(path_data, "runtime_roc_config_p40_erdos01.RData"))

load(paste0(path_data, "roc_config_p40_erdos01.RData")) 
dt_full <- roc_config_p40_erdos01

#### Merge in one graph
# Three sample sizes are used and the vector c(20,40,80) is replicated 50 times
# subset the dataframe in three parts corresponding to the relevant sample sizes
index <- seq(1, 150, by = 3)
roc_dt20 <- dt_full[index]

index <- seq(2, 150, by = 3)
roc_dt40 <- dt_full[index]

index <- seq(3, 150, by = 3)
roc_dt80 <- dt_full[index]

# compute the mean over the 50 ROC curves
roc_dt20 <- get_mean_ROC_stat(roc_dt20)
roc_dt40 <- get_mean_ROC_stat(roc_dt40)
roc_dt80 <- get_mean_ROC_stat(roc_dt80)

# restructure the list result in a matrix for plot
roc_dt20 <- reformat_roc_res_for_ggplot(roc_dt20)
roc_dt20$np <- 0.5 # add a ratio n over p variable
roc_dt40 <- reformat_roc_res_for_ggplot(roc_dt40)
roc_dt40$np <- 1
roc_dt80 <- reformat_roc_res_for_ggplot(roc_dt80)
roc_dt80$np <- 2

roc_dtf_erdos <- rbind(roc_dt20, roc_dt40, roc_dt80)

save(roc_dtf_erdos, file = paste0(path_data, "roc_dtf_erdos.RData"))
