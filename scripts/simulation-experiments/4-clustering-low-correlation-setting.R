## Correlation set to 0.25 -------------------------------------------------
### Settings ----------------------------------------------------------------
p         <- 40
seq_n     <- c(20, 40, 80) 
alpha     <- rep(1/5, 5)
seq_rho   <- c(0.25, 0.95)
seq_dnsty <- c(0.75)
type      <- NA     #1 ## unused to do: delete in configs_simu parameters
ngroup    <- length(alpha)
pi        <- diag(0.75, ngroup)

n_simu      <- 100
list_ii_rho <- configs_simu(n_simu, seq_rho, seq_dnsty, seq_n, type)
mc_cores    <- min(detectCores() - 2, length(list_ii_rho)) ## update according to the 
path_extended <- path_data
## number of available cores 
RNGkind("L'Ecuyer-CMRG")

### Stochastic Block Diagonal -----------------------------------------------
runtime_rand100_config_p40_bdiagflip001_cor025 <- system.time(
  rand100_config_p40_bdiagflip001_cor025 <- mclapply(
    list_ii_rho, 
    FUN = one_simu_extended, 
    model = "block_diagonal",
    p=p,
    pi=pi,
    alpha=alpha,
    path_extended=path_extended,
    mc.cores = mc_cores)
)

save(rand100_config_p40_bdiagflip001_cor025, 
     file = paste0(path_data, "rand100_config_p40_bdiagflip001_cor025.RData"))

## Rand index  
load(paste0(path_data, "rand100_config_p40_bdiagflip001_cor025.RData")) 
dt <- rand100_config_p40_bdiagflip001_cor025

# Calculate clusters partitions with thresh_fuse as the required difference
# threshold for merging two regression vectors
list_res <-
  lapply(dt, function(e) {
    get_perf_from_raw("rand", e, thresh_fuse = 1e-6)
  })
dt_rand <- do.call(rbind, list_res)

save(dt_rand,
     file = paste0(path_data, "rand_dt_lower_cor_sbm.RData"))
