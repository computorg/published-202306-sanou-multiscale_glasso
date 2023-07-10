## Correlation parameter set to 0.95 ---------------------------------------

### Settings ----------------------------------------------------------------
p         <- 40
seq_n     <- c(20, 40, 80)
alpha     <- rep(1 / 5, 5)
seq_rho   <- c(0.95)
seq_dnsty <- c(0.75)
type      <-
  NA     #1 ## unused to do: delete in configs_simu parameters
ngroup    <- length(alpha)
pi        <- diag(0.75, ngroup)

n_simu      <- 100
list_ii_rho <- configs_simu(n_simu, seq_rho, seq_dnsty, seq_n, type)
mc_cores    <- min(80, length(list_ii_rho))
RNGkind("L'Ecuyer-CMRG")

### Stochastic Block Diagonal -----------------------------------------------
runtime_rand100_config_p40_bdiagflip001_cor95 <- system.time(
  rand100_config_p40_bdiagflip001_cor95 <- mclapply(
    list_ii_rho, 
    FUN = one_simu_extended, 
    model = "block_diagonal",
    p=p,
    pi=pi,
    alpha=alpha,
    path_extended=path_extended,
    mc.cores = mc_cores)
)

save(rand100_config_p40_bdiagflip001_cor95, 
     file = paste0(path_data, "rand100_config_p40_bdiagflip001_cor95.RData"))

# Clustering  
# Settings:  
# - Inter-clusters edge probability $0.01$ (flip on all the missing edges)  

## Rand index  
load(paste0(path_data, "rand100_config_p40_bdiagflip001_cor95.RData")) 
dt <- rand100_config_p40_bdiagflip001_cor95

# Calculate clusters partitions with thresh_fuse as the required difference
# threshold for merging two regression vectors
list_res <-
  lapply(dt, function(e) {
    get_perf_from_raw("rand", e, thresh_fuse = 1e-6)
  })
dt_rand <- do.call(rbind, list_res)

save(dt_rand,
     file = paste0(path_data, "rand_dt_higher_cor_sbm.RData"))
