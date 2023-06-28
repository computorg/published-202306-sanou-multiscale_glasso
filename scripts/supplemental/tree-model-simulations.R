library(ape)
library(phytools)
library(huge)
library(kernlab)
library(parallel)
library(cvxclustr)
library(mclust)
library(ggplot2)
library(huge)
library(mglasso)

source("./scripts/supplemental/tree-model-functions.R")

# Settings ----------------------------------------------------------------
## Model -------------------------------------------------------------------
path_cvx <- "./data/"

pp <- c(40)
nn <- c(20, 40, 80)
model <- c("tree")
nsimu <- 1:35

simu_settings = expand.grid(nsimu = nsimu, pp = pp, nn = nn, model = model)
simu_settings

list_ii_rho <- split(simu_settings, seq(nrow(simu_settings)))

## Simulation --------------------------------------------------------------
mc_cores    <- min(45, length(list_ii_rho))
RNGkind("L'Ecuyer-CMRG")


system.time(rand_config_tree_set5_30simus_tvmax1 <- mclapply(
  list_ii_rho, 
  FUN = one_simu_extended2, 
  tvmax = 1,
  mc.cores = mc_cores))

save(rand_config_tree_set5_30simus_tvmax1, 
     file = paste0(path_cvx, "rand_config_tree_set5_30simus_tvmax1.rds"))

rand_config_tree_set5_xxsimus_tvmax1 <- rand_config_tree_set5_30simus_tvmax1
obj_len <- sapply(rand_config_tree_set5_xxsimus_tvmax1, length)
rand_config_tree_set5_xxsimus_tvmax1 <- rand_config_tree_set5_xxsimus_tvmax1[which(obj_len>1)]

save(rand_config_tree_set5_xxsimus_tvmax1,
     file = paste0(path_cvx, "rand_config_tree_set5_xxsimus_tvmax1.rds"))


