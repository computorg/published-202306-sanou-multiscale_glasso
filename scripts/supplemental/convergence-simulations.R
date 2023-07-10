library(parallel)
library(mglasso)
source("./scripts/supplemental/convergence-functions.R")

library(reticulate)
extra_pack = c("scipy == 1.7.1", "scikit-learn", "numpy == 1.22.4",
               "six", "matplotlib")
reticulate::py_install(packages = c("pylearn-parsimony", extra_pack),
                       envname = "test",
                       method = "conda",
                       conda = "auto",
                       python_version = 3.8,
                       pip = TRUE,
                       pip_options = 'git+https://github.com/desanou/pylearn-parsimony.git@mglasso_integration')

reticulate::use_condaenv("test", required = TRUE)
reticulate::py_config()

reticulate::source_python('./scripts/supplemental/optimization-algorithms.py')

# Settings ----------------------------------------------------------------

path_conv <- "./data/"

alpha     <- rep(1/5, 5)
ngroup    <- length(alpha)
pi        <- diag(0.75, ngroup)

pp <- c(10, 20, 40)
nn <- c(20)
model <- c("erdos", "scale_free", "block_diagonal")
eps <- c(1e0, 1e-1, 1e-2)
nsimu <- c(1)
rho_sbm <- c(0.95)

simu_settings = expand.grid(pp = pp, nn = nn, model = model, eps = eps, nsimu = nsimu, rho_sbm = rho_sbm)
simu_settings

list_ii_rho <- split(simu_settings, seq(nrow(simu_settings)))

no_cores    <- min(75, length(list_ii_rho))

runtime_conv_test_set1 <- system.time(
  conv_test_set1 <- mclapply(
    list_ii_rho,
    FUN = one_simu_convergence,
    max_iter_if_no_conv = 1e4,
    mc.cores = no_cores,
    pi = pi,
    alpha = alpha,
    path_conv = path_conv)
)

#rds
saveRDS(conv_test_set1,
     file = paste0(path_conv, "conv_test_set1.rds"))

