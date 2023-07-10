
library(simone)
library(ghibli)
library(phyloseq)
library(rmarkdown)
library(knitr)
library(dplyr)
library(tables)
library(capushe)
library(mclust)
library(parallel)
library(huge)
library(SpiecEasi)
library(mglasso)
mglasso::install_pylearn_parsimony(envname = "rmglasso", method = "conda") ## if pylearn_parsimony not installed
reticulate::use_condaenv("rmglasso", required = TRUE)
reticulate::py_config()
path_data <- "./data/"

