library(mclust)
library(latex2exp)
library(ghibli)

path_cvx <- "./data/"

load("./data/rand_config_p40_n20_tree_set1.RData")
load("./data/rand_config_tree_set2_28simus_tvmax10.RData")
load("./data/rand_config_tree_set3_27simus_tvmax5.RData")
load("./data/rand_config_tree_set4_35simus_tvmax2_5.RData")
load("./data/rand_config_tree_set5_30simus_tvmax1.RData")

dt <- c(rand_config_p40_n20_tree_set1,
        rand_config_tree_set2_28simus_tvmax10,
        rand_config_tree_set3_27simus_tvmax5,
        rand_config_tree_set4_35simus_tvmax2_5,
        rand_config_tree_set5_30simus_tvmax1)

list_res1 <- lapply(dt, function(e){get_perf_from_raw2("rand", e, thresh_fuse = 1e-4, cah_kmeans_available = TRUE)})

dt_rand <- do.call(rbind, list_res1)

plot_res2(dt_rand, crit_ = "rand", c(3, 5, 6, 10, 22), np_ = c(0.5, 1, 2))
