library(magrittr)
library(dplyr)
library(ggplot2)

##
path_conv <- "./data/"
conv_test_set1 <- readRDS(paste0(path_conv, "conv_test_set1.rds")) ## computed from ./convergence-simulations.R
conv_res <- do.call(rbind, conv_test_set1)

## subset according to the tolerance levels
by_eps1 <- conv_res %>% 
  group_by(method, prec, model)
grouped_eps1_dat <- by_eps1 %>% 
  summarise(costf = costf, time_cumu = time_cumu/60, converged = converged)

## 
prec.label <- c("epsilon == 0.01", "epsilon == 0.1", "epsilon == 1")
names(prec.label) <- c("0.01", "0.1", "1")
model.label <- c("Stochastic block diagonal", "Erdös-Rényi", "Scale-free")
names(model.label) <- c("block_diagonal", "erdos", "scale_free")

plt <- ggplot(
  grouped_eps1_dat, aes(x=time_cumu, y=costf, color = method)) +
  geom_line()+
  facet_grid(model ~ prec, 
             labeller = labeller(prec = as_labeller(prec.label, label_parsed),
                                 model = as_labeller(model.label, label_value))) +
  ylab("Cost function value") +
  xlab("Running time (minutes)") +
  scale_color_manual(
    name = "Method",
    labels = c("ADMM", "CONESTA", "Subgradient"),
    values = c("brown1", "forestgreen", "deepskyblue3")
  )
plt
