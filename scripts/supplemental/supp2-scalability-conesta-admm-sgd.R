##

##
path_conv <- "./data/"
conv_test_set1 <- readRDS(paste0(path_conv, "conv_test_set1.rds"))
conv_res <- do.call(rbind, conv_test_set1)

##
by_vars <- conv_res %>% group_by(simu_label, method, p, prec, model)
grouped_conv <- by_vars %>% summarise(
  time_cumu = max(time_cumu)/60,
  converged = unique(converged)
)

##
prec.label <- c("epsilon == 0.01", "epsilon == 0.1", "epsilon == 1")
names(prec.label) <- c("0.01", "0.1", "1")
model.label <- c("Stochastic block diagonal", "Erdös-Rényi", "Scale-free")
names(model.label) <- c("block_diagonal", "erdos", "scale_free")

plt <- ggplot(grouped_conv, 
              aes(x=p, y=time_cumu, color = method)) +
  geom_line()+
  facet_grid(model ~ prec, labeller = labeller(model = as_labeller(model.label),
                                               prec = as_labeller(prec.label, label_parsed)))+
  ylab("Running time (minutes)") +
  xlab("Number of variables") +
  scale_color_manual(
    name = "Method",
    labels = c("ADMM", "CONESTA", "Subgradient"),
    values = c("brown1", "forestgreen", "deepskyblue3")
  )
plt