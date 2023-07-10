
library(ape)
library(phytools)
library(Matrix)

tree <- rcoal(40)

plot(tree, show.tip.label = FALSE)


var_ou <- vcvPhylo(tree = tree,       #
                   anc.nodes = FALSE, # 
                   model = "OU")      # Ornstein-Uhlenbeck process


image(as(var_ou, "sparseMatrix"), sub = "", xlab = "", ylab = "")

image(as(solve(var_ou), "sparseMatrix"), sub = "", xlab = "", ylab = "")

png("cov-mat.png")
image(as(var_ou, "sparseMatrix"), sub = "", xlab = "", ylab = "")
dev.off()

png("coalescent-tree.png")
plot(tree, show.tip.label = FALSE)
dev.off()