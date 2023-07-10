#' simulate data with given graph structure
#'
#' @param p 
#' @param np_ratio 
#' @param structure 
#' @param alpha 
#' @param prob_mat 
#' @param rho 
#' @param g 
#' @param inter_cluster_edge_prob 
#' @param p_erdos 
#' @param verbose 
#'
#' @return A list:
#' graph : precision
sim_data <- function(p = 20,
                     np_ratio = 2,
                     structure = c("block_diagonal", "hub", "scale_free", "erdos", "tree"),
                     alpha,
                     prob_mat,
                     rho,
                     g,
                     inter_cluster_edge_prob = 0.01,
                     p_erdos = 0.1,
                     verbose = FALSE){
  
  structure <- match.arg(structure)
  n = round(np_ratio * p)
  htree = NULL
  
  switch (structure,
          erdos = {
            network <- rNetwork(p, p_erdos, 1)
            correlation <- solve(network$Theta)
            X           <- mvtnorm::rmvnorm(n, sigma = correlation)
            graph <- network$Theta
            
            if(verbose) message("Data, precision matrix and graph generated from Erdos structure.")
          },
          
          hub = {
            L = huge.generator(graph = "hub", g = 5, d = p, n = n, verbose = FALSE)
            
            X=L$data
            graph = L$omega
            correlation = L$sigma
            
            if(verbose) message("Data, precision matrix and graph generated from Hub structure.")
          },
          
          scale_free = {
            L = huge.generator(graph = "scale-free", d = p, n = n, verbose = FALSE) ## d edges graph
            
            X=L$data
            graph = L$omega
            correlation = L$sigma
            
            if(verbose) message("Data, precision matrix and graph generated from Scale free structure.")
          },
          
          tree = {
            min_eigen = -1 # Concentration matrix is not positive definite
            while(min_eigen < 0) {
              htree <- simulate_tree("coalescent", nleaves = p)
              
              var_ou <- vcvPhylo(tree = htree,       #
                                 anc.nodes = FALSE, # 
                                 model = "OU",      # Ornstein-Uhlenbeck process
                                 alpha = 1)   
              correlation = var_ou/var_ou[1,1]
              
              min_eigen <- min(eigen(correlation)$values)
            }
            
            htree <- as.hclust.phylo(htree)
            graph <- solve(correlation)
            X = mvtnorm::rmvnorm(n, sigma=correlation)
            colnames(X) <- colnames(correlation)
            
            if(verbose) message("Data, precision matrix and graph generated from Scale free structure.")
          },
          
          block_diagonal = { 
            
            if(inter_cluster_edge_prob != 0) { # inter-clusters edges added in the flipped way
              flag <- TRUE
              
              while(flag) { # To do: treat the flag case instead of generating a new model
                K <- bloc_diag(p, prob_mat, alpha, rho)
                
                target_indices <- which(K[upper.tri(K)] == 0) # Random selection of edges to be set to 1
                
                select_len <- round(length(target_indices) * inter_cluster_edge_prob)
                selected_indices <- sample(target_indices, select_len)
                
                precision_level <- unique(K[upper.tri(K)])
                precision_level <- max(precision_level[precision_level != 0]) # je prends la precision max à l'intérieur des blocs
                # et je l'affecte dans les arêtes inter-blocs
                
                K[upper.tri(K)][selected_indices] <- precision_level
                K <- as.matrix(forceSymmetric(K, uplo = "U"))
                
                flag <- any(eigen(K)$values <= 0) # Control of positive definiteness
              }
              
              correlation <- solve(K)
              graph = K
              X           <- mvtnorm::rmvnorm(n, sigma = correlation)
              
              if(verbose) message("Data, precision matrix and graph generated from block-diagonal 
                                  structure with inter-clusters edges.")
            }
            else { # Only intra-cluster edges while approximately controlling correlation level
              K <- bloc_diag(p, prob_mat, alpha, rho)
              correlation <- solve(K)
              graph = K
              X           <- mvtnorm::rmvnorm(n, sigma = correlation)
              
              if(verbose) message("Data, precision matrix and graph generated from block-diagonal 
                                  structure with only intra-clusters edges.")
            }
          }
  )
  
  return(list(X=X, graph=graph, correlation=correlation, htree = htree))
}


#' simulate a tree
simulate_tree <- function(type = "purebirth",
                          nleaves,
                          birthrate = 0.2,
                          pbiased = 0.3){
  if(type == "purebirth"){
    tree <- pbtree(b = birthrate,            # birth rate
                   n = nleaves,              # number of tips/leaves
                   type = "continuous",      # tree is simulated in continuous time
                   scale = 1,
                   tip.label = 1:nleaves)
  }else
    if(type == "coalescent"){
      tree <- rcoal(n = nleaves,
                    tip.label = 1:nleaves)
    }else
      if(type == "pda"){
        tree <- as.phylo(rtreeshape(n = nleaves,
                                    tip.number = nleaves,
                                    model = "pda")[[1]],
                         use.labels = FALSE)
      }else
        if(type == "aldous"){
          tree <- as.phylo(rtreeshape(n = nleaves,
                                      tip.number = nleaves,
                                      model = "aldous")[[1]],
                           use.labels = FALSE)
        }else
          if(type == "biased"){
            tree <- as.phylo(rtreeshape(n = nleaves,
                                        tip.number = nleaves,
                                        model = "biased",
                                        p = pbiased)[[1]],
                             use.labels = FALSE)
          }
  return(tree)
}

one_simu_extended2 <- function(list_ii_rho, verbose = FALSE, model = "tree", tvmax){
  ii = list_ii_rho$nsimu
  n = list_ii_rho$nn
  p = list_ii_rho$pp
  model = as.character(list_ii_rho$model)
  
  config <- one_config(n, p)
  
  data <- sim_data(p, n/p, model)
  #beta_true <- precision_to_regression(data$graph)
  #true_clusters <- as.numeric(colnames(data$correlation))
  tree <- data$htree
  
  mb_out      <- neighbor_select(data$X, config, lambda_min_ratio = 1e-3, 
                                 nlambda = 50, nresamples = 50, model = model, estim_var = 0.05) 
  lambda1     <- mb_out$lambda_opt
  #lambda1 = 0
  if(verbose) cat("selection of lambda1 done \n")
  
  ## CAH 
  time_cah <- system.time({
    cah_data <- hclust(dist(t(scale(data$X))), method = "ward.D")
    cah_out <- lapply(2:(p-1), function(ncl){cutree(cah_data, k = ncl)})
  }
  )
  
  if(verbose) cat("cah_out done \n")
  
  ## Spectral clustering 
  time_spec <- system.time(
    spectral_out <- lapply(2:(2*p/3), function(ncl) {specc(t(data$X), centers=ncl)@.Data})
  )
  cat("spectral clustering done \n")
  
  ### MGL
  pen_params <- seq_l2_l1_fixed(dt = data$X, l1 = lambda1, nl2 = 20, l2_max = tvmax) #set nl2 to 2 for testing purpose
  
  time_mgl <- system.time({
    #RNGkind("L'Ecuyer-CMRG")
    mgl <- lapply(pen_params,
                  FUN = mglasso_pair_param, X_ = data$X, type = "initial")
  }
  )
  if(verbose) cat("mglasso done \n")
  
  
  ## K-MEANS
  time_kmeans <- system.time(
    kmeans_out <- lapply(2:(p-1), function(ncl) {kmeans(scale(t(data$X)), centers = ncl)$cluster})
  )
  if(verbose) cat("kmeans_out done \n")
  
  ## Convex clustering
  time_cvx <- system.time({
    #RNGkind("L'Ecuyer-CMRG")
    cvx <- lapply(pen_params,
                  FUN = cvx_pair_param, X_ = data$X)
  }
  )
  if(verbose) cat("convex clustering done \n")
  
  running_times <- list(time_cah, time_kmeans, time_mgl, time_cvx, time_spec)
  
  res <- list(simu_label = ii, config = config, data = data, 
              tree = tree, cah_data = cah_data, spectral_out = spectral_out,
              pen_params = pen_params, mgl = mgl, cah_out = cah_out, kmeans_out = kmeans_out, cvx = cvx, running_times = running_times)
  
  # save
  file_name <- paste0("simu", ii, "_p", p, "_n", n, "_", model, ".rds")
  path <- paste0(path_cvx, file_name)
  assign(file_name, res)
  do.call(saveRDS, c(lapply(file_name, as.name), file = path))
  
  return(res)
}

plot_res2 <- function(dt, crit_, ncluster_, np_, method_ = NULL, main = ""){
  if (is.null(method_))
    method_ <- unique(dt$method)
  
  nclusters.labs <- c("3 clusters", "5 clusters", "6 clusters", "10 clusters", "22 clusters")
  names(nclusters.labs) <- c("3", "5", "6", "10", "22")
  
  np.labs <- c("frac(n, p) == 0.5", "frac(n, p) == 1", "frac(n, p) == 2")
  names(np.labs) <- c("0.5", "1", "2")
  
  bp_error <- ggplot(
    subset(x      = dt, 
           subset = (crit == crit_ & 
                       ncluster %in% ncluster_ & 
                       np %in% np_ & 
                       method %in% method_
           )), 
    
    aes(x     = factor(np), 
        y     = perfs, 
        fill = method
    )) + 
    geom_boxplot() + 
    facet_grid(np ~ ncluster, labeller = labeller(ncluster = nclusters.labs,  np = as_labeller(np.labs, label_parsed) )) + 
    ggtitle(main) + 
    ylab(TeX(r'($|| \beta - \hat{\beta} ||_F$)')) +
    xlab("Ratio n/p") + 
    scale_fill_manual(name = "Method", labels = c("HAC", "Convex Clustering", "K-means", "MGLasso", "Spectral Clustering"), values = rev(ghibli_palette("MarnieMedium1")))
  
  if(crit_ == "rand"){
    bp_error <- bp_error + 
      ylab("Adjusted Rand Index")
  }
  
  bp_error
}

get_perf_from_raw2 <- function(criterion, out, thresh_fuse = 1e-3, cah_kmeans_available = FALSE) {
  # mgl outs
  mgl       <- out$mgl
  cvx <- out$cvx
  # MGL clusters
  nclusters_mgl   <- sapply(mgl, function(e){length(unique(get_clusters_mgl(e$selected_Theta, fuse_thresh = thresh_fuse, p = out$config$p)))})
  
  # Cvx clusters
  nclusters_cvx <- sapply(cvx, function(e){length(unique(get_clusters_cvx(e$centroids_mat, fuse_thresh = thresh_fuse)))})
  
  #### RAND
  if(criterion == "rand") {
    method <- c("cah", "kmeans", "mglasso", "cvx", "spec")
    
    htree <- out$tree
    
    ncluster_cah <- ncluster_kmeans <- 2:(out$config$p-1)
    ncluster_spec <- 2:(length(out$spectral_out)+1)
    
    if(!cah_kmeans_available) {
      cah   <- lapply(2:(out$config$p -1), function(ncl){cutree(out$cah_data, k = ncl)})
      kmean <- lapply(2:(out$config$p -1), function(ncl) {kmeans(scale(t(out$data$X)), centers = ncl)$cluster})
      
    } else {
      cah   <- out$cah_out
      kmean <- out$kmeans_out
      spec <- out$spectral_out
    }
    
    nclusters <- c(ncluster_cah, ncluster_kmeans, nclusters_mgl, nclusters_cvx, ncluster_spec)
    vec_rep_num <- c(length(ncluster_cah), length(ncluster_kmeans), length(nclusters_mgl),
                     length(nclusters_cvx), length(ncluster_spec))
    
    methods <- rep(method, vec_rep_num)
    
    
    perfs_cah   <- sapply(cah, function(e){mclust::adjustedRandIndex(e, cutree(htree, length(unique(e))))})
    perfs_kmean <- sapply(kmean, function(e){mclust::adjustedRandIndex(e, cutree(htree, length(unique(e))))})
    
    perfs_mgl   <- sapply(mgl, 
                          function(e){
                            mglasso_cut <- get_clusters_mgl(e$selected_Theta, fuse_thresh = thresh_fuse, p = out$config$p)
                            return(mclust::adjustedRandIndex(mglasso_cut, cutree(htree, length(unique(e)))))
                          })
    
    perfs_cvx <- sapply(cvx, 
                        function(e){
                          cvx_cut <- get_clusters_cvx(e$centroids_mat, fuse_thresh = thresh_fuse)
                          return(mclust::adjustedRandIndex(cvx_cut, cutree(htree, length(unique(e)))))
                        })
    
    perfs_spec <- sapply(spec, function(e){mclust::adjustedRandIndex(e, cutree(htree, length(unique(e))))})
    
    perfs       <- c(perfs_cah, perfs_kmean, perfs_mgl, perfs_cvx, perfs_spec)
    
    res             <- data.frame(perfs)
    res$simu_label  <- out$simu_label
    res$crit        <- "rand"
    res$np          <- out$config$n/out$config$p
    
    res$method    <- methods
    res$ncluster  <- nclusters
    
    return(res)
  }
  
}

get_clusters_cvx<- function(U, fuse_thresh = 1e-3) {
  p <- ncol(U)
  clusters <- 1:p
  diffs <- as.matrix(dist(t(U))) ## Update distance matrix
  diag(diffs) <- NA
  pairs_to_merge <- which(diffs <= fuse_thresh, arr.ind = TRUE)
  
  if(nrow(pairs_to_merge) != 0)
    clusters <- mglasso:::merge_clusters(pairs_to_merge, clusters)  # merge clusters
  
  return(clusters)
}