#' Discover significant latent variation in data.
#' 
#' Calculate statistical significance of latent variation using jackstraw.
#' @importFrom jackstraw permutationPA
#' @importFrom graphics par
#' @import corpcor
#' @param mat Gene expression matrix, columns are cells and rows are genes.
#' @param B A number (a positive integer) of resampling iterations.
#' @param threshold A numeric value between 0 and 1 to threshold p-values.
#' @param log Comppute PCA of logarithm.
#' @param scale Scale the gene expression.
#' @return Number of significant pcs by permutation parallen analysis.
#' @export
sig.pcs <- function(mat, B = 20, threshold = 0.05, log = TRUE, scale = TRUE){
  if (log){
    mat <- log2(mat + 1)
  }
  if (scale){
    mat <- t(scale(t(mat)))
  }
  PA = permutationPA(mat, B = B, threshold = threshold)
  par(mfrow = c(1, 3))
  plot(PA$p, pch = 20, main = "Permutation Parallel Analysis P-values",
       ylab = "P-values", xlab = "Principal Component")
  svd.out = fast.svd(mat)
  plot(svd.out$d^2/sum(svd.out$d^2), pch = 20, main = "The scree plot",
       xlab = "PC", ylab = "Percent Variance Explained")
  plot(svd.out$d[1] * svd.out$v[, 1], pch = 20, main = "1st PC",
       xlab = "Observation", ylab = "Magnitude")
  
  cat("\npcs estimate: ", PA$r)
  
}

#' Discover high discriminant genes in unsupervised manner
#' 
#' Calculate statistical significance of variables driving latent variation using jackstraw
#' @importFrom jackstraw jackstraw_PCA
#' @importFrom stats p.adjust
#' @import corpcor
#' @param mat Gene expression matrix, columns are cells and rows are genes.
#' @param r A number (a positive integer) of significant principal components.
#' @param s A number (a positive integer) of “synthetic” null variables. Out of m variables, s variables are independently permuted.
#' @param B A number (a positive integer) of resampling iterations. There will be a total of s*B null statistics.
#' @param adj.method Correction method for p-values, see \code{\link[stats]{p.adjust}}.
#' @param log Comppute PCA of logarithm.
#' @param scale Scale the gene expression.
#' @param verbose	A logical specifying to print the computational progress.
#' @param ... Additional arguments passed on to \code{\link[jackstraw]{jackstraw_PCA}}.
#' @return Statistical results of jackstraw.
#' @export
sig.var <- function(mat,
                    r = NULL,
                    s = 100,
                    B = 400, 
                    adj.method = "fdr", 
                    log = TRUE, 
                    scale = TRUE, 
                    verbose = FALSE,
                    ...){
  mat <- as.matrix(mat)
  if (log){
    mat <- log2(mat + 1)
  }
  if (scale){
    mat <- t(scale(t(mat)))
  }
  # cat(head(mat))
  js.pca <- jackstraw_PCA(mat, r = r, s = s, B = B, verbose = verbose, ...)
  padj <- p.adjust(js.pca$p.value, method = adj.method)
  names(padj) <- rownames(mat)
  js.pca$padj <- padj
  js.pca$method <- adj.method
  
  return(js.pca)
  
}

#' SNN-Cliq clustering.
#' 
#' Clustering cells using SNN-Cliq algorithm.
#' @importFrom networkD3 forceNetwork JS saveNetwork
#' @import magrittr
#' @param dmat Paire-wised distance of cells.
#' @param knn Number of nearest neighbors to calculate shared nearest neighbors as secondary distance.
#' @param merge_cutoff merge_cutoff threshold.
#' @param r_cutoff r_cutoff threshold.
#' @return Clucstering labels
#' @export
snn_cliq <- function(dmat,
                     knn = 3,
                     merge_cutoff = 0.5,
                     r_cutoff = 0.7){
  # require(rPython)
  dmat <- as.matrix(dmat)
  rpy <- system.file("rpy", "snncliq.py", package = "sctools")
  rPython::python.load(rpy)
  rPython::python.assign("dmat_vector", as.vector(dmat))
  rPython::python.assign("shape",dim(dmat))
  rPython::python.assign("k", knn)
  rPython::python.assign("merge_cutoff", merge_cutoff)
  rPython::python.assign("r_cutoff", r_cutoff)
  rPython::python.exec("dmat = np.array(dmat_vector);dmat = np.reshape(dmat,shape)")
  
  # SNN-Cliq clustering
  rPython::python.exec("snn = SNN(dmat=dmat,k=k)")
  rPython::python.exec("cell_labels = snn_cliq(snn, merge_cutoff=merge_cutoff,r_cutoff=r_cutoff)")
  cell_labels <- rPython::python.get("cell_labels")
  snn <- matrix(unlist(rPython::python.get("snn")), ncol = 3, byrow = TRUE, dimnames = list(NULL,c("source","target","value")))

  # Visulize network
  links <- as.data.frame(snn)
  nodes <- data.frame(name = colnames(dmat), group = cell_labels)
  forceNetwork(Links = links, Nodes = nodes,
               Source = "source", Target = "target",
               colourScale = JS("d3.scaleOrdinal(d3.schemeCategory10);"),
               Value = "value", NodeID = "name",
               Group = "group", opacity = 0.8,
               fontSize = 20,zoom = TRUE,legend=TRUE) %>%
  saveNetwork(file = 'SNN-Network.html')
  cat("save SNN-Network to SNN-Network.html")
  
  return(list(cell_labels=cell_labels,snn=links))
  
}

#' Spectral clustering.
#' 
#' Clustering cells based on adaptive gauss kernel.
#' @importFrom SNFtool spectralClustering estimateNumberOfClustersGivenGraph
#' @param Dist Paire-wised distance matrix.
#' @param numk Number of clusters number. If NULL, estimate number will be used.
#' @param ka To equalize the number of neighbors we set the value sigma(i) for each cell i to the distance to its kath nearest neighbor. Default is 10.
#' @param kNN Additional affinities beyond Kth nearest neighbor are set to zero. If NULL, kNN will be set the 3ka.
#' @param type The variants of spectral clustering to use, see \code{\link[SNFtool]{spectralClustering}}.
#' @param estimateRage Estimate range.
#' @param scaleDist scale distance.
#' @return Spectral clustering labels, and affinity matrix.
#' @references David van Dijk et al. MAGIC: A diffusion-based imputation method reveals gene-gene interactions in single-cell RNA-sequencing data. BioRxiv 2017.
#' @references Bo Wang et al. Similarity network fusion for aggregating data types on a genomic scale. Nature Methods 2014.
#' @export
aspecc <- function(Dist, numk = NULL, ka=10, kNN=NULL, type = 3, estimateRage = 3:8, scaleDist = TRUE){
  
  d <- as.matrix(Dist)
  if (scaleDist){
    d <- scale(d, center = FALSE)
  }
  N = nrow(d)
  if (is.null(kNN)){
    kNN <- ka * 3
  }
  
  # Adapt kernel, the kernel is wider in sparse areas and smaller in dense areas.
  cat("Calculate adative kernel")
  sigma = matrix(apply(d, 1, sort)[ka+1,], nrow = N, ncol = N)
  W = exp(-1 * as.matrix(d^2)/sigma)
  if (kNN < N-1){
    Kthr = apply(W, 1, sort, decreasing = TRUE)[kNN+1,]
    for (i in 1:N){
      nn = W[i,]
      nn[nn < Kthr[i]] = 0
      W[i,] = nn
    }
  }
  W = (W + t(W)) / 2
  estimationResult = estimateNumberOfClustersGivenGraph(W, estimateRage)
  print(as.data.frame(estimationResult,
                      col.names = c("K1","K12","K2","K22"),
                      row.names = "EstimateK"))
  if (is.null(numk)){
    numk <- estimationResult[[1]]
  }
  labels = spectralClustering(W, numk)

  return(list(labels = labels, affinity = W))
  
}

#' Extract level of reprogramming cells.
#' 
#' Extract level of reprogramming cells.
#' @param mat Expression matrix.
#' @return factor represents timepoints
#' @export
annofactor <- function(mat){
  cells <- colnames(mat)
  days <- gsub("(mef|osk_d0|osk_d1|osk_d2|osk_d3|osk_d4|osk_d5|osk_d6|osk_d7|osk_d8|ips|ovsvk_d0|ovsvk_d1|ovsvk_d2|ovsvk_d3|ovsvk_d4|esc).*",
               "\\1",
               cells)
  lv <- c("mef","osk_d0","osk_d1","osk_d2","osk_d3","osk_d4","osk_d5","osk_d6","osk_d7","osk_d8","ips","ovsvk_d0","ovsvk_d1","ovsvk_d2","ovsvk_d3","ovsvk_d4","esc")
  
  return(factor(days,levels = lv[lv %in% unique(days)]))
  
}

#' Detect outliers.
#' 
#' Use Isolation Forest algorithm to detect outliers of cells.
#' @description The implement of iforest is dependent on numpy, scipy and sklearn in python, make sure that they have been installed in your enviroment.
#' @param mat Gene expression matrix, columns are cells and rows are genes.
#' @param n_estimators The number of base estimators in the ensemble.
#' @param outliers_fraction float in (0, 0.5). The amount of contamination of the data set, i.e. the proportion of outliers in the data set. Used when fitting to define the threshold on the decision function.
#' @param random_state RandomState instance.
#' @param n_jobs The number of jobs to run in parallel for both `fit` and `predict`. If -1, then the number of jobs is set to the number of cores.
#' @param verbose Controls the verbosity of the tree building process.
#' @return list of outlier score and indice.
#' @export
#' @references Liu, Fei Tony, Ting, Kai Ming and Zhou, Zhi-Hua. "Isolation forest." Data Mining, 2008. ICDM'08. Eighth IEEE International Conference on.
#' @references Liu, Fei Tony, Ting, Kai Ming and Zhou, Zhi-Hua. "Isolation-based anomaly detection." ACM Transactions on Knowledge Discovery from Data (TKDD) 6.1 (2012): 3.
iforest <- function(mat,
                    n_estimators = 1000,
                    outliers_fraction = 0.1,
                    random_state = 1,
                    n_jobs = 1){
  # require(rPython)
  mat <- as.matrix(mat)
  cells <- colnames(mat)
  rPython::python.exec("import numpy as np;from sklearn.ensemble import IsolationForest")
  rPython::python.assign("mat_vector", as.vector(mat))
  rPython::python.assign("shape",dim(mat))
  rPython::python.assign("n_estimators", n_estimators)
  rPython::python.assign("outliers_fraction", outliers_fraction)
  rPython::python.assign("random_state", random_state)
  rPython::python.assign("n_jobs", n_jobs)
  rPython::python.exec("mat = np.array(mat_vector);mat = np.reshape(mat,shape)")
  rPython::python.exec("iForest = IsolationForest(max_samples='auto',\
                       n_estimators = n_estimators,\
                       contamination = outliers_fraction,\
                       random_state = np.random.RandomState(random_state),\
                       n_jobs = n_jobs)")
  cat("construct iTrees and fit the data...")
  rPython::python.exec("iForest.fit(mat.T)")
  rPython::python.exec("scores_pred = list(iForest.decision_function(mat.T))")
  rPython::python.exec("y_pred = list(iForest.predict(mat.T))")
  scores_pred <- rPython::python.get("scores_pred")
  y_pred <- rPython::python.get("y_pred")
  pred <- data.frame(idx = 1:length(cells),
                     score = scores_pred, 
                     is_outlier = y_pred,
                     row.names = cells,
                     stringsAsFactors = FALSE)
  call <- list(n_estimators = n_estimators,
               outliers_fraction = outliers_fraction,
               random_state = random_state)
  
  return(list(pred = pred, call = call))
  
}
