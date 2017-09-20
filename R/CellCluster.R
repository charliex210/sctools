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
#' @importFrom jackstraw jackstraw.PCA
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
#' @param ... Additional arguments passed on to \code{\link[jackstraw]{jackstraw.PCA}}.
#' @return Statistical results of jackstraw.
#' @export
sig.var <- function(mat,
                    r = NULL,
                    s = 100,
                    B = 400, 
                    adj.method = "fdr", 
                    log = TRUE, 
                    scale = TRUE, 
                    verbose = TRUE,
                    ...){
  mat <- as.matrix(mat)
  if (log){
    mat <- log2(mat + 1)
  }
  if (scale){
    mat <- t(scale(t(mat)))
  }
  cat(head(mat))
  js.pca <- jackstraw.PCA(mat, r = r, s = s, B = B, verbose = TRUE, ...)
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
  
  return(list(cell_labels=cell_labels,snn=links))
  
}
