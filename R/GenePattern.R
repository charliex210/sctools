#' Heatmap for gene expression pattern.
#' 
#' Plot the clustering result of cells and genes 
#' @importFrom grDevices colorRampPalette
#' @importFrom pheatmap pheatmap
#' @param mat Expression matrix, columns are cells and rows are genes.
#' @param filename Filename to save your heatmap.
#' @param vmax The maximum value to plot.
#' @param vmin The minimum value to plot.
#' @param log log2 tranform the data.
#' @param standard Standardize the data, result in mean = 0 and std = 1 of each gene.
#' @param hlight The genes you want to highlight in the heatmap.
#' @param cmap The pallette of the heatmap, one can generate it using colorRampPalette().
#' @param cluster_rows boolean values determining if rows should be clustered or hclust object.
#' @param cluster_cols boolean values determining if columns should be clustered or hclust object.
#' @param gene_labels A vector of gene clustering labels.
#' @param cell_labels A vector of cell clustering labels.
#' @param gene_order Specify a sequence of gene clusters order, the values must be in gene_labels.
#' @param cell_order Specify a sequence of cell clusters order, the values must be in cell_labels.
#' @param gaps_row A vector indicates gaps of row, if NULL the gaps will be setting based on the gene_labels.
#' @param gaps_col A vector indicates gaps of column, if NULL the gaps will be setting based on the cell_labels.
#' @param plot_gaps_row Make gaps on row.
#' @param plot_gaps_col Make gaps on column.
#' @param ... Additional arguments passed on to pheatmap.
#' @return Ordering cluster results.
#' @export
lpheatmap <- function(mat,
                      filename = NA,
                      vmax = 2,
                      vmin = -1, 
                      log = FALSE,
                      standard = FALSE,
                      hlight = NULL,
                      cmap = NULL,
                      cluster_rows = TRUE,
                      cluster_cols = TRUE,
                      gene_labels = NULL,
                      cell_labels = NULL,
                      gene_order =NULL,
                      cell_order =NULL, 
                      gaps_row = NULL,
                      gaps_col = NULL,
                      plot_gaps_row = T,
                      plot_gaps_col = T,
                      ...){
  
  # require(nnet)
  # require(apcluster)
  if (!is.numeric(vmax) || !is.numeric(vmin)){
    stop("The expression boundary must be numeric")
  }
  
  # Order genes
  makeGeneOrder <- 1:nrow(mat)
  if (!is.null(gene_labels)){
    cluster_rows <- FALSE
    if (!is.null(gene_order)){
      makeGeneOrder <- alter_label(gene_labels, gene_order)
      if (plot_gaps_row){
        if (is.null(gaps_row)){
          gaps_row <- head(cumsum(table(gene_labels)[gene_order]),-1)
        }
      }
    }
    else {
      makeGeneOrder <- order(gene_labels)
      if (plot_gaps_row){
        gaps_row <- head(cumsum(table(gene_labels)),-1)
      }
    }
    mat <- mat[makeGeneOrder,]
    anno_row = data.frame(Module=factor(gene_labels[makeGeneOrder]))
    rownames(anno_row) <- rownames(mat)
    output_gene_labels <- data.frame(labels=gene_labels[makeGeneOrder],row.names = rownames(mat))
  }
  else{
    gaps_row <- NULL
    anno_row <- NA
    output_gene_labels <- NULL
  }
  
  # Order cells

  makeCellOrder <- ncol(mat)
  if (!is.null(cell_labels)){
    cluster_cols <- FALSE
    if (!is.null(cell_order)){
      makeCellOrder <- alter_label(cell_labels, cell_order)
      if (plot_gaps_col){
        if (is.null(gaps_col)){
          gaps_col <- head(cumsum(table(cell_labels)[cell_order]),-1)
        }
      }
    }
    else {
      makeCellOrder <- order(cell_labels)
      if (plot_gaps_col){
        gaps_col <- head(cumsum(table(cell_labels)),-1)
      }
    }
    mat <- mat[,makeCellOrder]
    anno_col <- data.frame(Cluster = factor(cell_labels[makeCellOrder])) # bestcells = bestcells EGFP=EGFP
    rownames(anno_col) <- colnames(mat)
    output_cell_labels <- data.frame(labels=cell_labels[makeCellOrder],row.names = colnames(mat))
  }
  else{
    gaps_col <- NULL
    anno_col <- NA
    output_cell_labels <- NULL
  }
  
  if (log){
    mat <- log2(mat + 1)
  }
  if (standard){
    mat <- t(scale(t(mat)))
  }
  
  mat[mat>vmax] <- vmax
  mat[mat<vmin] <- vmin
  
  # hightlight genes
  psuGene = rep(c(""),dim(mat)[1])
  
  # hlight = tfsinRep
  if (!is.null(hlight)){
    for (g in hlight){
      idx = which(rownames(mat) == g)
      psuGene[idx] = g
    }
  }
  
  if (is.null(cmap)){
    cmap <- colorRampPalette(c("#7431A6","#000000","#F2F208"), bias = 1.7)(255)
  }
  pheatmap(mat,
           filename = filename,
           annotation_col = anno_col,
           annotation_row = anno_row,
           gaps_col = gaps_col,
           gaps_row = gaps_row,
           col = cmap,
           cluster_rows = cluster_rows,
           cluster_cols = cluster_cols,
           labels_row = psuGene,
           ...)
  
  
  
  return(list(heat=mat,cell_labels=output_cell_labels,gene_labels=output_gene_labels))
  
}

#' Re-ordering clustering labels
#' 
#' Re-ordering clustering labels for heatmpa
#' @param ori_labels Origin labels
#' @param re_labels Re-ordered labels
#' @return New labels after ordering
#' @export
alter_label <- function(ori_labels,re_labels){
  if (min(ori_labels)==0){
    ori_labels = ori_labels + 1
  }
  if (min(re_labels)==0){
    re_labels = re_labels + 1
  }
  if (length(intersect(unique(ori_labels),unique(re_labels)))!=length(unique(ori_labels))){
    stop("The given labels are not consistant with thr origin labels")
  }
  in_order=seq(length(ori_labels))
  names(in_order)=ori_labels
  tre_labelspe=sort(unique(ori_labels))
  
  l=list()
  for(i in 1:length(tre_labelspe))
  {
    l[[i]]=in_order[which(names(in_order) %in% tre_labelspe[i])]
  }
  
  new_order=c()
  for(i in 1:length(re_labels))
  {
    new_order=c(new_order,l[[re_labels[i]]])
  }
  
  return(new_order)
  
}

ap_helper <- function(mat,cell_labels=NULL,k=7,method="pearson",
                      check_outlier_pattern=TRUE,
                      check_min_size=5,check_cor_thr=-0.5,
                      maxits=1000, convits=100,lam=0.9){
  
  orig <- rownames(mat)
  negDist <- corSimMat(mat, method = "pearson")
  apres <- apclusterK(negDist,K=k,maxits=maxits, convits=maxits,lam=lam)
  
  cluster <- apres@clusters
  clres <- data.frame(idx = unlist(cluster))
  # head(clres)
  
  exemplars <- apres@exemplars
  exem_idx <- rep(FALSE, dim(clres)[1])
  exem_idx[rownames(clres)%in%names(exemplars)] <- TRUE
  labels <- rep(1:length(exemplars),unlist(lapply(cluster,length)))
  
  clres$cluster <- labels
  clres$exemplars <- exem_idx
  
  if (!check_outlier_pattern){
    
    cat("AP clustering results:")
    print(table(clres$cluster))
    
    return(clres[orig,])
    
  }
  else{
    if (is.null(cell_labels) || length(cell_labels)!=dim(mat)[2]){
      stop("Cell cluster labels are missing or wrong")
    }
    # sort the gene expression matrix
    mat <- mat[rownames(clres),]
    
    # I want to remove the gene cluster that negatively correlate with outliers 
    emexp <- mat[names(exemplars),]
    
    dum <- class.ind(cell_labels)
    clsize <- colSums(dum)
    
    # find the clusters that with size smaller than check_min_size
    check.clus <- as.numeric(which(clsize <= check_min_size))
    
    # Gene cluster that correlated with Cell cluster
    ExCorWithDum <- cor(x=t(emexp),y=dum)[,check.clus]
    outlierIdx <- rowSums(ExCorWithDum < check_cor_thr) > 0
    outlierG <- names(outlierIdx[outlierIdx==TRUE])
    
    target.cluster <- clres[outlierG,"cluster"]
    fil_mask <- clres$cluster%in%target.cluster == FALSE
    mat.fil <- mat[fil_mask,]
    clres.fil <- clres[fil_mask,]
    new_cl <- clres.fil$cluster
    
    clres.fil$cluster <- rep(1:length(unique(new_cl)),table(new_cl))
    
    cat("AP clustering and filter outlier expression pattern:")
    print(table(clres.fil$cluster))
    
    orig <- orig[orig%in%rownames(clres.fil)]
    
    return(clres.fil[orig,])
    
  }
  
}
