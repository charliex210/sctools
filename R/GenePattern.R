#' Variable selction based on elasticnet penalty.
#' 
#' Fits a multinomial model via penalized maximum likelihood.
#' @import glmnet
#' @importFrom graphics plot
#' @importFrom stats coef
#' @param mat Gene expression matrix, columns are cells and rows are genes.
#' @param cell_labels cell clustering labels.
#' @param alpha The elastic-net mixing parameter.
#' @param lam Determine the lambda that according to cross-validation error - default is lambda.min. lambda.1se is also available.
#' @param family Response type. see \code{\link[glmnet]{glmnet}}.
#' @param type.multinomial If "grouped" then a grouped lasso penalty is used on the multinomial coefficients for a variable.
#' @param nfolds number of folds - default is 10.
#' @param cv.plot Save cross-validation plot.
#' @param ... Additional arguments passed on to \code{\link[glmnet]{glmnet}}.
#' @return Expression matrix with selected genes.
#' @export
glm.enet <- function(mat,
                     cell_labels,
                     alpha = 0.5,
                     lam = "lambda.min",
                     family = "multinomial",
                     type.multinomial = "grouped",
                     nfolds = 10,
                     cv.plot = FALSE,
                     ...){
  
  if (ncol(mat) != length(cell_labels)){
    stop("Cell number is not consistent with the length of cell_labels")
  }
  
  cvfit=cv.glmnet(x = t(mat),
                  y = cell_labels, 
                  alpha = alpha,
                  family = family,
                  type.multinomial = type.multinomial,
                  nfolds = nfolds,
                  ...)
  
  if (cv.plot){
    pdf(paste0(nfolds, "_cross_validation.pdf"))
    plot(cvfit)
    dev.off()
  }

  # select beta
  beta <- coef(cvfit, s = lam)  #lambda.1se
  gset <- beta$`1`@i
  print(length(gset))
  wholeSet <- as.numeric(names(table(gset))[-1]) # remove intercept
  
  return(mat[wholeSet,])
  
}

#' Clustering genes to k clusters using affinity propagation algorithm.
#' 
#' Runs affinity propagation clustering.
#' @importFrom apcluster negDistMat corSimMat apclusterK
#' @importFrom nnet class.ind
#' @param mat Gene expression matrix, columns are cells and rows are genes.
#' @param k Desired number of gene clusters.
#' @param metric Metric used to compute the paired distance - default is pearson correlation distance.
#' @param negDist The dissimilarity matrix given by user.
#' @param cell_labels cell clustering labels, it must be given when check_outlier_pattern = TRUE.
#' @param check_outlier_pattern Detect and remove negtive expression pattern dominated by rare cells.
#' @param check_min_size Check the correlation of the cell cluster smaller than check_min_size with all gene clusters.
#' @param check_cor_thr Correlation threshold to remove outlier pattern.
#' @param ... Additional arguments passed on to apclusterK.
#' @return Gene clustering results.
#' @export
apclusK <- function(mat,
                    k = 7,
                    metric = "pearson",
                    negDist = NULL,
                    cell_labels = NULL,
                    check_outlier_pattern = FALSE,
                    check_min_size = 5,
                    check_cor_thr = -0.5,
                    ...){

  mat <- as.data.frame(mat)
  if (is.null(negDist)){
    if (metric == "pearson"){
      negDist <- corSimMat(mat, method = metric)
    }
    else {
      negDist <- negDistMat(mat, method = metric)
    }
  }
  
  orig <- rownames(mat)
  
  apres <- apclusterK(negDist,K=k,...)
  clres <- data.frame(idx = unlist(apres@clusters))
  
  exem_idx <- rep(FALSE, dim(clres)[1])
  exem_idx[rownames(clres) %in% names(apres@exemplars)] <- TRUE
  labels <- rep(1:length(apres@exemplars),unlist(lapply(apres@clusters,length)))
  
  clres$cluster <- labels
  clres$exemplars <- exem_idx
  
  if (check_outlier_pattern){
    if (is.null(cell_labels) || length(cell_labels)!=dim(mat)[2]){
      stop("Cell cluster labels are missing or wrong when check_outlier_pattern")
    }

    # sort the gene expression matrix
    mat <- mat[rownames(clres),]
    
    # I want to remove the gene cluster that negatively correlate with rare cells 
    emexp <- mat[names(apres@exemplars),]
    
    dum <- class.ind(cell_labels)
    clsize <- colSums(dum)
    
    # find the clusters that with size smaller than check_min_size
    check.clus <- as.numeric(which(clsize <= check_min_size))
    
    # Gene cluster that correlated with Cell cluster
    ExCorWithDum <- cor(x=t(emexp),y=dum)[,check.clus]
    outlierIdx <- rowSums(ExCorWithDum < check_cor_thr) > 0
    outlierG <- names(outlierIdx[outlierIdx==TRUE])
    target.cluster <- clres[outlierG,"cluster"]
    fil_mask <- clres$cluster %in% target.cluster == FALSE
    mat.fil <- mat[fil_mask,]
    clres.fil <- clres[fil_mask,]
    new_cl <- clres.fil$cluster
    clres.fil$cluster <- rep(1:length(unique(new_cl)),table(new_cl))
    
    cat("AP clustering and filter outlier expression pattern:")
    print(table(clres.fil$cluster))
    orig <- orig[orig %in% rownames(clres.fil)]
    
    return(clres.fil[orig,])
    
  }
  
  else{
    cat("AP clustering results:")
    print(table(clres$cluster))
    
    return(clres[orig,])
    
  }
}

#' Calculate the p-value of correlation coefficient.
#' 
#' Calculate the significant level of correlation coefficient.
#' @importFrom stats cor pt p.adjust
#' @importFrom grDevices pdf dev.off
#' @param x A numeric vector, matrix or data frame.
#' @param y NULL (default) or a vector, matrix or data frame with compatible dimensions to x. The default is equivalent to y = x (but more efficient).
#' @param method A character string indicating which correlation coefficient (or covariance) is to be computed. One of "pearson" (default), "kendall", or "spearman": can be abbreviated.
#' @param correct Correction method for p-value, see \code{\link[stats]{p.adjust}}.
#' @param ... Additional arguments passed on to \code{\link[stats]{cor}}.
#' @return The matrix of the p-value of correlation coefficient.
#' @export
p.cor <- function(x,
                  y = NULL,
                  method = "pearson",
                  correct = "BH",
                  ...){
  
  p <- as.integer(ncol(x))
  if (p < 1) 
    stop("must have >1 column")
  n <- as.integer(nrow(x))
  if (n < 5) 
    stop("must have >4 observations")
  n <- nrow(x)
  corr <- cor(x, y, method = method)
  p <- matrix(2 * (1 - pt(abs(corr) * sqrt(n - 2)/sqrt(1 - corr * corr), n - 2)), ncol = ncol(corr))
  dimnames(p) <- dimnames(corr)
  # p[p == 0] <- NA
  p.adj <- apply(p,2,function(x) p.adjust(x,method = correct,...))
  
  return(list(r=corr, p=p, p.adj=p.adj))
  
}

#' Heatmap for gene expression pattern.
#' 
#' Plot the clustering result of cells and genes.
#' @importFrom grDevices colorRampPalette
#' @importFrom pheatmap pheatmap
#' @importFrom utils head
#' @param mat Expression matrix, columns are cells and rows are genes.
#' @param filename Filename to save your heatmap.
#' @param vmax The maximum value to plot.
#' @param vmin The minimum value to plot.
#' @param log log2 tranform the data.
#' @param scale Scale the data, result in mean = 0 and std = 1 of each gene.
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
#' @param anno_row Data frame that specifies the annotations shown on left side of the heatmap. Each row defines the features for a specific row. The rows in the data and in the annotation are matched using corresponding row names. Note that color schemes takes into account if variable is continuous or discrete.
#' @param anno_col Similar to annotation_row, but for columns.
#' @param plot_gaps_row Make gaps on row.
#' @param plot_gaps_col Make gaps on column.
#' @param re_write_gene_labels Re-write gene labels to an ordered sequeence.
#' @param re_write_cell_labels Re-write cell labels to an ordered sequeence.
#' @param ... Additional arguments passed on to pheatmap.
#' @return Ordering cluster results.
#' @export
cheatmap <- function(mat,
                     filename = NA,
                     vmax = 2,
                     vmin = -1, 
                     log = FALSE,
                     scale = FALSE,
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
                     anno_row = NA,
                     anno_col = NA,
                     plot_gaps_row = TRUE,
                     plot_gaps_col = TRUE,
                     re_write_gene_labels = FALSE,
                     re_write_cell_labels = FALSE,
                     ...){
  
  if (!is.numeric(vmax) || !is.numeric(vmin)){
    stop("The expression boundary must be numeric")
  }
  if (length(gene_order) == 1){
    gene_order <- as.numeric(strsplit(as.character(gene_order),"")[[1]])
  }
  if (length(cell_order) == 1){
    cell_order <- as.numeric(strsplit(as.character(cell_order),"")[[1]])
  }
  if (!is.null(gene_order) && !is.null(gene_labels) && length(intersect(gene_order, unique(gene_labels)))!=max(gene_labels)){
    stop("The elements of gene_order must be consistent with gene_labels: ",unique(gene_labels))
  }
  if (!is.null(cell_order) && !is.null(cell_labels) && length(intersect(cell_order, unique(cell_labels)))!=max(cell_labels)){
    stop("The elements of cell_order must be consistent with cell_labels: ",unique(cell_labels))
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
      if (re_write_gene_labels){
        gene_labels <- re_order(gene_labels, gene_order)
      }
    }
    else {
      makeGeneOrder <- order(gene_labels)
      if (plot_gaps_row){
        gaps_row <- head(cumsum(table(gene_labels)),-1)
      }
    }
    mat <- mat[makeGeneOrder,]
    output_gene_labels <- data.frame(labels=gene_labels[makeGeneOrder],row.names = rownames(mat))
    if (is.na(anno_row)){
      anno_row = data.frame(Module=factor(output_gene_labels$labels, levels = unique(output_gene_labels$labels)),
                            row.names = rownames(mat))
    }
  }
  else{
    gaps_row <- NULL
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
      if (re_write_cell_labels){
        cell_labels <- re_order(cell_labels, cell_order)
      }
    }
    else {
      makeCellOrder <- order(cell_labels)
      if (plot_gaps_col){
        gaps_col <- head(cumsum(table(cell_labels)),-1)
      }
    }
    mat <- mat[,makeCellOrder]
    output_cell_labels <- data.frame(labels=cell_labels[makeCellOrder],row.names = colnames(mat))
    if (is.na(anno_col)){
      anno_col <- data.frame(Cluster = factor(output_cell_labels$labels, levels = unique(output_cell_labels$labels)), 
                             row.names = colnames(mat)) # bestcells = bestcells EGFP=EGFP
    }
  }
  else{
    gaps_col <- NULL
    output_cell_labels <- NULL
  }
  
  if (log){
    mat <- log2(mat + 1)
  }
  if (scale){
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
           annotation_row = anno_row,
           annotation_col = anno_col,
           gaps_col = gaps_col,
           gaps_row = gaps_row,
           color = cmap,
           cluster_rows = cluster_rows,
           cluster_cols = cluster_cols,
           labels_row = psuGene,
           ...)
  
  return(list(heat=mat,cell_labels=output_cell_labels,gene_labels=output_gene_labels))
  
}

#' Re-ordering clusters.
#' 
#' Re-ordering clusters for heatmap.
#' @param ori_labels Origin labels.
#' @param re_labels Re-ordered labels.
#' @return New labels after ordering.
#' @export
alter_label <- function(ori_labels,re_labels){
  if (length(intersect(unique(ori_labels),unique(re_labels)))!=length(unique(ori_labels))){
    stop("The given labels are not consistent with the origin labels")
  }
  in_order=seq(length(ori_labels))
  names(in_order)=ori_labels
  tre_labelspe=sort(unique(ori_labels))
  
  l=list()
  for(i in 1:length(tre_labelspe)){
    l[[i]]=in_order[which(names(in_order) %in% tre_labelspe[i])]
  }
  
  new_order=c()
  for(i in 1:length(re_labels)){
    new_order=c(new_order,l[[re_labels[i]]])
  }
  
  return(new_order)
  
}

#' Re-write clustering labels.
#' 
#' Re-write clustering labels.
#' @param labs Clustering labels.
#' @param ord_id Cluster order idx.
#' @return A vector of re-write labels.
#' @export
#' 
re_order <- function(labs,ord_id){
  new_ord <- 1:max(labs)
  names(new_ord) <- ord_id
  new_labs <- new_ord[as.character(labs)]
  
  return(new_labs)
  
}

#' Find clustering order using TSP algorithm.
#' 
#' Determine the order of by solving travel saleman problem.
#' @importFrom stats aggregate dist
#' @importFrom TSP as.TSP insert_dummy solve_TSP cut_tour
#' @param mat Expression matrix, columns are cells and rows are genes.
#' @param labels Cluster labels on the row of mat.
#' @param reverse Reverse the path.
#' @param metric Metric used to compute the paired distance - default is pearson correlation distance.
#' @param method Method to solve the TSP, see \code{\link[TSP]{solve_TSP}}.
#' @param seed Random seed.
#' @return Path indice of the clusters.
#' @export
find_order <- function(mat,
                       labels,
                       reverse = FALSE,
                       metric = "pearson",
                       method = "farthest_insertion",
                       seed = 1){
  if (nrow(mat) != length(labels)){
    stop("The rows are not consistent with the length of labels")
  }
  mat <- as.data.frame(mat)
  mat$labels <- labels
  m <- aggregate(. ~ labels, data = mat, FUN = "mean")[,-1]
  set.seed(seed)
  if (metric == "pearson"){
    tsp <- insert_dummy(as.TSP(1-cor(t(m))), label = "cut")
  }
  else{
    tsp <- insert_dummy(as.TSP(dist(m, method = metric)), label = "cut")
  }
  tour <- solve_TSP(tsp, method = method)
  path <- as.vector(cut_tour(tour, "cut"))
  if (reverse) path <- rev(path)
  return(path)
}

#' Gene ontology enrichment analysis of gene clusters.
#' 
#' Perform GO analysis using topGO.
#' @import org.Mm.eg.db
#' @import ggplot2
#' @import ggthemes
#' @import GO.db
#' @import topGO
#' @importFrom AnnotationDbi as.list select keys
#' @importFrom methods new
#' @importFrom utils write.csv
#' @param gene_labels Data frame of gene names with cluster labels.
#' @param method Test statistics to deal with the GO graph structure - default is weight01, a mixture between the elim and weight algorithms. Alternative algorithms are classic, elim, weight, lea, parentchild. See \code{\link[topGO]{topGO-package}} for more detail.
#' @param save_dir Where GO results to be saved. If NA results will be save in the folder called 'GO_results'.
#' @param nodeSize The nodes with less than nodeSize annotated genes are removed from the GO hierarchy
#' @param showNodes How many term to show in the bar plot.
#' @param width Width of the barplot.
#' @param height Height of the barplot.
#' @return go_enrich returns a list of: 
#' \itemize{
#'   \item the statistic results of each clusters;
#'   \item the genes in each term of each cluster.
#' }
#' @export
go_enrich <- function(gene_labels,
                      method = "weight01",
                      save_dir = NA,
                      nodeSize = 30, 
                      showNodes = 15,
                      topNodes=50,
                      width = 12,
                      height = 10){
  if (is.na(save_dir)){
    save_dir = "GO_results"
  }
  if (!dir.exists(save_dir)){
    dir.create(save_dir)
  }
  colnames(gene_labels) <- "labels"
  cl_genes <- rownames(gene_labels)
  keys <- keys(org.Mm.eg.db,"ENSEMBL")
  anno <- select(org.Mm.eg.db, keys=keys, columns = c("SYMBOL"),keytype = "ENSEMBL")
  anno <- anno[!duplicated(anno$ENSEMBL) & !duplicated(anno$SYMBOL),]
  GOID2TERM <- as.list(GOTERM)
  
  cl_go <- list()
  stat_res <- list()
  for (cl in unique(gene_labels$labels)){
    gsym <- cl_genes[gene_labels == cl]
    
    sub_ensg <- anno$ENSEMBL[!is.na(match(anno$SYMBOL, gsym))]
    geneList <- factor(as.integer(anno$ENSEMBL %in% sub_ensg))
    names(geneList) <- anno$ENSEMBL
    
    #construct a topGOdata object
    GOdata <- new("topGOdata",
                  ontology = "BP",
                  allGenes = geneList,
                  geneSel = topDiffGenes,
                  nodeSize = nodeSize, 
                  annot = annFUN.org,
                  mapping = "org.Mm.eg.db",
                  ID = "Ensembl")
    Fish.stat <- new(paste0(method,"Count"),testStatistic=GOFisherTest,name="Fisher test") #elim algorithm is more conservation than classic 
    resultFisher <- getSigGroups(GOdata,Fish.stat)
    allRes <- GenTable(GOdata, Fis = resultFisher, topNodes = topNodes)
    # pvalFis <- score(resultFisher)
    # allRes <- GenTable(GOdata,Fis=resultFisher,orderBy="elimFis",ranksOf="Fis",topNodes=50)
    
    whole_term <- c()
    for (id in allRes$GO.ID){
      whole_term <- c(whole_term, GOID2TERM[[id]]@Term)
    }
    allRes$Term <- whole_term
    stat_res[[as.character(cl)]] <- allRes
    write.csv(allRes, file.path(save_dir, paste0("gene_cl", cl, ".csv")))
    
    #make bar plot for significant terms
    topRes <- allRes[1:showNodes,]
    topRes$Fis[topRes$Fis == "< 1e-30"] <- 1e-30
    goRes <- data.frame(term = paste0(topRes$Term," (", topRes$GO.ID, ")"),
                        pval = topRes$Fis,
                        logpval = -log10(as.numeric(topRes$Fis)))
    goRes$term <- factor(goRes$term, levels = rev(goRes$term))
    g <- ggplot(goRes, aes(term, logpval))
    g + geom_col() + coord_flip() + theme_gdocs() +
      theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + 
      labs(title = paste0("Gene cluster ",cl), y = expression(paste(-log["10"]," p-value"), x = "Terms"))
    ggsave(file.path(save_dir,paste0("gene_cl",cl,".pdf")),width = width, height = height)
    
    #Pick out the genes annotated to GO terms
    cat("Find genes in each term of cluster ",cl," ...\n")
    topGOid <- allRes$GO.ID
    if (!is.null(topGOid)){
      GeneInTerms <- list()
      row.names(anno) <- anno$ENSEMBL
      for (term in topGOid){
        ann.genes <- genesInTerm(GOdata,term)
        GeneInTerm <- intersect(sub_ensg, ann.genes[[1]])
        GeneSymbol <- anno[GeneInTerm,"SYMBOL"]
        names(GeneSymbol) <- GeneInTerm
        GeneInTerms[[term]] <- GeneSymbol
      }
      cl_go[[as.character(cl)]] <- GeneInTerms
    }
  } 
  
  return(list(stat_res=stat_res,cl_go=cl_go))
  
}

#' Fetch genes in GO term.
#' 
#' Fetch ensembl ids of specific GO term id.
#' @import org.Mm.eg.db
#' @import org.Hs.eg.db
#' @param term GO term id, forexample "GO:0031386".
#' @param species Mouse or human database to used. Mouse: mMus, Human: Hs
#' @return Ensembl id of specifc GO term
#' @export
FetchGO = function (term, species = "mMus"){
  if (!(species %in% c("mMus", "Hs"))) {
    stop("'species' needs to be either 'mMus' or 'Hs'")
  }
  if (species == "mMus"){
    xxGO <- AnnotationDbi::as.list(org.Mm.egGO2EG)
    x <- org.Mm.egENSEMBL
  }
  else{
    xxGO <- AnnotationDbi::as.list(org.Hs.egGO2EG)
    x <- org.Hs.egENSEMBL
  }
  BPEG <- unlist(xxGO[term])
  mapped_genes <- AnnotationDbi::mappedkeys(x)
  xxE <- as.list(x[mapped_genes])
  ens_ids <- unlist(xxE[BPEG])
  
  return(ens_ids)
  
}

#' Differential expreesed analysis of two conditions.
#' 
#' DE analysis based on DESeq2.
#' @importFrom DESeq2 DESeqDataSetFromMatrix estimateSizeFactors DESeq results
#' @importFrom BiocGenerics counts
#' @param mat Expression matrix.
#' @param cond Vector indicate the conditions of cells.
#' @param qval Significant level to select genes.
#' @param contrast How to compare.
#' @return data.frame of DESeq2 results.
#' @export
de <- function(mat, cond, qval = 0.1, contrast=NULL){
  
  if (ncol(mat) != length(cond)){
    stop("number of cells is not equal to the length of condition")
  }
  if (length(unique(cond)) > 2){
    stop("de only consider comparison of two conditions")
  }
  if (is.null(contrast)){
    contrast <- c("condition", unique(cond))
  }
  else if (length(contrast) != 2){
    stop("a least two condition should be given")
  }
  else if (!all(contrast %in% cond)){
    stop("contrast is not in condition")
  }
  else{
    contrast <- c("condition", contrast)
  }
  countData = round(mat)
  colData <- data.frame(condition = cond,row.names = colnames(countData))
  
  # construct the DESeqDataSet object from the matrix of counts
  # and the sample information table
  dds <- DESeqDataSetFromMatrix(countData = countData,
                                colData = colData,
                                design = ~ condition)
  # Pre-filtering the dataset
  dds <- dds[ rowSums(counts(dds)) > 1, ]
  dds <- estimateSizeFactors(dds)
  # Differential expression analysis (General Linear Model)
  dds <- DESeq(dds)
  # You may want to specify the comparison of any two levels of condition
  cat("contrast conditions: ", contrast[1], " vs ", contrast[2])
  res <- results(dds, contrast=contrast, alpha = qval)
  
  return(res[complete.cases(res),])
  
}

#' Violin plot.
#' 
#' Violin plot for specific gene in expression matrix.
#' @import ggplot2
#' @import ggthemes
#' @param mat Expression matrix.
#' @param gene Gene id in rownames of mat.
#' @param level Level in colnames of mat.
#' @param levels Factor, custom level define by user.
#' @param color Vector of rgb color indices for violin.
#' @param log Boolen value, log transform the expression or not.
#' @param jitter Boolen value, plot dots or not.
#' @return ggplot2 object
#' @export
violinPlot <- function(mat, gene, level=NULL, levels=NULL, color = NULL, log=TRUE, jitter=FALSE){
  if (!gene %in% rownames(mat)){
    stop("Gene ", gene, " is not in mat")
  }
  if (is.null(levels)){
    mask <- gsub(paste0("(",paste0(level, collapse = "|"),").*"),"\\1",colnames(mat))
    mat <- mat[,mask %in% level]
    days <- factor(mask[mask %in% level], levels = level)
  }
  else{
    days <- levels
  }
  if (log){
    df <- data.frame(exp = log2(as.numeric(mat[gene,]) + 1), level = days)
  }
  else{
    df <- data.frame(exp = as.numeric(mat[gene,]), level = days)
  }
  # xticklabs <- levels(days)
  p <- ggplot(df, aes(factor(level), exp))
  p <- p + geom_violin(aes(fill = factor(days)),
                       scale = "width",
                       trim = TRUE, 
                       adjust = 1) + 
           labs(title = paste0(gene," Expression"),
                fill="Days",
                x = "Days", 
                y = expression(paste(log["2"]," Expression"))) +
           theme_stata(scheme = "s1mono") +
           theme(axis.ticks = element_blank(),
                 axis.title.x = element_text(margin = margin(t = 8)),
                 axis.title.y = element_text(margin = margin(r = 10)))
  
  if (jitter){
    p <- p + geom_jitter(height = 0,
                         width = 0.2)
  }
  if (!is.null(color)){
    p <- p + scale_colour_manual(values = color)
  }
  else{
    p <- p + scale_fill_ptol()
  }
  
  p
  
}
