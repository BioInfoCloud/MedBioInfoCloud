
#' RNAseqDataTrans
#'
#' @param data a data.frame or matrix for raw count or FPKM of RNAseq
#' @param tfun tfun only support 'Counts2FPKM','FPKM2TPM' or 'Counts2TPM'
#' @param species 'hsa' or 'mus'
#' @param gtype 'Ensemble' or 'Symbol'
#'
#' @return data.frame
#' @export RNAseqDataTrans
#'
RNAseqDataTrans <- function(data,tfun,species,gtype){
  if(!(species %in% c("hsa","mus"))){
    message("The value of species should be 'hsa' or 'mus'")
    stop()
  }else if(!(gtype %in% c("Ensembl","Symbol"))){
    message("Gene types(gtype) only support 'Ensemble' or 'Symbol'")
    stop()
  }else if(!(tfun %in% c("Counts2FPKM","FPKM2TPM","Counts2TPM"))){
    message("tfun only support 'Counts2FPKM','FPKM2TPM' or 'Counts2TPM'")
    stop()
  }
  if(species == 'hsa'){
    conid = intersect(rownames(data),hsaGeneInfo[,gtype])
  }else{
    conid = intersect(rownames(data),musGeneInfo[,gtype])
  }
  if(length(conid)>0){
    data = data[conid,]
    effLen = hsaGeneInfo[match(conid,hsaGeneInfo[,gtype]),"effLen"]
  }
  if(tfun %in% c("Counts2FPKM",'Counts2TPM','FPKM2TPM')){
    if(tfun == "Counts2FPKM"){
      data <- apply(data, 2, Counts2FPKM, effLen = effLen)
    }else if(tfun == 'Counts2TPM'){
      data <- apply(data, 2, Counts2TPM, effLen = effLen)
    }else if(tfun == 'FPKM2TPM'){
      data <- apply(data,2,FPKM2TPM)
    }
    return(data)
  }else{stop("tfun only support 'Counts2FPKM','FPKM2TPM' or 'Counts2TPM'")}

}


#' A function for differential expression analysis
#'
#' @param data for matrix input: a matrix of non-negative integers
#' @param group for matrix input: a DataFrame or data.frame with at least a single column. Rows of group correspond to columns of countData
#' @param comparison A string linked by "-" that represents the grouping information in the group. Such as "tumor-normal".
#' @param cutFC >1
#' @param cutFDR 0 < cutFDR < 0.05
#' @param method One of DESeq2, edgeR, and limma.
#'
#' @return A data.frame
#' @export geneDEAnalysis
#'
geneDEAnalysis <- function (data, group, comparison,
                              method = "DESeq2",cutFC=2,cutFDR=0.05, filter = TRUE){
  dge = DGEList(counts = data)
  keep <- rowSums(cpm(dge) > 1) >= 0.5 * length(group)
  if (method == "DESeq2") {
    coldata <- data.frame(group)
    dds <- DESeqDataSetFromMatrix(countData = data, colData = coldata, design = ~group)
    dds$group <- factor(dds$group, levels = rev(strsplit(comparison, "-", fixed = TRUE)[[1]]))
    if (filter == TRUE) {
      dds <- dds[keep, ]
    }
    dds <- DESeq(dds)
    res <- results(dds)
    DEGAll <- data.frame(res)
    colnames(DEGAll) <- c("baseMean", "logFC", "lfcSE", "stat", "PValue", "FDR")
  }
  else if (method %in% c("edgeR", "limma")) {
    group <- factor(group)
    design <- model.matrix(~0 + group)
    colnames(design) <- levels(group)
    contrast.matrix <- makeContrasts(contrasts = comparison,
                                     levels = design)
    if (filter == TRUE) {
      dge <- dge[keep, , keep.lib.sizes = TRUE]
    }
    dge <- calcNormFactors(dge)
    if (method == "edgeR") {
      dge <- estimateDisp(dge, design)
      fit <- glmFit(dge, design)
      lrt <- glmLRT(fit, contrast = contrast.matrix)
      DEGAll <- lrt$table
      DEGAll$FDR <- p.adjust(DEGAll$PValue, method = "fdr")
    }
    else if (method == "limma") {
      v <- voom(dge, design = design, plot = FALSE)
      fit <- lmFit(v, design)
      fit2 <- contrasts.fit(fit, contrast.matrix)
      fit2 <- eBayes(fit2)
      DEGAll <- topTable(fit2, coef = 1, n = Inf)
      colnames(DEGAll) <- c("logFC", "AveExpr", "t", "PValue", "FDR", "B")
    }
  }
  DEGAll$symbol <- rownames(DEGAll)
  DEGAll <- dplyr::select(DEGAll, symbol, everything())

  DEGAll <- DEGAll  %>%
    dplyr::mutate(direction = factor(ifelse(FDR < cutFDR & abs(logFC) > cutFC,#添加direction一列
                                            ifelse(logFC > cutFC, "Up", "Down"),"Ns"),
                                     levels=c('Up','Down','Ns')))
  return(DEGAll)
}

#' arrayDataDEA_limma
#'
#' @param data matrix
#' @param group a vector
#' @param comparison character
#' @param cutFC >1
#' @param cutFDR 0< cutFDR <0.05
#'
#' @return data.frame
#' @export arrayDataDEA_limma
arrayDataDEA_limma <- function(data, group,comparison,cutFC = 2,cutFDR = 0.05){
  glist <- group %>% factor(., levels = unique(group), ordered = F)
  head(glist)
  glist <- model.matrix(~factor(glist)+0)  #把group设置成一个model matrix
  colnames(glist) <- unique(group)
  df.fit <- lmFit(data, glist)  ## 数据与list进行匹配
  df.matrix <- makeContrasts(comparison, levels = glist)
  fit <- contrasts.fit(df.fit, df.matrix)
  fit <- eBayes(fit)
  tempOutput <- topTable(fit,n = Inf, adjust = "fdr")
  colnames(tempOutput) <- c("logFC", "AveExpr","t","PValue","FDR","B")
  tempOutput <- tempOutput  %>%
    dplyr::mutate(direction = factor(ifelse(FDR < cutFDR & abs(logFC) > cutFC,#添加direction一列
                                            ifelse(logFC > cutFC, "Up", "Down"),"Ns"),
                                     levels=c('Up','Down','Ns')))
  return(tempOutput)
}
