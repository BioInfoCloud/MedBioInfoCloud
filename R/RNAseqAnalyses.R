
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
#' @param method One of DESeq2, edgeR, and limma.
#' @param filter Boolean value. Whether to filter out genes with low expression.
#'
#' @return A data.frame
#' @export geneDEAnalysis
#'
geneDEAnalysis <- function (data, group, comparison,
                            method = "DESeq2", filter = FALSE){
  dge = edgeR::DGEList(counts = data)
  keep <- rowSums(edgeR::cpm(dge) > 1) >= 0.5 * length(group)
  if (method == "DESeq2") {
    coldata <- data.frame(group)
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = data, colData = coldata, design = ~group)
    dds$group <- factor(dds$group, levels = rev(strsplit(comparison, "-", fixed = TRUE)[[1]]))
    if (filter == TRUE) {
      dds <- dds[keep, ]
    }
    dds <- DESeq2::DESeq(dds)
    res <- DESeq2::results(dds)
    DEGAll <- data.frame(res)
    DEGAll$FDR <- p.adjust(DEGAll$pvalue, method = "fdr")
    colnames(DEGAll) <- c("baseMean", "log2FC", "lfcSE", "stat", "pValue","pAdj", "FDR")
  }
  else if (method %in% c("edgeR", "limma")) {
    group <- factor(group)
    design <- model.matrix(~0 + group)
    colnames(design) <- levels(group)
    contrast.matrix <- limma::makeContrasts(contrasts = comparison,levels = design)
    if (filter == TRUE) {
      dge <- dge[keep, , keep.lib.sizes = TRUE]
    }
    dge <- edgeR::calcNormFactors(dge)
    if (method == "edgeR") {
      dge <- edgeR::estimateDisp(dge, design)
      fit <- edgeR::glmFit(dge, design)
      lrt <- edgeR::glmLRT(fit, contrast = contrast.matrix)
      DEGAll <- lrt$table
      DEGAll$pAdj <- p.adjust(DEGAll$PValue, method = "BH")
      DEGAll$FDR <- p.adjust(DEGAll$PValue, method = "fdr")
      colnames(DEGAll) <- c("log2FC", "logCPM", "LR", "pValue","pAdj", "FDR")
    }
    else if (method == "limma") {
      v <- limma::voom(dge, design = design, plot = FALSE)
      fit <- limma::lmFit(v, design)
      fit2 <- limma::contrasts.fit(fit, contrast.matrix)
      fit2 <- limma::eBayes(fit2)
      DEGAll <- limma::topTable(fit2, coef = 1, n = Inf)
      DEGAll$FDR <- p.adjust(DEGAll$P.Value, method = "fdr")
      colnames(DEGAll) <- c("log2FC", "AveExpr", "t", "pValue", "pAdj", "B","FDR")
    }
  }
  DEGAll$symbol <- rownames(DEGAll)
  DEGAll <- dplyr::select(DEGAll, symbol, everything())
  return(DEGAll)
}

#' arrayDataDEA_limma
#'
#' @param data matrix
#' @param group a vector
#' @param comparison character
#'
#' @return data.frame
#' @export arrayDataDEA_limma
arrayDataDEA_limma <- function(data, group,comparison){
  glist <- group %>% factor(., levels = unique(group), ordered = F)
  head(glist)
  glist <- limma::model.matrix(~factor(glist)+0)  #把group设置成一个model matrix
  colnames(glist) <- unique(group)
  df.fit <- limma::lmFit(data, glist)  ## 数据与list进行匹配
  df.matrix <- limma::makeContrasts(comparison, levels = glist)
  fit <- limma::contrasts.fit(df.fit, df.matrix)
  fit <- limma::eBayes(fit)
  DEGAll <- limma::topTable(fit, coef = 1, n = Inf)
  DEGAll$FDR <- p.adjust(DEGAll$P.Value, method = "fdr")
  colnames(DEGAll) <- c("log2FC", "AveExpr", "t", "pValue", "pAdj", "B","FDR")
  return(DEGAll)
}

#' Screening for differentially expressed genes
#'
#' @param DEAR a data.frame. The results of differential expression analysis,
#'   must contain the column specified by `pMethod` (e.g., "FDR", "pAdj", "pValue")
#'   and a "log2FC" column.
#' @param cutFC The threshold absolute value for log2 fold change.
#'   The default is 1.
#' @param pMethod Metrics for statistical testing: "FDR", "pAdj" or "pValue".
#' @param cutP Significance screening threshold. Default value: 0.05
#' @param gene_col 基因名列名，默认"label"
#' @param top_n 顶部基因的数量（上调和下调分别标记此数量），默认3
#' @returns data.frame with an additional "direction" column indicating
#'   Up/Down/Ns (Not significant)
#' @export selectDEG
#'
selectDEG <- function(DEAR, pMethod = "FDR",cutFC = 1,cutP = 0.05,
                      gene_col = "symbol", top_n = 0) {
  # 检查必要列
  required_cols <- c(pMethod, "log2FC", gene_col)
  if (!all(required_cols %in% colnames(DEAR))) {
    stop("缺少必要列：", paste(setdiff(required_cols, colnames(DEAR)), collapse = ", "))
  }

  # 检查top_n是否为正整数
  if (!is.numeric(top_n) || top_n < 0 || top_n != as.integer(top_n)) {
    stop("top_n必须是正整数（如0、1、2、3...）")
  }

  # 第一步：添加差异方向（Up/Down/Ns）
  df <- DEAR %>%
    dplyr::mutate(
      direction = dplyr::case_when(
        .data[[pMethod]] < cutP & .data[["log2FC"]] >= cutFC ~ "Up",
        .data[[pMethod]] < cutP & .data[["log2FC"]] <= -cutFC ~ "Down",
        TRUE ~ "Ns"
      )
    )

  # 第二步：按log2FC排序标记顶部基因（修复管道符）
  if(top_n != 0){
    # 上调基因：log2FC从大到小排序，取前top_n个
    top_up <- df %>%
      dplyr::filter(direction == "Up") %>%
      dplyr::arrange(dplyr::desc(log2FC)) %>%  # 修复：添加管道符
      dplyr::slice_head(n = top_n) %>%
      dplyr::pull({{gene_col}})

    # 下调基因：log2FC从小到大排序（最负的在前），取前top_n个
    top_down <- df %>%
      dplyr::filter(direction == "Down") %>%
      dplyr::arrange(log2FC) %>%  # 修复：添加管道符
      dplyr::slice_head(n = top_n) %>%
      dplyr::pull({{gene_col}})

    # 第三步：添加label列
    df %>%
      dplyr::mutate(
        label = dplyr::case_when(
          .data[[gene_col]] %in% top_up ~ .data[[gene_col]],
          .data[[gene_col]] %in% top_down ~ .data[[gene_col]],
          TRUE ~ ""
        )
      )
  }else{df$label = ""}
  return(df)

}
