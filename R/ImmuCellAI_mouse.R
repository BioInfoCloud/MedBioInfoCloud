#' Predict Immune Cell Abundance in Mouse Samples Using Hierarchical Strategy
#'
#' This function predicts the abundance of multiple immune cell types (24 subtypes in 3 layers) in mouse transcriptome data,
#' simulating flow cytometry's hierarchical classification process. It supports group comparison and custom reference files.
#'
#' @param sample Gene expression profile matrix. Rows are gene symbols (mouse), columns are sample IDs.
#'               Accepts FPKM/TPM format for RNA-seq or log2-transformed signal for microarray.
#'               If `group_tag = 1`, the first row must be group labels (e.g., "control" vs "tumor").
#' @param data_type Character string specifying the type of input data. Must be one of "rnaseq" (RNA-seq data) or "microarray" (microarray data).
#' @param group_tag Numeric indicator for group comparison. 1 = perform group comparison (requires group labels in the first row of `sample`), 0 = no group comparison.
#' @param customer Numeric indicator for custom reference files. 1 = use user-provided reference files, 0 = use built-in reference files (default).
#' @param sig_file List of custom gene signatures (R list format). Each element is a gene vector corresponding to a cell type. Required only when `customer = 1`.
#' @param exp_file Custom reference expression matrix (tab-separated). Rows are genes, columns are cell types. Required only when `customer = 1`.
#'
#' @return A list with two components:
#'         \item{abundance}{Data frame of immune cell abundance. Rows are samples, columns are 24 immune cell types + total infiltration score. Values are rounded to 4 decimal places.}
#'         \item{group_result}{Data frame of group comparison results (only if `group_tag = 1`). Rows are groups + "p value", columns are immune cell types + infiltration score.}
#' @export
#'
#' @examples
#' # Load example data (built-in after package installation)
#' data("example_mouse_expr") # Example RNA-seq TPM matrix (rows: genes, columns: samples)
#'
#' # Run prediction (no group comparison, use built-in references)
#' result <- ImmuCellAI_mouse(
#'   sample = example_mouse_expr,
#'   data_type = "rnaseq",
#'   group_tag = 0,
#'   customer = 0
#' )
#'
#' # View cell abundance results
#' head(result$abundance)
ImmuCellAI_mouse = function (sample, data_type, group_tag, customer, sig_file = NULL,
                             exp_file = NULL)
{
  data("l1_marker")
  data("l2_marker")
  data("l3_marker")
  data("marker_exp")
  data("l1_cell_correction_matrix_new")
  data("l2_cell_correction_matrix_new")
  data("l3_cell_correction_matrix_new")
  T_sub_ls <- c("CD4_T_cell", "CD8_T_cell", "NKT", "Tgd")
  B_sub_ls <- c("B1_cell", "Follicular_B", "Germinal_center_B",
                "Marginal_Zone_B", "Memory_B", "Plasma_cell")
  DC_sub_ls <- c("cDC1", "cDC2", "MoDC", "pDC")
  Gra_sub_ls <- c("Basophil", "Eosinophil", "mast_cell", "Neutrophils")
  Macro_sub_ls <- c("M1_macrophage", "M2_macrophage")
  CD4_T_sub_ls <- c("CD4_Tm", "Naive_CD4_T", "T_helper_cell",
                    "Treg")
  CD8_T_sub_ls <- c("CD8_Tc", "CD8_Tcm", "CD8_Tem", "CD8_Tex",
                    "Naive_CD8_T")
  group_fre <- c()
  marker_exp_raw <- marker_exp
  group_index = 0
  if (group_tag) {
    group_column <- sample[1, ]
    group_content <<- sample[1, ]
    sample = sample[-1, ]
  }
  sam = apply(sample, 2, as.numeric)
  row.names(sam) = row.names(sample)
  sam_exp = as.matrix(sam)
  colnames(sam_exp) = colnames(sample)
  sample_exp <- sam_exp
  paper_marker <- l1_marker

  #修改后的代码===========
  ssgsea_param <- GSVA::ssgseaParam(as.matrix(marker_exp_raw[, names(paper_marker)]), paper_marker)
  ref_pre <- GSVA::gsva(ssgsea_param)

  # 原来的代码
  # ref_pre <- gsva(as.matrix(marker_exp_raw[, names(paper_marker)]),
  #     paper_marker, method = "ssgsea", ssgsea.norm = TRUE)

  cell_cor_new_mat <- matrix(rep(0, length(names(l1_marker)) *
                                   length(names(l1_marker))), ncol = length(names(l1_marker)))
  row.names(cell_cor_new_mat) = names(l1_marker)
  colnames(cell_cor_new_mat) = names(l1_marker)
  diag(cell_cor_new_mat) <- 1
  sample = as.matrix(sample_exp)
  sam = apply(sample, 2, as.numeric)
  row.names(sam) = row.names(sample_exp)
  tt = intersect(row.names(sam), as.vector(unlist(paper_marker)))
  genes = intersect(tt, row.names(marker_exp))
  sam_exp = as.matrix(sam[genes, ])
  row.names(sam_exp) = genes
  colnames(sam_exp) = colnames(sample)
  tt = intersect(row.names(sam_exp), as.vector(unlist(paper_marker)))
  genes = intersect(tt, row.names(marker_exp))
  marker_exp_partial <- marker_exp[genes, ]
  marker_tag_mat = c()
  for (cell in names(paper_marker)) {
    tag = marker_tag(genes, as.vector(unlist(paper_marker[cell])))
    marker_tag_mat = cbind(marker_tag_mat, tag)
  }
  row.names(marker_tag_mat) = row.names(marker_exp_partial)
  colnames(marker_tag_mat) = names(paper_marker)
  marker_reccur <- marker_times_count(paper_marker)
  gene_weight_mt <- gene_ref_weight_mt(paper_marker, marker_reccur)
  gene_weight_mt <- gene_weight_mt[genes, ]
  immune_deviation_sample = apply(sam_exp, 2, function(x) sample_ratio(x,
                                                                       marker_exp_partial[, names(paper_marker)], marker_tag_mat,
                                                                       data_type, gene_weight_mt, ref_pre, cell_cor_new_mat))
  infil_marker <- genes
  infil_mt <- c()
  for (cell in names(l1_marker)[-c(5, 6)]) {
    genes_new <- intersect(l1_marker[[cell]], row.names(immune_deviation_sample))
    infil_marker_cell <- list()
    infil_exp <- data.frame(as.vector(unlist(immune_deviation_sample[genes_new,
    ])))
    row.names(infil_exp) <- as.vector(unlist(sapply(colnames(immune_deviation_sample),
                                                    function(x) paste(x, genes_new, sep = "_"))))
    colnames(infil_exp) <- "Ratio"
    for (sam in colnames(immune_deviation_sample)) {
      infil_marker_cell[[sam]] <- row.names(infil_exp)[grep(sam,
                                                            row.names(infil_exp))]
    }


    # 修改代码
    infil_score <- GSVA::gsva(GSVA::ssgseaParam(as.matrix(infil_exp), infil_marker_cell))
    # 老代码
    # infil_score <- gsva(as.matrix(infil_exp), infil_marker_cell,
    #     method = "ssgsea", ssgsea.norm = TRUE)
    infil_mt <- cbind(infil_mt, infil_score)
  }
  infil_mt[which(infil_mt < 0)] = 0
  infil_score <- apply(infil_mt, 1, function(x) sum(x))
  tt = intersect(row.names(sample_exp), as.vector(unlist(paper_marker)))
  genes = intersect(tt, row.names(marker_exp))
  sam_exp <- sample_exp[genes, ]
  layer1_pre <- getResult(sam_exp, data_type, marker_exp_raw[,
                                                             names(paper_marker)], paper_marker, ref_pre, cell_cor_new_mat,
                          l1_cell_correction_matrix_new)
  layer1_result <- round(t(apply(layer1_pre, 1, function(x) x/sum(x))),
                         4)
  marker_exp <- marker_exp_raw
  paper_marker <- l2_marker
  cell_cor_mat = matrix(rep(0, length(names(paper_marker)) *
                              length(names(paper_marker))), ncol = length(names(paper_marker)))
  row.names(cell_cor_mat) = names(paper_marker)
  colnames(cell_cor_mat) = names(paper_marker)

  #修改后的代码===========
  ssgsea_pre_param <- GSVA::ssgseaParam(as.matrix(marker_exp_raw[, names(paper_marker)]), paper_marker)
  ref_pre <- GSVA::gsva(ssgsea_pre_param)
  # 旧代码
  # ref_pre <- gsva(as.matrix(marker_exp_raw[, names(paper_marker)]),
  #     paper_marker, method = "ssgsea", ssgsea.norm = TRUE)

  tt = intersect(row.names(sample_exp), as.vector(unlist(paper_marker)))
  genes = intersect(tt, row.names(marker_exp))
  sam_exp <- sample_exp[genes, ]
  layer2_pre <- getResult(sam_exp, data_type, marker_exp_raw[,
                                                             names(paper_marker)], paper_marker, ref_pre, cell_cor_new_mat,
                          l2_cell_correction_matrix_new)
  layer2_pre_trans <- data.frame(t(layer2_pre), check.names = F)
  tt <- t(layer_norm(layer2_pre_trans, T_sub_ls))
  tmp <- round(t(apply(tt, 1, function(x) x/sum(x))), 4)
  layer2_T_norm <- tmp * layer1_result[, "T_cell"]
  tt <- t(layer_norm(layer2_pre_trans, B_sub_ls))
  tmp <- round(t(apply(tt, 1, function(x) x/sum(x))), 4)
  layer2_B_norm <- tmp * layer1_result[, "B_cell"]
  tt <- t(layer_norm(layer2_pre_trans, DC_sub_ls))
  tmp <- round(t(apply(tt, 1, function(x) x/sum(x))), 4)
  layer2_DC_norm <- tmp * layer1_result[, "Dendritic_cells"]
  tt <- t(layer_norm(layer2_pre_trans, Gra_sub_ls))
  tmp <- round(t(apply(tt, 1, function(x) x/sum(x))), 4)
  layer2_Gra_norm <- tmp * layer1_result[, "Granulocytes"]
  tt <- t(layer_norm(layer2_pre_trans, Macro_sub_ls))
  tmp <- round(t(apply(tt, 1, function(x) x/sum(x))), 4)
  layer2_macro_norm <- tmp * layer1_result[, "Macrophage"]
  layer2_cell_all <- cbind(layer2_T_norm, layer2_B_norm, layer2_DC_norm,
                           layer2_Gra_norm, layer2_macro_norm)
  marker_exp <- marker_exp_raw
  paper_marker <- l3_marker

  #修改后的代码===========
  ssgsea_pre_param <- GSVA::ssgseaParam(as.matrix(marker_exp_raw[, names(paper_marker)]), paper_marker)
  ref_pre <- GSVA::gsva(ssgsea_pre_param)

  # 旧代码========
  # ref_pre <- gsva(as.matrix(marker_exp_raw[, names(paper_marker)]),
  #     paper_marker, method = "ssgsea", ssgsea.norm = TRUE)

  tt = intersect(row.names(sample_exp), as.vector(unlist(paper_marker)))
  genes = intersect(tt, row.names(marker_exp))
  sam_exp <- sample_exp[genes, ]
  layer3_pre <- getResult(sam_exp, data_type, marker_exp_raw[,
                                                             names(paper_marker)], paper_marker, ref_pre, cell_cor_new_mat,
                          l3_cell_correction_matrix_new)
  layer3_pre_trans <- data.frame(t(layer3_pre), check.names = F)
  tt <- t(layer_norm(layer3_pre_trans, CD4_T_sub_ls))
  tmp <- round(t(apply(tt, 1, function(x) x/sum(x))), 4)
  layer3_CD4_sub <- tmp * layer2_T_norm[, "CD4_T_cell"]
  tt <- t(layer_norm(layer3_pre_trans, CD8_T_sub_ls))
  tmp <- round(t(apply(tt, 1, function(x) x/sum(x))), 4)
  layer3_CD8_sub <- tmp * layer2_T_norm[, "CD8_T_cell"]
  result_final <- cbind(layer1_result, layer2_cell_all, layer3_CD4_sub,
                        layer3_CD8_sub)
  result_final <- round(result_final, 4)
  result_final[grep(TRUE, is.na(result_final))] = 0
  result_mat <- result_final
  if (group_tag) {
    group_name <- sort(unique(as.vector(unlist(group_content))))
    p_va = c()
    group_column <- as.numeric(as.factor(as.vector(unlist(group_content))))
    Infiltration_score <- round(infil_score, 3)
    result_group = cbind(result_mat, Infiltration_score,
                         group_column)
    result_tt = apply(result_group, 2, as.numeric)
    if (length(group_name) > 2) {
      for (cell in colnames(result_group)) {
        result_group_new <- result_group[, c(cell, "group_column")]
        t = aov(group_column ~ ., data.frame(result_group_new))
        p_va = c(p_va, round(summary(t)[[1]][["Pr(>F)"]],
                             2))
      }
    }
    else {
      g1_index = grep(1, group_column)
      g2_index = grep(2, group_column)
      result_mat1 <- cbind(result_mat, Infiltration_score)
      for (cell in colnames(result_mat)) {
        c_ = wilcox.test(result_mat[g1_index, cell],
                         result_mat[g2_index, cell])
        p_va = c(p_va, round(c_$p.value, 2))
      }
    }
    row.names(result_tt) = row.names(result_group)
    result_tt = data.frame(result_tt)
    exp_median = aggregate(. ~ group_column, data = result_tt,
                           median)
    exp_median = rbind(exp_median[, -1], p_va)
    row.names(exp_median) = c(group_name, "p value")
    group_fre <- exp_median
  }
  Infiltration_score <- round(infil_score, 3)
  T_FRE <- cbind(result_mat, Infiltration_score)
  return(list(abundance = T_FRE, group_result = group_fre))
}

#' Calculate Immune Cell Abundance with SSGSEA and Compensation
#'
#' Internal helper function for `ImmuCellAI_mouse`. It computes initial immune cell abundance using SSGSEA,
#' then adjusts the result with a pre-defined compensation matrix to reduce cross-cell interference.
#'
#' @param sample Gene expression profile matrix (rows: genes, columns: samples) after filtering.
#' @param data_type Type of input data. One of "rnaseq" or "microarray".
#' @param marker_exp Reference expression matrix of cell type-specific markers (rows: genes, columns: cell types).
#' @param paper_marker List of cell type-specific marker genes (names = cell types, values = gene vectors).
#' @param ref_pre Pre-computed SSGSEA result of the reference marker expression matrix.
#' @param cell_cor_new_mat Cell correlation matrix (rows/columns = cell types) for initial correction.
#' @param compensation_matrix Pre-defined compensation matrix to adjust final abundance (reduces cross-cell interference).
#'
#' @return Transposed matrix of adjusted immune cell abundance. Rows = cell types, columns = samples.
#' @export
#'
#' @examples
#' # Internal use only; called by ImmuCellAI_mouse()
#' # Load built-in reference data
#' data("l1_marker")
#' data("marker_exp")
#' data("l1_cell_correction_matrix_new")
#'
#' # Pre-compute ref_pre with SSGSEA
#' ssgsea_param <- GSVA::ssgseaParam(as.matrix(marker_exp[, names(l1_marker)]), l1_marker)
#' ref_pre <- GSVA::gsva(ssgsea_param)
#'
#' # Example filtered sample matrix
#' sample_filtered <- marker_exp[1:100, 1:3] # Mock sample data
#'
#' # Run getResult
#' initial_result <- getResult(
#'   sample = sample_filtered,
#'   data_type = "rnaseq",
#'   marker_exp = marker_exp,
#'   paper_marker = l1_marker,
#'   ref_pre = ref_pre,
#'   cell_cor_new_mat = matrix(1, nrow = length(l1_marker), ncol = length(l1_marker)),
#'   compensation_matrix = l1_cell_correction_matrix_new
#' )
getResult = function(sample,data_type,marker_exp,paper_marker,ref_pre,cell_cor_new_mat,compensation_matrix){
  #marker_exp <- read.csv("/project/xiamx/immunecell/data/GEO/main/prepare_for_train_set.txt",sep = "\t",row.names = 1)
  #marker_exp <-
  sample=as.matrix(sample)
  sam = apply(sample,2,as.numeric)
  row.names(sam) = row.names(sample)
  tt = intersect(row.names(sam),as.vector(unlist(paper_marker)))
  genes = intersect(tt,row.names(marker_exp))
  sam=data.frame(sam)
  sam_exp = as.matrix(sam[genes,])
  colnames(sam_exp) = colnames(sample)
  marker_exp = marker_exp[genes,]


  marker_tag_mat = c()
  # binary tag matrix
  for(cell in names(paper_marker)){
    tag = marker_tag(genes,as.vector(unlist(paper_marker[cell])))
    marker_tag_mat = cbind(marker_tag_mat,tag)
  }
  row.names(marker_tag_mat) = row.names(marker_exp)
  colnames(marker_tag_mat) = names(paper_marker)
  marker_reccur<-marker_times_count(paper_marker)
  # gene_weight_mt<-gene_ref_weight_mt(paper_marker,marker_reccur)
  # gene_weight_mt<-gene_weight_all_ls[[num_layer]]
  # genes<-intersect(genes,row.names(gene_weight_mt))
  gene_weight_mt<-c()

  exp_new = apply(sam_exp,2,function(x) sample_ratio(x,marker_exp[,names(paper_marker)],marker_tag_mat,data_type,gene_weight_mt,ref_pre,cell_cor_new_mat))
  #print(head(exp_new))
  #print(paper_marker)
  #修改后的代码===========
  ssgsea_result_param <- GSVA::ssgseaParam(exp_new,paper_marker)
  result <- GSVA::gsva(ssgsea_result_param)
  # 原来的代码
  # result = gsva(exp_new,paper_marker,method="ssgsea",ssgsea.norm=TRUE)
  #exp_raio_new<-apply(sam_exp,2,function(x) sample_ratio_correction(x,marker_exp[,names(paper_marker)],marker_tag_mat,data_type,gene_weight_mt))
  #print(head(result))
  #
  result<-apply(result,1,function(x) round(x-min(x)+min(abs(x)),3))

  # if(ncol(result)<3){
  #   result[which(result<0)]=0
  # }else{
  #   result = result - apply(result,1,min)
  #   #result[which(result<0)]=0
  # }
  compensation_matrix_num = apply(compensation_matrix,2,as.numeric)
  # progress$set(value = 20,detail = "Adjusting result by Compensation matrix")
  # incProgress(0.2, detail = "Immune infiltration calculating")
  # Sys.sleep(0.5)
  # shinyWidgets::updateProgressBar(title = 'Immune infiltration calculating',id = "pb2",value=70,session=getDefaultReactiveDomain())
  # if(customer==0){
  row.names(compensation_matrix_num) = row.names(compensation_matrix)
  result_mat = compensation(t(result),compensation_matrix_num)
  if(ncol(result_mat)>1){
    result_mat=apply(result_mat,1,function(x) round(x,3))
  }else{
    result_mat=t(round(result_mat,3))
  }
  return(result_mat)
}

#' Adjust Immune Cell Abundance with Compensation Matrix
#'
#' Internal helper function to adjust raw SSGSEA scores using a compensation matrix,
#' minimizing cross-reactivity between different cell types via constrained least squares.
#'
#' @param raw_score Raw SSGSEA score matrix (rows = samples, columns = cell types).
#' @param compensation_matrix Compensation matrix (rows/columns = cell types) for adjusting cross-cell interference.
#'
#' @return Adjusted abundance matrix. Rows = cell types, columns = samples. Values are non-negative.
#' @export
#'
#' @examples
#' # Example raw SSGSEA scores (2 samples, 3 cell types)
#' raw_score <- matrix(c(0.8, 0.6, 0.2, 0.7, 0.5, 0.3), nrow = 2, ncol = 3,
#'                     dimnames = list(c("Sample1", "Sample2"), c("Tcell", "Bcell", "DC")))
#'
#' # Example compensation matrix (3x3)
#' comp_matrix <- matrix(c(1, 0.1, 0.05, 0.1, 1, 0.03, 0.05, 0.03, 1), nrow = 3, ncol = 3,
#'                       dimnames = list(colnames(raw_score), colnames(raw_score)))
#'
#' # Run compensation
#' adjusted_score <- compensation(raw_score, comp_matrix)
compensation = function(raw_score,compensation_matrix){
  raw_score=as.matrix(raw_score)
  compensation_matrix = compensation_matrix
  diag(compensation_matrix) = 1
  rows <- rownames(raw_score)[rownames(raw_score) %in%  rownames(compensation_matrix)]
  #print(rows)
  if(ncol(raw_score)==1){
    scores <- as.matrix(pracma::lsqlincon(compensation_matrix[rows,rows], raw_score, lb = 0))
  }else{
    scores <- apply(raw_score[rows,], 2, function(x) pracma::lsqlincon(compensation_matrix[rows,rows], x, lb = 0))

  }
  scores<-apply(scores,1,function(x) round(x-min(x)+min(abs(x)),3))
  #scores[scores < 0] = 0
  colnames(scores) <- rows
  return(t(scores))
}

#' Normalize Immune Cell Abundance by Cell Subset
#'
#' Internal helper function to filter and normalize abundance of specific cell subsets (e.g., T cell subtypes)
#' based on a pre-defined cell list.
#'
#' @param pre_result Pre-computed abundance matrix (rows = cell types, columns = samples).
#' @param cell_ls Character vector of target cell types to filter and normalize.
#'
#' @return Filtered and normalized abundance matrix (rows = target cell types, columns = samples).
#' @export
#'
#' @examples
#' # Example pre-result matrix (5 cell types, 2 samples)
#' pre_result <- matrix(c(0.3, 0.2, 0.1, 0.05, 0.35, 0.25, 0.2, 0.15, 0.1, 0.3),
#'                      nrow = 5, ncol = 2,
#'                      dimnames = list(c("CD4_Tcell", "CD8_Tcell", "Bcell", "DC", "Macrophage"),
#'                                      c("Sample1", "Sample2")))
#'
#' # Target T cell subsets
#' t_cell_ls <- c("CD4_Tcell", "CD8_Tcell")
#'
#' # Run layer_norm
#' normalized_tcell <- layer_norm(pre_result, t_cell_ls)
layer_norm<-function(pre_result,cell_ls){
  pre_result%>%
    dplyr::mutate(cellType=row.names(.))%>%
    dplyr::filter(cellType%in%cell_ls)->sub_fra
  row.names(sub_fra)<-sub_fra$cellType
  sub_fra<-sub_fra[,-ncol(sub_fra)]
  return(sub_fra)
}

#' Count Marker Gene Occurrences Across Cell Types
#'
#' Internal helper function to count how many cell types each marker gene is associated with,
#' based on the input marker list.
#'
#' @param paper_marker List of cell type-specific marker genes (names = cell types, values = gene vectors).
#'
#' @return List where each element is a vector of cell types associated with the corresponding gene (name = gene symbol).
#' @export
#'
#' @examples
#' # Example marker list (2 cell types, 3 genes)
#' marker_list <- list(Tcell = c("Cd3e", "Cd4"), Bcell = c("Cd19", "Cd4"))
#'
#' # Count gene occurrences
#' gene_count <- marker_times_count(marker_list)
marker_times_count=function(paper_marker){
  #all_genes=unique(as.vector(unlist(paper_marker)))
  gene_count=list()
  for(cell in names(paper_marker)){
    for(gene in paper_marker[[cell]]){
      gene_count[[gene]]<-c(gene_count[[gene]],cell)
      # if(length(which(names(gene_count)==gene))==0){
      #   gene_count[[gene]]=1
      # }else{
      #   gene_count[[gene]]=gene_count[[gene]]+1
      # }
    }
  }
  return(gene_count)
}

#' Calculate Gene Weight Matrix Based on Reference Expression
#'
#' Internal helper function to compute a gene weight matrix, where weights are proportional to
#' the expression of each marker gene in its associated cell types.
#'
#' @param marker_ls List of cell type-specific marker genes (names = cell types, values = gene vectors).
#' @param marker_reccur List from `marker_times_count()` (gene names as keys, associated cell types as values).
#'
#' @return Gene weight matrix (rows = marker genes, columns = cell types). Values are normalized weights.
#' @export
#'
#' @examples
#' # Load built-in reference expression
#' data("marker_exp")
#'
#' # Example marker list and recurrence count
#' marker_ls <- list(Tcell = c("Cd3e", "Cd4"), Bcell = c("Cd19", "Cd4"))
#' marker_reccur <- marker_times_count(marker_ls)
#'
#' # Compute gene weight matrix
#' weight_mat <- gene_ref_weight_mt(marker_ls, marker_reccur)
gene_ref_weight_mt<-function(marker_ls,marker_reccur){
  gene_weight_mt<-c()
  for(gene in names(marker_reccur)){
    tmp<-rep(1,length(names(marker_ls)))
    names(tmp)<-names(marker_ls)
    cells<-marker_reccur[[gene]]
    if(length(cells)>1){
      cell_exp_sum<-sum(as.numeric(as.vector(unlist(marker_exp[gene,cells]))))
      for(cell in cells){
        tmp[cell]=marker_exp[gene,cell]/cell_exp_sum
      }
    }
    gene_weight_mt<-rbind(gene_weight_mt,tmp)
  }
  row.names(gene_weight_mt)<-names(marker_reccur)
  colnames(gene_weight_mt)<-names(marker_ls)
  return(gene_weight_mt)
}

#' Calculate Correlation Between Cell Type Marker Expressions
#'
#' Internal helper function to compute Pearson correlation coefficients between the expression profiles
#' of different cell types based on their marker genes.
#'
#' @param marker_exp Reference expression matrix (rows = marker genes, columns = cell types).
#' @param marker_ls List of cell type-specific marker genes (names = cell types, values = gene vectors).
#'
#' @return Cell-cell correlation matrix (rows/columns = cell types). Values are Pearson correlation coefficients.
#' @export
#'
#' @examples
#' # Load built-in data
#' data("marker_exp")
#' data("l1_marker")
#'
#' # Compute cell correlation
#' cell_cor <- cell_type_correlation(marker_exp, l1_marker)
cell_type_correlation<-function(marker_exp,marker_ls){
  exp_cor_mt<-c()
  for(c1 in names(marker_ls)){
    cor_va<-c()
    c1_exp<-marker_exp[marker_ls[[c1]],c1]
    for(c2 in names(marker_ls)){
      c2_exp<-marker_exp[marker_ls[[c1]],c2]
      cor_tmp<-cor.test(as.vector(unlist(c1_exp)),as.vector(unlist(c2_exp)))
      cor_va<-c(cor_va,cor_tmp$estimate)
    }
    exp_cor_mt<-rbind(exp_cor_mt,cor_va)
  }
  colnames(exp_cor_mt)<-names(marker_ls)
  row.names(exp_cor_mt)<-names(marker_ls)
  return(exp_cor_mt)
}

#' Calculate Total Immune Infiltration Score
#'
#' Helper function to compute a total immune infiltration score by summing normalized marker gene expressions
#' across all immune cell types.
#'
#' @param exp Gene expression vector for a single sample (names = gene symbols).
#' @param gene_name Character vector of gene symbols corresponding to `exp`.
#' @param data_type Type of input data. One of "rnaseq" or "microarray".
#' @param pre_result Pre-computed abundance result (not used in current implementation).
#'
#' @return Numeric value of total immune infiltration score for the sample.
#' @export
#'
#' @examples
#' # Example gene expression vector (5 marker genes)
#' exp_vec <- c(10, 8, 5, 3, 6)
#' names(exp_vec) <- c("Cd3e", "Cd4", "Cd19", "Cd8a", "Foxp3")
#'
#' # Compute infiltration score
#' inf_score <- immune_infiltate_calculate(
#'   exp = exp_vec,
#'   gene_name = names(exp_vec),
#'   data_type = "rnaseq",
#'   pre_result = NULL
#' )
immune_infiltate_calculate=function(exp,gene_name,data_type,pre_result){
  inf = 0
  names(exp) = gene_name
  for (cell in names(immune_infiltate_marker)){
    abun = 0
    markers = as.vector(unlist(immune_infiltate_marker[cell]))
    for (gene in markers){
      if(data_type == "microarray"){
        abun = abun + as.numeric(exp[gene])/marker_exp_T[gene,cell]
      }else{
        abun = abun + as.numeric(log2(exp[gene]+1))/marker_exp_T[gene,cell]
      }
    }
    inf = inf+abun/length(as.vector(unlist(immune_infiltate_marker[cell])))
  }
  return(inf)
}

#' Create Binary Marker Gene Tag Matrix
#'
#' Internal helper function to generate a binary matrix indicating whether each gene is a marker for each cell type.
#'
#' @param comgenes Character vector of genes to include in the tag matrix (rows of the output).
#' @param tag_gene Character vector of marker genes for a specific cell type.
#'
#' @return Numeric vector (length = length(comgenes)) where 1 = gene is a marker, 0 = not a marker.
#' @export
#'
#' @examples
#' # All genes in the expression matrix
#' all_genes <- c("Cd3e", "Cd4", "Cd19", "Cd8a", "Foxp3")
#'
#' # Marker genes for T cells
#' tcell_markers <- c("Cd3e", "Cd4", "Cd8a")
#'
#' # Create binary tag vector
#' tag_vec <- marker_tag(all_genes, tcell_markers)
marker_tag = function(comgenes,tag_gene){
  a = comgenes
  a[which(comgenes%in%tag_gene)] = 1
  a[which(a!=1)] = 0
  a = as.numeric(a)
  return(a)
}

#' Calculate Sample-Specific Marker Gene Ratio
#'
#' Internal helper function to compute the ratio of sample gene expression to reference marker expression,
#' weighted by a binary tag matrix (indicating marker genes for each cell type).
#'
#' @param data Gene expression vector for a single sample (numeric).
#' @param marker_exp Reference marker expression matrix (rows = genes, columns = cell types).
#' @param marker_tag_mat Binary tag matrix (rows = genes, columns = cell types) from `marker_tag()`.
#' @param data_type Type of input data. One of "rnaseq" or "microarray".
#' @param gene_weight_mt Gene weight matrix from `gene_ref_weight_mt()` (not used in current implementation).
#' @param ref_pre Pre-computed SSGSEA result of reference data (not used in current implementation).
#' @param cell_cor_new_mat Cell correlation matrix (not used in current implementation).
#'
#' @return Numeric vector of normalized gene ratios (length = number of genes in `data`).
#' @export
#'
#' @examples
#' # Example sample data (5 genes)
#' sample_data <- c(12, 9, 4, 7, 5)
#'
#' # Reference marker expression (5 genes, 2 cell types)
#' ref_marker <- matrix(c(10, 8, 0, 6, 0, 0, 0, 10, 0, 8), nrow = 5, ncol = 2,
#'                      dimnames = list(c("G1", "G2", "G3", "G4", "G5"), c("Tcell", "Bcell")))
#'
#' # Binary tag matrix (5 genes, 2 cell types)
#' tag_mat <- cbind(c(1,1,0,1,0), c(0,0,1,0,1))
#'
#' # Calculate sample ratio
#' ratio <- sample_ratio(
#'   data = sample_data,
#'   marker_exp = ref_marker,
#'   marker_tag_mat = tag_mat,
#'   data_type = "rnaseq",
#'   gene_weight_mt = NULL,
#'   ref_pre = NULL,
#'   cell_cor_new_mat = NULL
#' )
sample_ratio = function(data,marker_exp,marker_tag_mat,data_type,gene_weight_mt,ref_pre,cell_cor_new_mat){
  exp = 0
  #print(data[1:2])
  #Eprint(data_type)
  if(data_type == "microarray"){
    for (cell in colnames(marker_exp)){
      # cor_cell<-intersect(names(which(cell_cor_new_mat[cell,]==1)),names(which(ref_pre[,cell]>0)))
      # cor_cell_ratio<-ref_pre[cor_cell,cell]/sum(ref_pre[cor_cell,cell])
      # names(cor_cell_ratio)<-cor_cell
      # for(cc in cor_cell){
      exp<-exp+data/(marker_exp[,cell])*marker_tag_mat[,cell]
    }
    #}
  }else{
    for (cell in colnames(marker_exp)){
      #exp = exp+log2(data+1)*gene_weight_mt[,cell]/marker_exp[,cell]*marker_tag_mat[,cell]
      # cor_cell<-intersect(names(which(cell_cor_new_mat[cell,]==1)),names(which(ref_pre[,cell]>0)))
      # cor_cell_ratio<-ref_pre[cor_cell,cell]/sum(ref_pre[cor_cell,cell])
      # names(cor_cell_ratio)<-cor_cell
      # for(cc in cor_cell){
      exp<-exp+log2(data+1)/(marker_exp[,cell])*marker_tag_mat[,cell]
      #}
    }
  }
  #  print(exp[1:5])
  return(exp)
}
