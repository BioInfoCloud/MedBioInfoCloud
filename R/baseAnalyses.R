#' Get random values from permutation test,  the mean of which as baseline to calculate FUN score.
#'
#' @param input.list A list of variables containing scaled gene expression matrix, giving number of bins.
#' @param genes.dist.bins A matrix of giving expression bins within genes, which were binned according to their average expression across cells or samples.
#' @param b.sign A logical value to indicate the overlapped features between target features and the features derived from gene expressio matrix.
#' @param num.rounds  A integer value to indicate the permutation times of iteration; 1000 by default and 10000 will be better for reproducibility.
#'
#' @return Continuous numerical variables with values, named by cells or samples.
#' @export get_random_FUN
#'
#' @examples
#'
#' \dontrun{
#' r.scores <- get_random_FUN(input.list, input.list$genes.dist.bins, b.sign, num.rounds = num.rounds)
#' }
#'
get_random_FUN <- function(input.list, genes.dist.bins, b.sign, num.rounds = 1000){
  sign.bins <- as.matrix(table(genes.dist.bins[b.sign]))
  q <-rownames(sign.bins)
  bg <- matrix(data = F, nrow = length(genes.dist.bins), ncol = num.rounds)
  for (i in 1:nrow(sign.bins)){
    num.genes <- sign.bins[i]
    if(num.genes > 0){
      idx <- which(is.element(genes.dist.bins, q[i]))
      for (j in 1:num.rounds){
        idxj <- sample(idx, num.genes)
        bg[idxj, j] <- T
      }
    }
  }
  r.scores <- apply(bg, 2, function(x)colMeans(input.list$expr.scaled[x,]))
  r.scores <- rowMeans(r.scores)
  return(r.scores)
}



#' scoringSys
#'
#' @param expr TPM matrix derived from scRNAseq, or normalized gene expression matrix from bulk samples; rows are genes, cols are cells or samples.
#' @param geneset For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param TCGA For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param method For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param study.type A string value to indicate the study type for calculation. Allowed values contain c('scRNAseq', 'bulk_RNAseq').
#' @param save For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param folder For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#'
#' @return data.frame
#' @export scoringSys
#'
scoringSys <- function(expr,geneset
                       ,TCGA = FALSE,
                       method ="ssgsea",
                       filtergenenum = 10,
                       study.type = 'bulk_RNAseq',
                       save = TRUE,folder = "."){
  ifelse(dir.exists(folder),"",dir.create(folder,recursive = T))
  if(TCGA == TRUE){
    # expr <- expr[["tpm"]]
    expr<- filterGeneTypeExpr(expr =expr,
                               fil_col = "gene_type",
                               filter = "protein_coding")
    expr <- splitTCGAmatrix(data = expr,sample = "Tumor")
    expr <- delTCGA_dup_sample(expr,col_rename = TRUE)
  }
  if(is.data.frame(geneset)){
    geneset <- read.gmt.to.getGmt(gs)
  }else if(grep(".gmt$",geneset) & file.exists(geneset)){
    geneset <- GSEABase::getGmt(geneset,
                                geneIdType=GSEABase::SymbolIdentifier())
  }
  if(method %in% c("ssgsea","gsva","zscore","plage")){
    score = GSVA::gsva(as.matrix(expr),
                 geneset,
                 method = method,
                 kcdf='Gaussian',
                 abs.ranking=TRUE)
  }else if(method == "CRDscore"){
    genename <- names(geneset)
    score <- lapply(genename, function(x){
      calFUNscore(expr = expr
                  ,fungeneset = geneset[[x]]@geneIds
                  ,study.type = study.type
                  )
    })
    score <- do.call(rbind,score)
    rownames(score) <- genename
  }else{stop("")}
  if(save == TRUE){
    ifelse(dir.exists(folder),"",dir.create(folder,recursive = T))
    write.csv(score,file = paste0(folder,"/",method,"-score.csv"))
  }
  return(score)
}




