#' Counts2TPM
#'
#' @param counts a data.frame or matrix for raw count of RNAseq.
#' @param effLen The length of genes.
#'
#' @return a data.frame
#' @export
#'
Counts2TPM <- function(counts, effLen){
  rate <- log(counts) - log(effLen)
  denom <- log(sum(exp(rate)))
  exp(rate - denom + log(1e6))
}
#' Counts2FPKM
#'
#' @param counts a data.frame or matrix for raw count of RNAseq.
#' @param effLen The length of genes.
#'
#' @return a data.frame
#' @export
#'
Counts2FPKM <- function(counts, effLen){
  N <- sum(counts)
  exp( log(counts) + log(1e9) - log(effLen) - log(N) )
}
#' FPKM2TPM
#'
#' @param fpkm a data.frame or matrix for fpkm of RNAseq.
#'
#' @return a data.frame
#' @export
#'
FPKM2TPM <- function(fpkm){
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}


#' RNAseqDataTrans
#'
#' @param data a data.frame or matrix for raw count or FPKM of RNAseq
#' @param tfun tfun only support 'Counts2FPKM','FPKM2TPM' or 'Counts2TPM'
#' @param species 'hsa' or 'mus'
#' @param gtype 'Ensemble' or 'Symbol'
#'
#' @return
#' @export
#'
#' @examples
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
  if(tfun == "Counts2FPKM"){
    data <- apply(data, 2, Counts2FPKM, effLen = effLen)
  }else if(tfun == 'Counts2TPM'){
    data <- apply(data, 2, Counts2TPM, effLen = effLen)
  }else{
    data <- apply(data,2,FPKM2TPM)
  }
  return(data)
}


