
#' getGeneLenFromeGTF
#'
#' @param gtf gtf Gene annotation file in gtf format. https://www.gencodegenes.org/
#'
#' @return a data.fram
#' @seealso getGeneInfoFromeGTF,getGeneBaseInfo
#' @export getGeneLenFromeGTF
getGeneLenFromeGTF <- function(gtf){
  # reference 1: https://mp.weixin.qq.com/s/dwbpJ0nhzyIp9fDv7fEWEQ
  # reference 2: https://mp.weixin.qq.com/s/lazavD3jzRVO4QkxysHQcg
  txdb <- GenomicFeatures::makeTxDbFromGFF(gtf,format="gtf")
  exons.list.per.gene <- GenomicFeatures::exonsBy(txdb, by = "gene")
  exonic.gene.sizes <- lapply(exons.list.per.gene,
                              function(x){sum(width(reduce(x)))})
  eff_length <- do.call(rbind,lapply(exonic.gene.sizes, data.frame))
  eff_length <- data.frame(Ensembl = do.call(rbind,strsplit(rownames(eff_length),'\\.'))[,1],
                           effLen = eff_length[,1])
  return(eff_length)
}

#' getGeneInfoFromeGTF
#'
#' @param gtf gtf gtf Gene annotation file in gtf format. https://www.gencodegenes.org/
#'
#' @return a data.frame
#' @export getGeneInfoFromeGTF
#'
getGeneInfoFromeGTF <- function(gtf){
  if (is.character(gtf)) {
    if(!file.exists(gtf)) stop("Bad gtf file.")
    message("Treat gtf as file")
    gtf = data.table::fread(gtf, header = FALSE)
  } else {
    data.table::setDT(gtf)
  }
  gtf = gtf[gtf[[3]] == "gene", ]
  pattern_id = ".*gene_id \"([^;]+)\";.*"
  pattern_name = ".*gene_name \"([^;]+)\";.*"
  pattern_type = ".*gene_type \"([^;]+)\";.*"
  gene_id = sub(pattern_id, "\\1", gtf[[9]])
  gene_id = do.call(rbind,strsplit(gene_id,'\\.'))[,1]
  gene_name = sub(pattern_name, "\\1", gtf[[9]])
  gene_type = sub(pattern_type, "\\1", gtf[[9]])

  EnsemblTOGenename <- data.frame(Ensembl = gene_id,
                                  Symbol = gene_name,
                                  gene_type = gene_type,
                                  stringsAsFactors = FALSE)
  return(EnsemblTOGenename)
}

#' getGeneBaseInfo
#'
#' @param gtf gtf gtf Gene annotation file in gtf format. https://www.gencodegenes.org/
#'
#' @return a data.frame
#' @export getGeneBaseInfo
#'
getGeneBaseInfo <- function(gtf){
  ens2symInfo <- getGeneInfoFromeGTF(gtf)
  eff_length <- getGeneLenFromeGTF(gtf)
  GeneInfo <- merge(ens2symInfo,eff_length,by = "Ensembl")
  return(GeneInfo)
}

