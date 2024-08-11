
#' getGeneLenFromeGTF
#'
#' @param gtf The gtf parameter represents the path to a GTF (Gene Transfer Format) format file.
#'
#' @return data.frame
#' @export GenomicFeatures,GenomicRanges
#'
getGeneLenFromeGTF <- function(gtf){
  # reference 1: https://mp.weixin.qq.com/s/dwbpJ0nhzyIp9fDv7fEWEQ
  # reference 2: https://mp.weixin.qq.com/s/lazavD3jzRVO4QkxysHQcg
  # 加载必要的包
  # library(GenomicFeatures)
  # library(GenomicRanges)
  txdb <- GenomicFeatures::makeTxDbFromGFF(gtf,format="gtf")
  exons.list.per.gene <- GenomicFeatures::exonsBy(txdb, by = "gene")
  exonic.gene.sizes <- lapply(exons.list.per.gene,
                              function(x){sum(GenomicRanges::width(GenomicRanges::reduce(x)))})
  eff_length <- do.call(rbind,lapply(exonic.gene.sizes, data.frame))
  eff_length <- data.frame(Ensembl = do.call(rbind,strsplit(rownames(eff_length),'\\.'))[,1],
                           effLen = eff_length[,1])
  return(eff_length)
}


#' getGeneTypeInfoFromeGTF
#'
#' @param gtf The gtf parameter represents the path to a GTF (Gene Transfer Format) format file.
#'
#' @return data.frame
#' @export data.table
#'
getGeneTypeInfoFromeGTF <- function(gtf){
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
#' @param gtf The gtf parameter represents the path to a GTF (Gene Transfer Format) format file.
#'
#' @return data.frame
#' @export data.table,GenomicFeatures,GenomicRanges
#'
getGeneBaseInfo <- function(gtf){
  ens2symInfo <- getGeneInfoFromeGTF(gtf)
  eff_length <- getGeneLenFromeGTF(gtf)
  GeneInfo <- merge(ens2symInfo,eff_length,by = "Ensembl")
  return(GeneInfo)
}

#' RNAseqDataConversion
#'
#' @param data RNAseq expression data, either a matrix or a data frame
#' @param type must be one of "Counts2TPM", "Counts2FPKM", or "FPKM2TPM"
#' @param species must be one of "homo", "mus", or NULL.
#' @param gtf The gtf parameter represents the path to a GTF (Gene Transfer Format) format file.
#'
#' @return either a matrix or a data frame
#' @export dplyr
#'
RNAseqDataConversion <- function(data,type,species = "homo",gtf = NULL){
  if(type %in% c("Counts2TPM","Counts2FPKM","FPKM2TPM")){
    if(type == "FPKM2TPM"){
      return(FPKM2TPM(data))
    }else if(is.null(species) & file.exists(gtf)){
      ano <- getGeneBaseInfo(gtf)
    }else if(species == "homo"){
      ano <- MedBioInfoCloud::hsaGeneInfo
    }else if(species == "mus"){
      ano <- MedBioInfoCloud::musGeneInfo
    }
    ano <- dplyr::arrange(ano,Symbol,desc(effLen))
    ano <- ano[!duplicated(ano$Symbol),]
    rownames(ano) <- ano$Symbol
    congene <- intersect(rownames(data),rownames(ano))
    if(!is.numeric(congene)){
      ano <- ano[congene,]
      data <- data[congene,]
      if(type == "Counts2TPM"){
        return(Counts2TPM(counts = data,effLen = ano$effLen))
      }else if(type == "Counts2FPKM"){
        return(Counts2FPKM(counts = data,effLen = ano$effLen))
      }
    }else{message("Gene annotation information does not match the row names of the data.")}

  }
}




#' Counts2TPM
#'
#' @param counts a data.frame or matrix for raw count of RNAseq.
#' @param effLen The length of genes.
#'
#' @return a data.frame
#' @export Counts2TPM
#'
Counts2TPM <- function(counts, effLen){
  rate <- log(counts) - log(effLen)
  denom <- log(sum(exp(rate)))
  return(exp(rate - denom + log(1e6)))
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

#' outputGmtFile
#'
#' @param input list or data.frame
#' @param description NA
#' @param folder For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param filename filepath
#' @return NULL
#' @export
#'

outputGmtFile <- function(input,description = NA,folder= ".",filename){
  ifelse(dir.exists(folder),"",dir.create(folder,recursive = T))
  if(is.data.frame(input)){
    nms <- unique(input[,1])
    ls <- lapply(nms, function(nm){
      input[input[,1]== nms[nm],2]
    })
    names(ls) <- nms
    input = ls
  }
  if(is.list(input)){
    output <- file(paste0(folder,"/",filename), open="wt")
    nms <- names(input)
    lapply(1:length(nms),function(x){
      if(!is.na(description)){
        dsc = description[x]
      }else{dsc = description}
      outlines = paste0(c(nms[x], dsc, input[[x]]),collapse='\t')
      writeLines(outlines, con=output)
    })
    close(output)
  }
}

#' geneset2gmt
#'
#' @param geneset A gene vector or a txt file path. Note that the contents of the txt text file should be one gene in a row.
#' @param genesetname A string meaning the name of the gene set.
#' @param description A string used to introduce the gene set, default is NA.
#' @param return Return a value of "data frame" or "GeneSetCollection", said what kind of value, the equivalent of respectively with clusterProfiler: : read. GMT () function and GSEABase: : getGmt () function to read in the result of the GMT file, You can return the result as required.
#' @param folder For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param filename For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#'
#' @return "data frame" or "GeneSetCollection" object.
#' @export geneset2gmt
#'
geneset2gmt <- function(geneset,
                        genesetname,
                        description = NA,
                        return = "data.frame",
                        folder= ".",filename){
  if(length(geneset) == 1 ){
    if(grep(".txt$",geneset) & file.exists(geneset)){
      geneset <- readLines(geneset)
    }
  }
  if(is.vector(geneset) & length(geneset) > 1){
    input <- list(geneset = geneset)
    names(input) <- genesetname
  }
  if(exists("input") & grep("gmt$",filename)){
    outputGmtFile(input = input,description = description,folder,filename = filename)
    if(return == "data.frame"){
      gs <- clusterProfiler::read.gmt(paste0(folder,"/",filename))
    }else if(return == "GeneSetCollection"){
      gs <- GSEABase::getGmt(paste0(folder,"/",filename))
    }
    return(gs)
  }else{stop("Error:geneset or filename")}
}


#' tidy.gmt
#'
#' @param filepath For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param fun "stat" or "merge"
#' @param Source For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param termName For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param addTotal For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param save For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param folder For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param filename For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#'
#' @return data.frame
#' @export clusterProfiler
#'
tidy.gmt <-  function (filepath, fun = "stat", Source = "", termName = NULL,
                       addTotal = FALSE, save = TRUE, folder = ".", filename = "geneset") {
  ifelse(dir.exists(folder), "", dir.create(folder, recursive = T))
  if (dir.exists(filepath)) {
    fp <- dir(filepath, pattern = "gmt$", full.names = T)
  }else if (is.vector(filepath) & sum(grep("gmt$", filepath)) ==
            length(filepath)) {
    fp <- filepath
  }
  stat.gs <- data.frame()
  mergedf <- data.frame()
  gs <- c()
  for (x in fp) {
    ldf <- clusterProfiler::read.gmt(x)
    if(nrow(ldf) != 0){
      mergedf <- rbind(mergedf, ldf)
      gs <- unique(c(gs, ldf$gene))
      df <- data.frame(term = gsub("_", " ", as.character(unique(ldf$term))),
                       Source = Source, `Gene count` = nrow(ldf))
      stat.gs <- rbind(stat.gs, df)
    }

  }
  totallgene <- data.frame(term = "", Source = Source, `Gene count` = paste0(length(unique(gs)),
                                                                             "(unique)"))
  stat.gs <- rbind(stat.gs, totallgene)
  if (!is.null(termName)) {
    colnames(stat.gs)[1] <- termName
  }
  if (fun == "stat") {
    if (save == TRUE) {
      writeLines(unique(gs), con = paste0(folder, "/",
                                          filename, "-stat.unique.txt"))
      write.csv(stat.gs, file = paste0(folder, "/", filename,
                                       "-stat.csv"))
    }
    return(stat.gs)
  }
  else if (fun == "merge") {
    if (save == TRUE) {
      writeLines(unique(gs), con = paste0(folder, "/",
                                          filename, "-merge.unique.txt"))
      outputGmtFile(input = mergedf, description = NA,
                    folder = folder, filename = paste0(filename,
                                                       "-merge.gmt"))
      if (addTotal == TRUE) {
        addmergedf <- data.frame(term = ifelse(is.null(termName),
                                               "Total gene", termName), gene = unique(gs))
        mergedf <- rbind(mergedf, addmergedf)
        outputGmtFile(input = mergedf, description = NA,
                      folder = folder, filename = paste0(filename,
                                                         "-add-unique.gmt"))
      }
    }
    return(mergedf)
  }
}

# tidy.gmt(filepath = "G:/publicData/base_files/GeneSet/Cytoskeleton",
#          fun = "merge",
#          folder = "G:/myProject/Cytoskeleton/data/geneset")




#' tidyGene.fromeReactome
#'
#' @param filepath For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param fun For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param Source For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param termName For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param save For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param folder For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param filename For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#'
#' @return data.frame
#' @export stringr
#'
tidyGene.fromeReactome <- function(filepath
                                   ,fun = "stat"
                                   ,Source = "Reactome"
                                   ,addTotal = TRUE
                                   ,termName=NULL
                                   ,save = TRUE
                                   ,folder = "."
                                   ,filename = "geneset"){
  ifelse(dir.exists(folder),"",dir.create(folder,recursive = T))
  if(dir.exists(filepath)){
    fp <- dir(filepath,pattern = "tsv$",full.names = T)
  }else if(is.vector(filepath) & sum(grep("tsv$",filepath)) == length(filepath)){
    fp <- filepath
  }
  stat.gs <- data.frame()
  mergedf <- data.frame()
  gs <- c()
  for(x in fp){
    if(Source == "Reactome"){
      rc <- read.csv(x,header = T,sep = "\t")[,c("MoleculeType","MoleculeName")]
      rc <- tidyr::separate(rc,MoleculeName,c("UniProt","gene")," ")
      rc <- rc[rc$MoleculeType == "Proteins",]
      rc <- rc[!duplicated(rc$gene),]
      filenm <- unlist(strsplit(x,"/"))
      filenm <- filenm[length(filenm)]
      # 提取[号前的内容
      # str_extract(filenm, "^[^\\[]+")
      term = stringr::str_extract(filenm,"^[^\\[]*")
      PathwayID = stringr::str_extract(filenm, "(?<=\\[)[^\\]]+")
    }else{
      message("开发中")
    }
    gs <- c(gs,rc$gene)
    statdf <- data.frame(term = term
                         ,Source = Source
                         ,PathwayID = PathwayID
                         ,`Gene count` = nrow(rc))
    stat.gs <- rbind(stat.gs,statdf)
    gmtdf <- data.frame(term = rep(term,nrow(rc))
                        ,gene = rc$gene)
    mergedf <- rbind(mergedf,gmtdf)
  }
  totallgene <- data.frame(term = ""
                           ,Source = ""
                           ,PathwayID =""
                           ,`Gene count` = paste0(length(unique(gs)),"(unique)")
  )
  stat.gs <- rbind(stat.gs,totallgene)
  if (!is.null(termName)) {
    colnames(stat.gs)[1] <- termName
  }

  if(fun == "stat"){
    if(save == TRUE){
      writeLines(unique(gs),con = paste0(folder,"/",filename,"-stat.unique.txt"))
      write.csv(stat.gs,file = paste0(folder,"/",filename,"-stat.csv"))
    }
    return(stat.gs)
  } else if(fun == "merge"){
    if(save == TRUE){
      outputGmtFile(input = mergedf
                    ,description = NA
                    ,folder= folder
                    ,filename = paste0(filename,"-merge.gmt"))
      writeLines(unique(gs),con = paste0(folder,"/",filename,"-merge.unique.txt"))
      if(addTotal == TRUE){
        addmergedf <- data.frame(term = ifelse(is.null(termName),"Total gene",termName)
                                 ,gene = unique(gs))
        mergedf <- rbind(mergedf,addmergedf)
        outputGmtFile(input = mergedf
                      ,description = NA
                      ,folder= folder
                      ,filename = paste0(filename,"-add-unique.gmt"))
      }
    }
    return(mergedf)
  }
}

#' read.gmt.to.getGmt
#'
#' @param genesetdf The result of the read.gmt() function
#'
#' @return GeneSetCollection object.
#' @export GSEABase,read.gmt.to.getGmt
#'
read.gmt.to.getGmt <- function(genesetdf){
  # 首先，按照 term 分组
  gene_sets_list <- split(genesetdf$gene, genesetdf$term)
  # 创建 GeneSet 对象列表，过滤掉空字符串
  gene_sets <- lapply(names(gene_sets_list), function(term) {
    genes <- gene_sets_list[[term]]
    genes <- genes[genes != ""]  # 去掉空字符串
    GeneSet(
      geneIds = genes,
      setName = term,
      setIdentifier = term
    )
  })
  # 创建 GeneSetCollection 对象
  gsc <- GeneSetCollection(gene_sets)
  return(gsc)
}


#' GeneSetCollection.to.df
#'
#' @param GeneSetCollection GeneSetCollection object.
#'
#' @return data.frame
#' @export GeneSetCollection.to.df
#'
GeneSetCollection.to.df <- function(GeneSetCollection){
  genesetname <- names(GeneSetCollection)
  geneset = lapply(genesetname, function(x){
    data.frame(term = x,gene = GeneSetCollection[[x]]@geneIds)
  })
  geneset <- do.call(rbind,geneset)
  return(geneset)
}



