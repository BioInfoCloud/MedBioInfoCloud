
#' getTCGA_RNAseqData
#'
#' @param project A list of valid project (see list with TCGAbiolinks:::getGDCprojects()$project_id)]
#' @param save TRUE or FALSE
#' @param folder The path specifies the location of the folder where the data is saved when the 'save' parameter is set to TRUE.
#'
#' @return
#' @export
#'
#' @examples
getTCGA_RNAseqData <- function(project,save = FALSE,folder = "."){
  message(paste0("=====================",project," Starting====================="))
  query <- TCGAbiolinks::GDCquery(
    project = project,
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts"
  )
  TCGAbiolinks::GDCdownload(query = query,method = "api",files.per.chunk=10)
  expRda <- TCGAbiolinks::GDCprepare(query,save = F,summarizedExperiment = F)##
  count <- dplyr::select(expRda, matches("^unstranded|gene_id|gene_name|gene_type"))
  colnames(count) = gsub("unstranded_","",colnames(count))
  tpm <- dplyr::select(expRda, matches("^tpm_unstranded|gene_id|gene_name|gene_type"))
  colnames(tpm) = gsub("tpm_unstranded_","",colnames(tpm))
  fpkm <- dplyr::select(expRda, matches("^fpkm_unstranded|gene_id|gene_name|gene_type"))
  colnames(fpkm) = gsub("fpkm_unstranded_","",colnames(fpkm))
  STARdata <- list(count = count,tpm = tpm,fpkm = fpkm)
  if(save == TRUE){
    save(STARdata,file = paste0(folder,"/",project,"-STARdata.Rdata"))
  }
  return(STARdata)
}

#' getProteinExp
#'
#' @param project A list of valid project (see list with TCGAbiolinks:::getGDCprojects()$project_id)]
#' @param save TRUE or FALSE
#' @param folder The path specifies the location of the folder where the data is saved when the 'save' parameter is set to TRUE.
#'
#' @return data.frame
#' @export
#'
#' @examples
getProteinExp <- function(project,save = FALSE,folder = "."){
  message(paste0("=====================",project," Starting====================="))
  query <- TCGAbiolinks::GDCquery(project = project,
                                  data.category = "Proteome Profiling",
                                  data.type = "Protein Expression Quantification")
  GDCdownload(query)
  Proteome_data <- TCGAbiolinks::GDCprepare(query)
  if(save == TRUE){
    save(Proteome_data,file = paste0(folder,"/",project,"-Proteome.Rdata"))
  }
  return(Proteome_data)
}

#' getSNV_Masked_data
#'
#' @param project A list of valid project (see list with TCGAbiolinks:::getGDCprojects()$project_id)]
#' @param save TRUE or FALSE
#' @param folder The path specifies the location of the folder where the data is saved when the 'save' parameter is set to TRUE.
#'
#' @return
#' @export
#'
#' @examples
getSNV_Masked_data <- function(project,save = FALSE,folder = "."){
  query_SNV <- TCGAbiolinks::GDCquery(project = project,
                       data.category = "Simple Nucleotide Variation",
                       data.type = "Masked Somatic Mutation",
                       workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking")
  TCGAbiolinks::GDCdownload(query_SNV)
  snv <- GDCprepare(query_SNV)
  if(save == TRUE){
    save(snv,file = paste0(folder,"/",project,"-SNV_MaskedSomaticMutation.Rdata"))
  }
  return(snv)
}


###用于读入TCGA数据库miRNA数据(Isoform Expression Quantification)
#' dlTCGAmiRNAdata
#'
#' @param data data is a data frame where the first column represents the names of miRNAs, the second column contains the corresponding values, and the third column is the barcode.
#'
#' @return
#' @export
#'
#' @examples
dlTCGAmiRNAdata <- function(data){
  colnames(data) <- c("mirName","value","barcode")
  temp <- lapply(unique(data$barcode), FUN = function(r){
    assign(r,data[data$barcode == r,] %>% dplyr::group_by(mirName) %>%
             dplyr::summarize(value= sum(value)))

  })
  names(temp) <- unique(data$barcode)
  mirdata <- NULL
  for(id in names(temp)){
    data <- temp[[id]]
    colnames(data) <- c("mirName",id)
    if(is.null(mirdata)){
      mirdata <- data
    }else(mirdata <- merge(mirdata,data,by="mirName",all = TRUE))
  }
  rownames(mirdata) <- mirdata$mirName
  mirdata <- mirdata[,-1]
  mirdata[is.na(mirdata)] <- 0
  return(mirdata)
}

#' get_miRNA_IsoformEQ
#'
#' @param project A list of valid project (see list with TCGAbiolinks:::getGDCprojects()$project_id)]
#' @param save TRUE or FALSE
#' @param folder The path specifies the location of the folder where the data is saved when the 'save' parameter is set to TRUE.
#'
#' @return
#' @export
#'
#' @examples
get_miRNA_IsoformEQ <- function(project,save = FALSE,folder = "."){
  message(paste0("=====================",project," Starting====================="))
  query <- TCGAbiolinks::GDCquery(project = project,
                    data.category = "Transcriptome Profiling",
                    data.type = "Isoform Expression Quantification")
  TCGAbiolinks::GDCdownload(query = query,files.per.chunk=10)
  mir <- TCGAbiolinks::GDCprepare(query,save = F,summarizedExperiment = TRUE)##保存的数据对象加载后名称为data
  dim(mir)
  mir <- mir[grep("mature",mir$miRNA_region),]
  mir$miRNA_region <- gsub("mature,","",mir$miRNA_region)

  acc <- unique(mir$miRNA_region)
  ano <- AnnotationDbi::select(miRBaseVersions.db,
                               keys = acc,
                               keytype = "VW-MIMAT-21.0",
                               columns = c("ACCESSION","NAME"))
  mir$mirName <- ano$NAME[match(mir$miRNA_region,ano[,"ACCESSION"])]
  message("Processing data......")
  mirCount <- dlTCGAmiRNAdata(mir[,c("mirName","read_count","barcode")])
  mirRPM <- dlTCGAmiRNAdata(mir[,c("mirName","reads_per_million_miRNA_mapped","barcode")])#
  miRNA_data <- list(mirCount = mirCount,RPM = mirRPM)
  if(save == TRUE){
    save(miRNA_data,file = paste0(folder,"/",project,"-miRNA_IsoformExpressionQuantification.Rdata"))
  }
  return(miRNA_data)
}

#' get_miRNAEQ
#'
#' @param project A list of valid project (see list with TCGAbiolinks:::getGDCprojects()$project_id)]
#' @param save TRUE or FALSE
#' @param folder The path specifies the location of the folder where the data is saved when the 'save' parameter is set to TRUE.
#'
#' @return
#' @export
#'
#' @examples
get_miRNAEQ <- function(project,save = FALSE,folder = "."){
  message(paste0("=====================",project," Starting====================="))
  query <- TCGAbiolinks::GDCquery(project = project,
                    data.category = "Transcriptome Profiling",
                    data.type = "miRNA Expression Quantification")
  TCGAbiolinks::GDCdownload(query = query,files.per.chunk=10)
  mir <- TCGAbiolinks::GDCprepare(query,save = F,summarizedExperiment = TRUE)##保存的数据对象加载后名称为data
  rownames(mir) <- mir$miRNA_ID
  mirCount <- dplyr::select(mir, matches("^read_count_"))
  colnames(mirCount) = gsub("read_count_","",colnames(mirCount))

  mirRPM <- dplyr::select(mir, matches("^reads_per_million_miRNA_mapped_"))
  colnames(mirRPM) = gsub("reads_per_million_miRNA_mapped_","",colnames(mirRPM))
  miRNA_data <- list(mirCount = mirCount,RPM = mirRPM)
  if(save == TRUE){
    save(miRNA_data,file = paste0(folder,"/",project,"-miRNAExpressionQuantification.Rdata"))
  }
  return(miRNA_data)
}


#' getMetData
#'
#' @param project A list of valid project (see list with TCGAbiolinks:::getGDCprojects()$project_id)]
#' @param save TRUE or FALSE
#' @param folder The path specifies the location of the folder where the data is saved when the 'save' parameter is set to TRUE.
#'
#' @return
#' @export
#'
#' @examples
getMetData <- function(project,save = FALSE,folder = "."){
  message(paste0("=====================",project," Starting====================="))
  query <- TCGAbiolinks::GDCquery(project = project,
                                  data.category = "DNA Methylation",
                                  data.type = "Methylation Beta Value",
                                  workflow.type = "SeSAMe Methylation Beta Estimation")
  GDCdownload(query = query,files.per.chunk=10)

  MetData <- Met450prepare(query)##readMet450data.R
  if(save == TRUE){
    save(MetData,file = paste0(folder,"/",project,"-450MethylationBetaValue.Rdata"))
  }
  return(MetData)
}





# getCNV.GeneLevel <- function(project,save = FALSE,folder = "."){
#   message(paste0("=====================",project," Starting====================="))
#   query <- TCGAbiolinks::GDCquery(project = project,
#                     data.category = "Copy Number Variation",
#                     data.type = "Gene Level Copy Number",
#                     workflow.type = "ASCAT2")
#   TCGAbiolinks::GDCdownload(query = query,files.per.chunk=10)
#   CNVdata.gl <- TCGAbiolinks::GDCprepare(query)
# }





#' In the expression matrix of TCGA tumor samples, the data of duplicate patients were deleted.
#'
#' @param data  Gene expression matrix of TCGA tumor samples.
#' @param col_rename TRUE or FALSE
#'
#' @return data.frame
#' @export
#'
#' @examples
del_dup_sample <- function(data,col_rename =T){
  data <- data[,sort(colnames(data))]
  pid = gsub("(.*?)-(.*?)-(.*?)-.*","\\1-\\2-\\3",colnames(data))
  reps = pid[duplicated(pid)]##重复样本
  if(length(reps)== 0 ){
    message("There were no duplicate patient samples for this data")
  }else{
    data <- data[,!duplicated(pid)]
  }
  if(col_rename == T){
    colnames(data) <- gsub("(.*?)-(.*?)-(.*?)-.*","\\1-\\2-\\3",colnames(data))
  }
  return(data)
}

#' RNAseq data in TCGA database were filtered according to genotype.
#'
#' @param expr Gene expression data.
#' @param fil_col A string. Column name to be filtered.
#' @param filter A value in the fil_col column.
#'
#' @return
#' @export
#'
#' @examples
filterGeneTypeExpr <- function(expr,fil_col = "gene_type",filter = FALSE){
  #filter :"protein_coding","lncRNA","miRNA","misc_RNA","snRNA","miRNA","scRNA",....

  ##Delete unexpressed genes(rows)in all samples from expr.
  dat <- expr[apply(expr[,-c(1:3)], 1, var)!=0,]
  ##rowSums
  dat <- dplyr::mutate(dat,Sums = rowSums(dat[,-c(1:3)]),.before = 4)
  dat <- dplyr::arrange(dat,gene_name,desc(Sums))
  dat <- dat[!duplicated(dat$gene_name),]
  rownames(dat) <- dat$gene_name
  if (filter == FALSE) {
    return(dat[,-c(1:4)])
  }else{
    colnames(dat)[grep(fil_col,colnames(dat))] = "fil_col"
    dat = subset(dat,fil_col == "protein_coding") %>% as.data.frame()
    rownames(dat) <- dat$gene_name
    dat <- dat[,-c(1:4)]
    # dat <- dat[dat[,fil_col] == filter,-c(1:4)]
    return(dat)
  }
}


