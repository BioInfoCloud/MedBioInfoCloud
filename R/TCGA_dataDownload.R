
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

#' getTCGA_ProteinExp
#'
#' @param project A list of valid project (see list with TCGAbiolinks:::getGDCprojects()$project_id)]
#' @param save TRUE or FALSE
#' @param folder The path specifies the location of the folder where the data is saved when the 'save' parameter is set to TRUE.
#'
#' @return data.frame
#' @export
#'
#' @examples
getTCGA_ProteinExp <- function(project,save = FALSE,folder = "."){
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

#' getTCGA_SNV_Masked_data
#'
#' @param project A list of valid project (see list with TCGAbiolinks:::getGDCprojects()$project_id)]
#' @param save TRUE or FALSE
#' @param folder The path specifies the location of the folder where the data is saved when the 'save' parameter is set to TRUE.
#'
#' @return
#' @export
#'
#' @examples
getTCGA_SNV_Masked_data <- function(project,save = FALSE,folder = "."){
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

#' getTCGA_miRNA_IsoformEQ
#'
#' @param project A list of valid project (see list with TCGAbiolinks:::getGDCprojects()$project_id)]
#' @param save TRUE or FALSE
#' @param folder The path specifies the location of the folder where the data is saved when the 'save' parameter is set to TRUE.
#'
#' @return
#' @export
#'
#' @examples
getTCGA_miRNA_IsoformEQ <- function(project,save = FALSE,folder = "."){
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
  mirCount <- MedBioInfoCloud::dlTCGAmiRNAdata(mir[,c("mirName","read_count","barcode")])
  mirRPM <- MedBioInfoCloud::dlTCGAmiRNAdata(mir[,c("mirName","reads_per_million_miRNA_mapped","barcode")])#
  miRNA_data <- list(mirCount = mirCount,RPM = mirRPM)
  if(save == TRUE){
    save(miRNA_data,file = paste0(folder,"/",project,"-miRNA_IsoformExpressionQuantification.Rdata"))
  }
  return(miRNA_data)
}

#' getTCGA_miRNAEQ
#'
#' @param project A list of valid project (see list with TCGAbiolinks:::getGDCprojects()$project_id)]
#' @param save TRUE or FALSE
#' @param folder The path specifies the location of the folder where the data is saved when the 'save' parameter is set to TRUE.
#'
#' @return
#' @export
#'
#' @examples
getTCGA_miRNAEQ <- function(project,save = FALSE,folder = "."){
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


#' getTCGA_MethylationData
#'
#' @param project A list of valid project (see list with TCGAbiolinks:::getGDCprojects()$project_id)]
#' @param save TRUE or FALSE
#' @param folder The path specifies the location of the folder where the data is saved when the 'save' parameter is set to TRUE.
#'
#' @return
#' @export
#'
#' @examples
getTCGA_MethylationData <- function(project,save = FALSE,folder = "."){
  message(paste0("=====================",project," Starting====================="))
  query <- TCGAbiolinks::GDCquery(project = project,
                                  data.category = "DNA Methylation",
                                  data.type = "Methylation Beta Value",
                                  workflow.type = "SeSAMe Methylation Beta Estimation")
  GDCdownload(query = query,files.per.chunk=10)

  MetData <- Met450prepare(query)##readMet450data.R
  # MetData <- readMet450data(project = project)
  if(save == TRUE){
    save(MetData,file = paste0(folder,"/",project,"-450MethylationBetaValue.Rdata"))
  }
  return(MetData)
}
#' getTCGA_CNV.data
#' @param project A list of valid project (see list with TCGAbiolinks:::getGDCprojects()$project_id)]
#' @param save TRUE or FALSE
#' @param folder The path specifies the location of the folder where the data is saved when the 'save' parameter is set to TRUE.
#' @param data.type "Gene Level Copy Number" or "Gene Level Copy Number Scores"
#'
#' @importFrom purrr map2_dfr
getTCGA_CNV.data <- function(project,save = FALSE,folder = ".",data.type = "Gene Level Copy Number"){
  message(paste0("=====================",project," Starting====================="))
  query <- TCGAbiolinks::GDCquery(project = project,
                    data.category = "Copy Number Variation",
                    data.type = data.type,
                    workflow.type = "ASCAT2")
  TCGAbiolinks::GDCdownload(query = query,files.per.chunk=10)

  # CNVdata.gl <- TCGAbiolinks::GDCprepare(query)
  files <- file.path(
    query$results[[1]]$project,
    gsub(" ","_",query$results[[1]]$data_category),
    gsub(" ","_",query$results[[1]]$data_type),
    gsub(" ","_",query$results[[1]]$file_id),
    gsub(" ","_",query$results[[1]]$file_name)
  )

  files <- file.path("GDCdata", files)
  if(unique(query$results[[1]]$data_type) == "Gene Level Copy Number"){
    data <- read_gene_level_copy_number(files,cases = query$results[[1]]$cases)
  }else if(query$results[[1]]$data_type == "Gene Level Copy Number Scores"){
    data <- readGISTIC(files, cases = query$results[[1]]$cases)
  }else{
    data <- read_copy_number_variation(files = files, cases = query$results[[1]]$cases)
  }
  if(save == TRUE){
    save(MetData,file = paste0(folder,"/",project,"-",data.type,".Rdata"))
  }
  return(data)
}

#' getTCGA_ClinicalData
#'
#' @param project A list of valid project (see list with TCGAbiolinks:::getGDCprojects()$project_id)]
#' @param save TRUE or FALSE
#' @param folder The path specifies the location of the folder where the data is saved when the 'save' parameter is set to TRUE.
#' @param trim TRUE or FALSE
#'
#' @return A data.frame
#' @export
#'
#' @examples
getTCGA_ClinicalData <- function(project,save = FALSE,folder = ".",trim = TRUE){
  projects <- TCGAbiolinks::getGDCprojects()$project_id
  projects <- projects[grep("TCGA-",projects)]
  proj1 <- c("TCGA-READ","TCGA-COAD","TCGA-PAAD","TCGA-ESCA","TCGA-KIRP","TCGA-HNSC",
             "TCGA-BLCA","TCGA-STAD","TCGA-CHOL","TCGA-SKCM","TCGA-LUAD","TCGA-LIHC",
             "TCGA-KIRC","TCGA-KICH","TCGA-MESO","TCGA-LUSC","TCGA-GBM","TCGA-UVM",
             "TCGA-BRCA","TCGA-TGCT","TCGA-THCA")
  proj2 <- setdiff(projects,proj1)
  allInfo_clindata <- TCGAbiolinks::GDCquery_clinic(project, type = "clinical", save.csv = F)
  if(trim == TRUE){
    if(project %in% sub_proj){
      processed_clindata <- allInfo_clindata[,c("submitter_id","gender",
                                                "age_at_index","vital_status",
                                                "days_to_last_follow_up","days_to_death",
                                                "ajcc_pathologic_t" ,
                                                "ajcc_pathologic_n","ajcc_pathologic_m","ajcc_pathologic_stage")]
      colnames(processed_clindata) <- c("submitter_id","gender",
                                        "age","vital_status",
                                        "follow","death",
                                        "ajcc_t" ,
                                        "ajcc_n","ajcc_m","ajcc_stage")
      processed_clindata <- processed_clindata[!is.na(processed_clindata$vital_status),]
      processed_clindata$vitalStat <- ifelse(processed_clindata$vital_status=="Alive",0,1)
      processed_clindata$surTime <- ifelse(processed_clindata$vital_status=="Alive",
                                           processed_clindata$follow,
                                           processed_clindata$death)
      processed_clindata$ajcc_t <- substr(processed_clindata$ajcc_t,1,2)
      processed_clindata$ajcc_n <- substr(processed_clindata$ajcc_n,1,2)
      processed_clindata$ajcc_m <- substr(processed_clindata$ajcc_m,1,2)
      processed_clindata$ajcc_stage <- gsub("[A-F]","",processed_clindata$ajcc_stage)
      processed_clindata$ageGroup <- ifelse(processed_clindata$age >50,">50","≦50")
      ###============下面2行代码很重要
      processed_clindata <- dplyr::arrange(processed_clindata,submitter_id,desc(surTime))
      processed_clindata <- processed_clindata[!duplicated(processed_clindata$submitter_id),]
      if(save == TRUE){
        save(processed_clindata,file = paste0(folder,"/",project,"-clindata.Rdata"))
      }
      return(processed_clindata)
    }else if(project %in% proj2){
      colnm <-intersect(colnames(allInfo_clindata), c("submitter_id","gender",
                                                      "age_at_index","vital_status",
                                                      "days_to_last_follow_up","days_to_death",
                                                      "ajcc_pathologic_t" ,
                                                      "ajcc_pathologic_n","ajcc_pathologic_m","figo_stage"))
      processed_clindata <- allInfo_clindata[,colnm]

      c1 <- c("submitter_id","days_to_last_follow_up","ajcc_pathologic_t" ,"ajcc_pathologic_n",
              "gender","vital_status","age_at_index","days_to_death")
      c2 <- c("submitter_id","figo_stage","days_to_last_follow_up","gender","vital_status","age_at_index","days_to_death")
      c3 <- c("submitter_id", "days_to_last_follow_up","gender","vital_status","age_at_index","days_to_death")
      if(identical(colnm,c1)){
        colnames(processed_clindata) <- c("submitter_id","follow","ajcc_t" ,"ajcc_n","gender","vital_status","age","death")
      }else if(identical(colnm,c2)){
        colnames(processed_clindata) <- c("submitter_id","figo_stage","follow","gender","vital_status","age","death")
      }else if(identical(colnm,c3)){
        colnames(processed_clindata) <- c("submitter_id", "follow","gender","vital_status","age","death")
      }

      processed_clindata <- processed_clindata[!is.na(processed_clindata$vital_status),]
      processed_clindata$vitalStat <- ifelse(processed_clindata$vital_status=="Alive",0,1)
      processed_clindata$surTime <- ifelse(processed_clindata$vital_status=="Alive",
                                           processed_clindata$follow,
                                           processed_clindata$death)
      processed_clindata$ajcc_t <- substr(processed_clindata$ajcc_t,1,2)
      processed_clindata$ajcc_n <- substr(processed_clindata$ajcc_n,1,2)
      processed_clindata$ajcc_m <- substr(processed_clindata$ajcc_m,1,2)

      processed_clindata$figo_stage2 <- gsub("[A-F]","",processed_clindata$figo_stage)
      processed_clindata$figo_stage2 <- gsub("[0-9]","",processed_clindata$figo_stage2)
      processed_clindata$ageGroup <- ifelse(processed_clindata$age >50,">50","≦50")

      unique(processed_clindata$submitter_id)
      ###============下面2行代码很重要
      processed_clindata <- dplyr::arrange(processed_clindata,submitter_id,desc(surTime))
      processed_clindata <- processed_clindata[!duplicated(processed_clindata$submitter_id),]
      if(save == TRUE){
        save(processed_clindata,file = paste0(folder,"/",project,"-clindata.Rdata"))
      }
      return(processed_clindata)
    }
  }else{
    if(save == TRUE){
      save(allInfo_clindata,file = paste0(folder,"/",project,"-allInfo_clindata.Rdata"))
    }
    return(allInfo_clindata)
  }
}
