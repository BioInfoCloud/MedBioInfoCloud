
#' dlTCGAmiRNAdata
#'
#' @param data data is a data frame where the first column represents the names of miRNAs, the second column contains the corresponding values, and the third column is the barcode.
#'
#' @return data.frame
#' @export dlTCGAmiRNAdata
#'
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



#' readMet450data
#'
#' @param project A list of valid project (see list with TCGAbiolinks:::getGDCprojects()$project_id)]
#' @param dataPath "GDCdata"
#' @param json NULL
#'
#' @return data.frame
#' @export readMet450data
#'
readMet450data = function(project = NULL,dataPath = "GDCdata",json = NULL){

  message(paste0("===========================MedBioInfoCloud:Starting==========================="))

  ###从json文件获取信息
  if(is.null(project) & is.null(json))stop("Please set  parameter:project/json")
  if(!is.null(json)){
    metadata_json <- rjson::fromJSON(file=json)
    json_info <- do.call(rbind, lapply(1:length(metadata_json), function(i){
      TCGA_Barcode <- metadata_json[[i]][["associated_entities"]][[1]][["entity_submitter_id"]]
      json_File_Info <- data.frame(cases = TCGA_Barcode)
      rownames(json_File_Info) <- metadata_json[[i]][["file_name"]]
      return(json_File_Info)
    }))
  }else if(!is.null(project)){
    query <- GDCquery(
      project = project,
      data.category = "DNA Methylation",
      data.type = "Methylation Beta Value",
      workflow.type = "SeSAMe Methylation Beta Estimation"
    )
    isServeOK()
    if(missing(query)) stop("Please set query parameter")
    json_info <- query[,1][[1]]
    json_info <- json_info[,c("file_name","cases")]
    rownames(json_info) <- json_info$file_name
    source <- ifelse(query$legacy,"legacy","harmonized")
    # files <- file.path(
    #   query$results[[1]]$project, source,
    #   gsub(" ","_",query$results[[1]]$data_category),
    #   gsub(" ","_",query$results[[1]]$data_type),
    #   gsub(" ","_",query$results[[1]]$file_id),
    #   gsub(" ","_",query$results[[1]]$file_name)
    # )
    # files <- file.path("GDCdata", files)
  }

  filepath = dir(path = dataPath,
                 pattern = "methylation_array.sesame.level3betas.txt$",
                 full.names = T,
                 recursive = T)
  # wd <- filepath[1]
  if(length(filepath)!=0){
    beta <- lapply(filepath,function(wd){
      tempPath <- unlist(strsplit(wd,"/"))
      filename <- tempPath[length(unlist(strsplit(wd,"/")))]
      message(paste0("WeChat:MedBioInfoCloud:Reading:\n",filename))
      oneSampBeta <- read.table(wd,header = F,sep = "\t")
      SampBeta <- data.frame(value = oneSampBeta[,2])
      rownames(SampBeta) <- oneSampBeta[,1]
      colnames(SampBeta) <- json_info[filename,"cases"]
      return(SampBeta)
    })
    si <- unique(unlist(lapply(beta, function(x){nrow(x)})))
    if(length(si) == 1){
      data <- do.call(cbind,beta)
    }else{

      # message(paste0("部分文件函数不一样，处理需要更长的时间"))
      mergdat <- NULL
      for(f in 1:length(beta)){
        oneSampBeta <- dplyr::mutate(beta[[f]],cg = rownames(beta[[f]]),.before = 1)
        if(is.null(mergdat)){
          mergdat <- oneSampBeta
        }else{
          mergdat <- merge(mergdat,oneSampBeta,by = "cg",all.x = T)
        }
      }
      rownames(mergdat) <- mergdat$cg
      data <- mergdat[,-1]
    }
    message(paste0("=====================MedBioInfoCloud:Finish======================="))
    return(data)
  }else{message("No related files were maped")}
}

#' Met450prepare
#'
#' @param query The GDCquery function requests the return result of the methylation data.
#' @return data.frame
#' @export Met450prepare

Met450prepare <- function(query){
  isServeOK()
  if(missing(query)) stop("Please set query parameter")
  json_info <- query[,1][[1]]
  json_info <- json_info[,c("file_name","cases")]
  rownames(json_info) <- json_info$file_name
  # source <- ifelse(query$legacy,"legacy","harmonized")
  files <- file.path(
    query$results[[1]]$project,
    gsub(" ","_",query$results[[1]]$data_category),
    gsub(" ","_",query$results[[1]]$data_type),
    gsub(" ","_",query$results[[1]]$file_id),
    gsub(" ","_",query$results[[1]]$file_name)
  )
  filepath <- file.path("GDCdata", files)
  if(length(filepath)!=0){
    beta <- lapply(filepath,function(wd){
      tempPath <- unlist(strsplit(wd,"/"))
      filename <- tempPath[length(unlist(strsplit(wd,"/")))]
      message(paste0("WeChat:MedBioInfoCloud:Reading:\n",filename))
      oneSampBeta <- read.table(wd,header = F,sep = "\t")
      SampBeta <- data.frame(value = oneSampBeta[,2])
      rownames(SampBeta) <- oneSampBeta[,1]
      colnames(SampBeta) <- json_info[filename,"cases"]
      return(SampBeta)
    })
    si <- unique(unlist(lapply(beta, function(x){nrow(x)})))
    if(length(si) == 1){
      data <- do.call(cbind,beta)
    }else{

      message(paste0("Parts of the files have variable numbers of rows and take longer to process"))
      mergdat <- NULL
      for(f in 1:length(beta)){
        oneSampBeta <- dplyr::mutate(beta[[f]],cg = rownames(beta[[f]]),.before = 1)
        if(is.null(mergdat)){
          mergdat <- oneSampBeta
        }else{
          mergdat <- merge(mergdat,oneSampBeta,by = "cg",all.x = T)
        }
      }
      rownames(mergdat) <- mergdat$cg
      data <- mergdat[,-1]
    }

    message(paste0("=====================MedBioInfoCloud:Finish======================="))
    return(data)
  }else{message("No related files were maped")}

}


#' In the expression matrix of TCGA tumor samples, the data of duplicate patients were deleted.
#'
#' @param data  Gene expression matrix of TCGA tumor samples.
#' @param col_rename TRUE or FALSE
#'
#' @return data.frame
#' @export delTCGA_dup_sample
#'
delTCGA_dup_sample <- function(data,col_rename =T){
  if(ncol(data) == 1){
    colnames(data) <- gsub("(.*?)-(.*?)-(.*?)-.*","\\1-\\2-\\3",colnames(data))
  }else{
    data <- data[,sort(colnames(data))]
    pid = gsub("(.*?)-(.*?)-(.*?)-.*","\\1-\\2-\\3",colnames(data))
    reps = pid[duplicated(pid)]##重复样本
    if(length(reps)== 0 ){
      message("There were no duplicate patient samples for this data")
    }else{
      data <- data[,!duplicated(pid)]
    }
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
#' @param rownname For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param delcol For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#'
#' @return data.frame
#' @export filterGeneTypeExpr
#'
filterGeneTypeExpr <- function (expr,
                                rownname  = "gene_name",
                                fil_col = "gene_type",
                                filter = FALSE,delcol = c(1:3)) {
  expr <- as.data.frame(expr)
  dat <- expr[apply(expr[, -c(delcol)], 1, var) != 0, ]
  dat <- dplyr::mutate(dat, Sums = rowSums(dat[, -c(delcol)]))
  dat <- dplyr::arrange(dat, fil_col, desc(Sums))
  dat <- dat[!duplicated(dat[,rownname]), ]
  dat <- dat[!is.na(dat[,rownname]),]
  if(rownname %in% colnames(dat)){
    rownames(dat) <- dat[,rownname]
  }else{stop("error:rownname not found in column names or rownname is empty")}

  if (filter == FALSE) {
    return(dat[, -c(delcol,ncol(dat))])
  }else if((filter %in% unique(expr[,fil_col]) & length(filter) == 1)){
    dat = subset(dat, fil_col == fil_col) %>% as.data.frame()
    dat <- dat[, -c(delcol,ncol(dat))]
    return(dat)
  }else{stop("error:filter")}
}

#' read_gene_level_copy_number
#'
#' @param files file path.
#' @param cases For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param summarizedExperiment
#'
#' @return data.frame
#' @export read_gene_level_copy_number
#'
read_gene_level_copy_number <- function(files,cases,summarizedExperiment = FALSE){
  message("Reading Gene Level Copy Number files")
  gistic.df <- NULL
  gistic.list <- plyr::alply(files,1,.fun = function(file) {
    #message("Reading file: ", file)
    data <- read_tsv(
      file = file,
      col_names = TRUE,
      progress = TRUE,
      col_types = readr::cols()
    )
    colnames(data)[-c(1:5)] <- paste0(cases[match(file,files)],"_",  colnames(data)[-c(1:5)])

    return(data)
  })

  # Check if the data is in the same order
  stopifnot(all(unlist(gistic.list %>% map(function(y){all(y[,1:5] ==  gistic.list[[1]][,1:5])}) )))

  # need to check if it works in all cases
  # tested for HTSeq - FPKM-UQ and Counts only
  df <- bind_cols(
    gistic.list[[1]][,1:5], # gene info
    gistic.list %>%  map_dfc(.f = function(x) x[,6:8]) # copy number, min,max
  )
  if(summarizedExperiment) {
    se <- make_se_from_gene_level_copy_number(df, cases)
    return(se)
  }
  return(df)
}

#' readGISTIC
#'
#' @param files For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param cases For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#'
#' @return For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @export readGISTIC
#'
readGISTIC <- function(files, cases){
  message("Reading GISTIC file")
  gistic.df <- NULL
  gistic.list <- plyr::alply(files,1,.fun = function(file) {
    message("Reading file: ", file)
    data <- read_tsv(
      file = file,
      col_names = TRUE,
      progress = TRUE,
      col_types = readr::cols()
    )

    aliquot <- colnames(data)[-c(1:3)]
    info <- splitAPICall(
      FUN = getBarcodefromAliquot,
      step = 20,
      items = aliquot
    )

    barcode <- as.character(info$submitter_id)[match(aliquot,as.character(info$aliquot_id))]
    colnames(data)[-c(1:3)] <- barcode
    return(data)
  })
  gistic.df <- gistic.list %>%
    join_all(by =  c("Gene Symbol","Gene ID","Cytoband"), type='full')

  return(gistic.df)
}
#' read_copy_number_variation
#' @param files For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param cases For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#'
#' @importFrom purrr map2_dfr
#' @export read_copy_number_variation
read_copy_number_variation <- function(files, cases){
  message("Reading copy number variation files")

  col_types <- ifelse(any(grepl('ascat2', files)),"ccnnnnn","ccnnnd")

  purrr::map2_dfr(
    .x = files,
    .y = cases,
    .f = function(file,case) {
      data <- readr::read_tsv(file, col_names = TRUE, col_types = col_types);
      if(!missing(case)) data$Sample <- case
      data
    })
}

#' splitTCGARNAseqdata
#'
#' @param data data.frame or matrix, the columns should be the barcode of the TCGA patient sample
#' @param sample "Tumor" or "Normal"
#'
#' @return For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @export splitTCGAmatrix
#'
splitTCGAmatrix <- function(data,sample = "Tumor"){
  ##normal
  SamN <- TCGAbiolinks::TCGAquery_SampleTypes(barcode = colnames(data),
                                              typesample = c("NT","NB","NBC","NEBV","NBM"))
  ##tumor
  SamT <- setdiff(colnames(data),SamN)
  if(sample == "Tumor"){
    return(data[,SamT])
  }else if(sample == "Normal" & length(SamN)> 1 ){
    return(data[,SamN])
  }else if(length(SamN) == 1){
    dat1 = data.frame(data[,SamN])
    colnames(dat1) = SamN
    rownames(dat1) = rownames(data)
    return(dat1)
  }
  else{
    message("expr contains no normal samples or the parameter 'sample' is incorrect.")
    return(NULL)
  }
}


#' mergeSurExp
#'
#' @param expr For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param survival For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param Timeunit For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param survivalFrome For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param TCGA For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param TCGAfrome For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param feature For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param save For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param folder For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#'
#' @return a data.frame
#' @export mergeSurExp
#'
mergeSurExp <- function(expr
                        ,survival
                        ,Timeunit=1
                        ,survivalFrome = NULL
                        ,TCGA = FALSE
                        ,TCGAfrome = "MedBioInfoCloud"
                        ,feature = NULL
                        ,save = FALSE
                        ,folder = "."
){
  if(TCGA == TRUE & TCGAfrome == "MedBioInfoCloud"){
    # expr <- STARdata[["tpm"]]
    # table(expr$gene_type)
    expr <- filterGeneTypeExpr(expr = expr
                               ,fil_col = "gene_type"
    )
    expr <- splitTCGAmatrix(data = expr
                            ,sample = "Tumor"
    )
    expr <- delTCGA_dup_sample(expr,col_rename = T)
  }
  if(!is.null(feature)){
    conFeature <- intersect(rownames(expr),feature)
    if(length(conFeature)>=1){
      expr <- expr[conFeature,]
    }else{stop("Error:feature")}

  }
  if(survivalFrome == "UCSC2022"){
    survival = dplyr::arrange(survival,desc(OS.time))
    survival = survival[!duplicated(survival$X_PATIENT),]
    rownames(survival) <- survival$X_PATIENT
    survival <- survival[,c("X_PATIENT","OS","OS.time")]
    colnames(survival) <- c("submitter_id","vitalStat","surTime")
  }else if(survivalFrome == "GDCquery_clinic"){
    colnm <- intersect(c("submitter_id","vitalStat","surTime"),colnames(survival))
    survival <- survival[,colnm]
    survival <- dplyr::arrange(survival,desc(surTime))
    survival <- survival[!duplicated(survival$submitter_id),]
    rownames(survival) <- survival$submitter_id
  }
  conSample <- intersect(colnames(expr),rownames(survival))
  if(length(conSample)>0){
    expr <- expr[,conSample]
    survival <- survival[conSample,]
    mergdata <-cbind(survival,t(expr))
  }else{stop("Sample mismatch")}
  mergdata$surTime <- mergdata$surTime/Timeunit

  if(save == TRUE){
    if(!dir.exists(folder)){
      dir.create(folder,recursive = TRUE)
    }
    save(mergdata,file = paste0(folder,"/mergeSurExp.Rdata"))
    write.csv(mergdata,file = paste0(folder,"/mergeSurExp.csv"))
  }
  return(mergdata)
}


#' mergeClinExp
#'
#' @param expr For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param clinical For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param survivalFrome For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param TCGA For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param TCGAfrome For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param feature For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param save For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param folder For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#'
#' @return data.frame
#' @export mergeClinExp
#'
mergeClinExp <- function(expr
                         ,clinical
                         ,survivalFrome = NULL
                         ,TCGA = FALSE
                         ,TCGAfrome = "MedBioInfoCloud"
                         ,feature = NULL
                         ,save = FALSE
                         ,folder = "."
){
  if(TCGA == TRUE & TCGAfrome == "MedBioInfoCloud"){
    # expr <- STARdata[["tpm"]]
    # table(expr$gene_type)
    expr <- filterGeneTypeExpr(expr = expr
                               ,fil_col = "gene_type"
    )
    expr <- splitTCGAmatrix(data = expr
                            ,sample = "Tumor"
    )
    expr <- delTCGA_dup_sample(expr,col_rename = T)
  }
  if(!is.null(feature)){
    conFeature <- intersect(rownames(expr),feature)
    if(length(conFeature)>=1){
      expr <- expr[conFeature,]
    }else{stop("Error:feature")}

  }
  if(survivalFrome == "UCSC2022"){
    clinical = dplyr::arrange(clinical,desc(OS.time))
    clinical = clinical[!duplicated(clinical$X_PATIENT),]
    rownames(clinical) <- clinical$X_PATIENT
    clinical <- clinical[,c("X_PATIENT","OS","OS.time")]
    colnames(clinical) <- c("submitter_id","vitalStat","surTime")
  }
  conSample <- intersect(colnames(expr),rownames(clinical))
  if(length(conSample)>0){
    expr <- expr[,conSample]
    clinical <- clinical[conSample,]
    mergdata <-cbind(clinical,t(expr))
  }else{stop("Sample mismatch")}

  if(save == TRUE){
    if(!dir.exists(folder)){
      dir.create(folder,recursive = TRUE)
    }
    save(mergdata,file = paste0(folder,"/mergeSurExp.Rdata"))
    write.csv(mergdata,file = paste0(folder,"/mergeSurExp.csv"))
  }
  return(mergdata)
}

#' getInfiltDataOfTCGAsample
#'
#' @param expr For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param idtype For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param datatype For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param TIMER2 For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param method For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#'
#' @return For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @export getInfiltDataOfTCGAsample
#'
getInfiltDataOfTCGAsample <- function(expr,
                                      idtype = "barcode",
                                      datatype = "tumor",
                                      TIMER2,
                                      method){
  if(is.data.frame(immdata)){
    immdata = TIMER2
  }else if(file.exists(TIMER2)){
    immdata = read.csv(TIMER2,header = T,row.names = 1,check.names = F)
  }
  else{stop("The TIMER2 parameter should be the path to the immune-infiltrating data file downloaded from the TIMER2 database or a data.frame that has been read into the file.")}

  if(method %in% c("TIMER","CIBERSORT","CIBERSORT-ABS",
                   "QUANTISEQ","MCPCOUNTER","XCELL","EPIC")){
    method.immdata <- immdata[,colnames(immdata)[grep(method,colnames(immdata))]]
    colnames(method.immdata) <- gsub(paste0("_",method),"",colnames(method.immdata))

    process_data <- function(data, sort) {
      # 排序行
      if (sort == "row") {
        data <- data[sort(rownames(data)), ]
        id = gsub("(.*?)-(.*?)-(.*?)-.*","\\1-\\2-\\3",rownames(data))
        data = data[!duplicated(id),]
        rownames(data) =  gsub("(.*?)-(.*?)-(.*?)-.*","\\1-\\2-\\3",rownames(data))
      }else if(sort == "column"){
        data = data[ ,sort(colnames(data))]
        id = gsub("(.*?)-(.*?)-(.*?)-.*","\\1-\\2-\\3",colnames(data))
        data = data[,!duplicated(id)]
        colnames(data) =  gsub("(.*?)-(.*?)-(.*?)-.*","\\1-\\2-\\3",colnames(data))
      }
      return(data)
    }

    # expr
    expr.nor.id = TCGAbiolinks::TCGAquery_SampleTypes(barcode = colnames(expr),
                                                      typesample = c("NT","NB","NBC","NEBV","NBM"))
    tur.expr = expr[,setdiff(colnames(expr),expr.nor.id)]
    if(length(expr.nor.id) > 1){
      nor.expr = expr[,expr.nor.id]
    }else if(datatype == "normal"){stop("expr should contain at least 2 normal samples.")}

    # normal
    nor.method.immdata = method.immdata[TCGAbiolinks::TCGAquery_SampleTypes(barcode = rownames(method.immdata),
                                                                            typesample = c("NT","NB","NBC","NEBV","NBM")),]

    # tumor
    tur.method.immdata = method.immdata[setdiff(rownames(method.immdata),rownames(nor.method.immdata)),]

    if(idtype == "barcode" ){
      if(nchar(colnames(expr)[1])>=15){
        colnames(tur.expr) = sub("^(.{15}).*", "\\1", colnames(tur.expr))
        if(length(expr.nor.id) > 1){
          colnames(nor.expr) = sub("^(.{15}).*", "\\1", colnames(nor.expr))
        }

      }else{
        message("The parameter idtype = 'barcode' is invalid because the character length is less than 15.
                Forced use of short IDs is equivalent to idtype = 'patient'.")
      }
    }else if(idtype == "patient" | nchar(colnames(expr)[1]) < 15){
      nor.method.immdata <- process_data(nor.method.immdata, sort = "row")
      tur.method.immdata <- process_data(tur.method.immdata, sort = "row")
      if(length(expr.nor.id) > 1){
        nor.expr <- process_data(nor.expr, sort = "column")
      }

      tur.expr <- process_data(tur.expr, sort = "column")
    }else{stop("The idtype parameter should be one of 'barcode' or 'patient'.")}


    if(datatype == "tumor"){
      con.tur.id = intersect(colnames(tur.expr),rownames(tur.method.immdata))
      filtered.data <- list(expr = tur.expr[,con.tur.id],
                            immdata = as.data.frame(t(tur.method.immdata[con.tur.id,])))
    }else if(length(expr.nor.id) > 0 & datatype == "normal"){
      con.nor.id = intersect(colnames(nor.expr),rownames(nor.method.immdata))
      filtered.data <- list(expr = nor.expr[,con.nor.id],
                            immdata = as.data.frame(t(nor.method.immdata[con.nor.id,])))
    }
  }
  return(filtered.data)
}
