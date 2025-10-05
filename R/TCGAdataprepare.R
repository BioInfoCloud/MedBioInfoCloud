#' getGeneExpData.pancancer
#'
#' @param datafolder is the folder where the data is stored, downloaded using the getTCGA_RNAseqData().
#' @param geneSymbol is a vector of gene names to be analyzed.
#' @param projects is TCGA project.
#' @param geneType refers to the fil_col in the filterGeneTypeExpr() function
#' @param dataType is one of 'tpm', 'fpkm', or 'count'.
#' @param pattern is a pattern regular expression used to match the data files in the datafolder.
#' @param paired specifies whether to only retrieve data from paired samples.
#' @param nnorm indicates the minimum number of normal samples that must be included.
#'
#' @return a data.frame
#' @export getGeneExpData.pancancer
#'
getGeneExpData.pancancer <- function(datafolder,
                                     geneSymbol,
                                     projects = "All",
                                     geneType = "protein_coding",
                                     dataType = "tpm",
                                     pattern = "STARdata.Rdata$",
                                     paired = FALSE,
                                     nnorm = 3){
  FilePath <- dir(datafolder,pattern,full.names = T)
  if(length(FilePath) > 0){
    ###TCGA数据库中33中癌症类型
    dbproj = TCGAbiolinks::getGDCprojects()$project_id
    tcga.project <- dbproj[grep("TCGA-",dbproj)]

    if(length(projects) == 1 && projects == "All"){
      projects = tcga.project
    }else{
      projects = intersect(projects,tcga.project)
    }
    if(length(projects)==0){stop("project parameter error")}
    geneexpdata <- data.frame()
    for(project in projects){
      message(project)
      load(FilePath[grep(project,FilePath)])#STARdata

      data <- STARdata[[dataType]]
      data <- filterGeneTypeExpr(expr =data,
                                 fil_col = "gene_type",
                                 filter = geneType)
      geneSymbol <- intersect(geneSymbol,rownames(data))
      if(length(geneSymbol)>=1){
        # 肿瘤样本的表达数据
        turexp <- splitTCGAmatrix(data = data,sample = "Tumor")
        # 获取肿瘤样本的基因的表达数据，已经进行了log2转换
        tur_gsExp <- log2(data.frame(turexp[geneSymbol,],check.names = T) + 1) %>% t() %>%
          as.data.frame() %>% mutate(Sample = rep(paste0("Tumor(n=",ncol(turexp),")"),ncol(turexp)),
                                     .before = 1)%>%  reshape2::melt(id.vars = "Sample")
        norexp <- splitTCGAmatrix(data = data,sample = "Normal")
        #
        if(!is.null(ncol(norexp))){
          ##构建数据框
          if(ncol(norexp) >= nnorm){
            nor_gsExp <- log2(data.frame(norexp[geneSymbol,],check.names = T) + 1) %>% t() %>%
              as.data.frame() %>% mutate(Sample = rep(paste0("Normal(n=",ncol(norexp),")"),ncol(norexp)),
                                         .before = 1)%>%  reshape2::melt(id.vars = "Sample")
          }
        }else{nor_gsExp = data.frame()}
        if(!is.null(ncol(norexp)) & paired ){
          # 配对样本
          turexp <- delTCGA_dup_sample (turexp)
          norexp <- delTCGA_dup_sample (norexp)
          concol <- intersect(colnames(turexp),colnames(norexp))
          if(length(concol) >= 2){
            turexp <- turexp[,concol]
            norexp <- norexp[,concol]
            # 获取肿瘤样本的基因的表达数据，已经进行了log2转换，会覆盖前面的。
            tur_gsExp <- log2(data.frame(turexp[geneSymbol,],check.names = T) + 1) %>% t() %>%
              as.data.frame() %>% mutate(Sample = rep(paste0("Tumor(n=",ncol(turexp),")"),ncol(turexp)),
                                         .before = 1)%>%  reshape2::melt(id.vars = "Sample")
            nor_gsExp <- log2(data.frame(norexp[geneSymbol,],check.names = T) + 1) %>% t() %>%
              as.data.frame() %>% mutate(Sample = rep(paste0("Normal(n=",ncol(norexp),")"),ncol(norexp)),
                                         .before = 1)%>%  reshape2::melt(id.vars = "Sample")
          }else{
            message("The data lacks paired samples.")
          }

        }
        data = rbind(tur_gsExp,nor_gsExp)
        colnames(data) <- c("Sample","gene","exp")
        # compaired <- list(unique(data$Sample))
        data <- mutate(data,cancer = unlist(strsplit(project,"-"))[2],.before = 1)
        geneexpdata <- rbind(geneexpdata,data)
      }
    }
  }
  return(geneexpdata)
}

#' getTCGAgeneExpDat
#'
#' @param datafolder is the folder where the data is stored, downloaded using the getTCGA_RNAseqData().
#' @param geneSymbol is a vector of gene names to be analyzed.
#' @param projects is TCGA project.
#' @param geneType refers to the fil_col in the filterGeneTypeExpr() function
#' @param dataType is one of 'tpm', 'fpkm', or 'count'.
#' @param pattern is a pattern regular expression used to match the data files in the datafolder.
#'
#' @return data.frame
getTCGAgeneExpDat <- function(datafolder,
                              geneSymbol,
                              projects = "All",
                              geneType = "protein_coding",
                              dataType = "tpm",
                              pattern = "STARdata.Rdata$"
){
  FilePath <- dir(datafolder,pattern,full.names = T)
  dbproj = TCGAbiolinks::getGDCprojects()$project_id
  tcga.project <- dbproj[grep("TCGA-",dbproj)]

  if(length(projects) == 1 && projects == "All"){
    project = tcga.project
  }else{
    project = intersect(tcga.project,projects)
  }
  if(length(project)==0){stop("project parameter error")}

  expdata = data.frame()
  for(proj in project){
    message(proj)
    load(FilePath[grep(proj,FilePath)])#STARdata
    exp <- STARdata[[dataType]]
    exp <- filterGeneTypeExpr(expr = as.data.frame(exp),
                              fil_col = "gene_type",
                              filter = geneType)

    ##过滤不表达的基因
    exp <- exp[apply(exp,1,var) !=0,]
    ##判断基因是否在表达矩阵中
    gs <- intersect(geneSymbol,rownames(exp))
    if(length(gs) == 0){stop()}
    if(length(gs)< length(geneSymbol)){message("not exist:",setdiff(geneSymbol,gs))}

    ##正常组织样本ID
    SamN <- TCGAbiolinks::TCGAquery_SampleTypes(barcode = colnames(exp),
                                                typesample = c("NT","NB","NBC","NEBV","NBM"))

    # turexp <- splitTCGAmatrix(data = exp,sample = "Tumor")
    # norexp <- splitTCGAmatrix(data = exp,sample = "Normal")
    ##肿瘤组织样本ID
    SamT <- setdiff(colnames(exp),SamN)
    tur_exp <- delTCGA_dup_sample(exp[,SamT],col_rename = T)
    tur_gsExp <- log2(data.frame(tur_exp[gs,],check.names = T) + 1) %>% t() %>%
      as.data.frame() %>% mutate(Sample = rep(paste0("Tumor(n=",ncol(tur_exp),")"),ncol(tur_exp)),
                                 .before = 1)%>%  melt(id.vars = "Sample")
    if(length(SamN) !=0){
      ###去除重复样本
      if(length(SamN) == 1){
        nor_exp = data.frame(exp[,SamN])
        colnames(nor_exp) = SamN
        rownames(nor_exp) = rownames(exp)
      }else{nor_exp = exp[,SamN]}
      nor_exp <- delTCGA_dup_sample(nor_exp,col_rename = T)
      ##构建数据框
      nor_gsExp <- log2(data.frame(nor_exp[gs,],check.names = T) + 1) %>% t() %>%
        as.data.frame() %>% mutate(Sample = rep(paste0("Normal(n=",ncol(nor_exp),")"),ncol(nor_exp)),
                                   .before = 1)%>%  melt(id.vars = "Sample")
      data <- rbind(nor_gsExp,tur_gsExp)
    }else{data = tur_gsExp}
    colnames(data) <- c("Sample","gene","exp")
    data <- mutate(data,cancer = unlist(strsplit(proj,"-"))[2],.before = 1)
    expdata = rbind(expdata,data)
  }
  return(expdata)
}




#' ggplotGenePancancerExp
#'
#' @param data is obtained from the getGeneExpData.fancancer function, representing the gene expression data.
#' @param gene is a string representing a specific gene.
#' @param save TRUE or FALSE
#' @param folder The path specifies the location of the folder where the data is saved when the 'save' parameter is set to TRUE.
#' @param paired indicates whether the data consists of paired samples.
#'
#' @return ggolot object.
#' @export ggplotGenePancancerExp
#'
ggplotGenePancancerExp <- function(data,gene,save = FALSE,folder = ".",paired = FALSE){
  ifelse(dir.exists(folder),"",dir.create(folder,recursive = T))
  if(gene %in% data$gene){
    data <- data[data$gene == gene,]
    data$SampleType = gsub("[(].*?[)]","",data$Sample)
    ggplot2::ggplot(data, aes(x= cancer, y = exp,color = SampleType)) +
      ggplot2::geom_boxplot(aes(color = SampleType),alpha =1,
                            lwd = 0.1, outlier.size = 1,
                            outlier.colour = "white")+ #color = c("red", "blue"),

      # ggplot2::theme_bw()+
      ggplot2::theme_classic()+
      labs(y= expression(log[2](TPM + 1)))+
      ggpubr::stat_compare_means(label = "p.signif",paired = paired) +
      theme(legend.position = "top",
            axis.text.x = element_text(face = "bold",colour = "#1A1A1A"),
            axis.text.y = element_text(face = "bold",colour = "#1A1A1A"),
            axis.title = element_text(size = 12,face = "bold", colour = "#1A1A1A"),
            legend.title = element_blank(),
            axis.title.x = element_blank(),
            legend.text = element_text(size = 12, face = "bold",colour = "#1A1A1A")
      )
    if(save == TRUE){
      ggsave(filename = paste0(folder,"/",gene,"-all_cancer.pdf"),plot = p,height=5,width= 8)
    }
    return(p)
  }
}


#' getGeneExpData.pancancer2
#'
#' @param datafolder is the folder where the data is stored, downloaded using the getTCGA_RNAseqData().
#' @param geneSymbol is a vector of gene names to be analyzed.
#' @param projects is TCGA project.
#' @param geneType refers to the fil_col in the filterGeneTypeExpr() function
#' @param dataType is one of 'tpm', 'fpkm', or 'count'.
#' @param pattern is a pattern regular expression used to match the data files in the datafolder.
#' @return a data.frame
#' @export getGeneExpData.pancancer2
#'
getGeneExpData.pancancer2 <- function(datafolder,
                                      geneSymbol,
                                      projects = "All",
                                      geneType = "protein_coding",
                                      dataType = "tpm",
                                      pattern = "STARdata.Rdata$"){
  FilePath <- dir(datafolder,pattern,full.names = T)
  if(length(FilePath) > 0){
    ###TCGA数据库中33中癌症类型
    dbproj = TCGAbiolinks::getGDCprojects()$project_id
    tcga.project <- dbproj[grep("TCGA-",dbproj)]

    if(length(projects) == 1 && projects == "All"){
      projects = tcga.project
    }else{
      projects = intersect(projects,tcga.project)
    }
    if(length(projects)==0){stop("project parameter error")}
    geneexpdata <- data.frame()
    for(project in projects){
      message(project)
      load(FilePath[grep(project,FilePath)])#STARdata

      data <- STARdata[[dataType]]
      data <- filterGeneTypeExpr(expr =data,
                                 fil_col = "gene_type",
                                 filter = geneType)
      geneSymbol <- intersect(geneSymbol,rownames(data))
      if(length(geneSymbol)>=1){
        # 肿瘤样本的表达数据
        turexp <- splitTCGAmatrix(data = data,sample = "Tumor")
        turexp <- delTCGA_dup_sample(turexp)
        # 获取肿瘤样本的基因的表达数据，已经进行了log2转换
        tur_gsExp <- log2(data.frame(turexp[geneSymbol,],check.names = T) + 1) %>% t() %>%
          as.data.frame() %>% mutate(barcode = colnames(turexp),Sample = rep(paste0("Tumor(n=",ncol(turexp),")"),ncol(turexp)),
                                     .before = 1)%>%  tidyr::pivot_longer(cols = -c(Sample,barcode),
                                                                          names_to = "gene",
                                                                          values_to = "exp")
        norexp <- splitTCGAmatrix(data = data,sample = "Normal")
        #
        if(!is.null(ncol(norexp))){
          norexp <- delTCGA_dup_sample(norexp)
          nor_gsExp <- log2(data.frame(norexp[geneSymbol,],check.names = T) + 1) %>% t() %>%
            as.data.frame() %>% mutate(barcode = colnames(norexp),Sample = rep(paste0("Normal(n=",ncol(norexp),")"),ncol(norexp)),
                                       .before = 1)%>%  tidyr::pivot_longer(cols = -c(Sample,barcode),
                                                                            names_to = "gene",
                                                                            values_to = "exp")
        }else{nor_gsExp = data.frame()}

        data = rbind(tur_gsExp,nor_gsExp)
        concol <- intersect(colnames(turexp),colnames(norexp))
        data$paired = "unpaired"
        if(length(concol)>1){
          # all = data$barcode
          # unlist(purrr::map(concol, ~ which(all == .x)))
          data$paired[data$barcode %in% concol] = "paired"

        }
      }

      # compaired <- list(unique(data$Sample))
      data <- mutate(data,cancer = unlist(strsplit(project,"-"))[2],.before = 1)
      geneexpdata <- rbind(geneexpdata,data)
    }
    geneexpdata$paired = as.logical(geneexpdata$paired)
    return(geneexpdata)
  }
}

#
geneExpdatauseplot <- function(data,cancer,gene,paired = FALSE){
  projdat <- data[data$cancer == cancer, ]
  projdat <- data[data$gene == gene, ]
}

