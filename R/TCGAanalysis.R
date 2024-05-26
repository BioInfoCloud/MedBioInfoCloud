#' getGeneExpData.pancancer
#'
#' @param datafolder is the folder where the data is stored, downloaded using the getTCGA_RNAseqData().
#' @param geneSymbol is a vector of gene names to be analyzed.
#' @param geneType refers to the fil_col in the filterGeneTypeExpr() function
#' @param dataType is one of 'tpm', 'fpkm', or 'count'.
#' @param pattern is a pattern regular expression used to match the data files in the datafolder.
#' @param paired specifies whether to only retrieve data from paired samples.
#' @param nnorm indicates the minimum number of normal samples that must be included.
#'
#' @return a data.frame
#' @export
#'
#' @examples
getGeneExpData.pancancer <- function(datafolder,
                                     geneSymbol,
                                     geneType = "protein_coding",
                                     dataType = "tpm",
                                     pattern = "STARdata.Rdata$",
                                     paired = FALSE,
                                     nnorm = 10){
  FilePath <- dir(datafolder,pattern,full.names = T)
  if(length(FilePath) > 0){
    ###TCGA数据库中33中癌症类型
    projects <- TCGAbiolinks::getGDCprojects()$project_id
    projects <- projects[grep("TCGA-",projects)]
    geneexpdata <- data.frame()
    for(project in projects){
      message(project)
      load(FilePath[grep(project,FilePath)])#STARdata

      data <- STARdata[[dataType]]
      data <- MedBioInfoCloud::filterGeneTypeExpr(expr =data,
                                                  fil_col = "gene_type",
                                                  filter = geneType)
      geneSymbol <- intersect(geneSymbol,rownames(data))
      if(length(geneSymbol)>=1){
        turexp <- MedBioInfoCloud::splitTCGAmatrix(data = data,sample = "Tumor")
        norexp <- MedBioInfoCloud::splitTCGAmatrix(data = data,sample = "Normal")

        if(!is.null(ncol(norexp))){

          if(ncol(norexp) >= nnorm){
            if(paired){
              # 配对样本
              turexp <- MedBioInfoCloud::del_dup_sample(turexp)
              norexp <- MedBioInfoCloud::del_dup_sample(norexp)
              concol <- intersect(colnames(turexp),colnames(norexp))
              if(length(concol) >= 2){
                turexp <- turexp[,concol]
                norexp <- norexp[,concol]
              }else{
                message("The data lacks paired samples.")
              }
            }
            ###============非配对样本==
            ##构建数据框
            nor_gsExp <- log2(data.frame(norexp[geneSymbol,],check.names = T) + 1) %>% t() %>%
              as.data.frame() %>% mutate(Sample = rep(paste0("Normal(n=",ncol(norexp),")"),ncol(norexp)),
                                         .before = 1)%>%  reshape2::melt(id.vars = "Sample")
            tur_gsExp <- log2(data.frame(turexp[geneSymbol,],check.names = T) + 1) %>% t() %>%
              as.data.frame() %>% mutate(Sample = rep(paste0("Tumor(n=",ncol(turexp),")"),ncol(turexp)),
                                         .before = 1)%>%  reshape2::melt(id.vars = "Sample")
            data <- rbind(nor_gsExp,tur_gsExp)
            head(data)
            colnames(data) <- c("Sample","gene","exp")
            # compaired <- list(unique(data$Sample))
            data <- mutate(data,cancer = unlist(strsplit(project,"-"))[2],.before = 1)
            geneexpdata <- rbind(geneexpdata,data)
          }
        }
      }
    }
    return(geneexpdata)
  }

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
#' @export
#'
#' @examples
ggplotGenePancancerExp <- function(data,gene,save = FALSE,folder = ".",paired = FALSE){
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


