


geneSymbol = "FBXO22"
datafolder = "G:/DatabaseData/TCGA/new/processedTCGAdata/TCGA-STAR_Exp"
getGeneExpData.fancancer <- function(datafolder,
                                     geneSymbol,
                                     geneType = "protein_coding",
                                     dataType = "tpm",
                                     pattern = "STARdata.Rdata$",
                                     nnorn = 10){
  FilePath <- dir(datafolder,pattern,full.names = T)
  if(length(FilePath) > 0){
    ###TCGA数据库中33中癌症类型
    projects <- getGDCprojects()$project_id
    projects <- projects[grep("TCGA-",projects)]
    geneexpdata <- data.frame()
    for(project in projects){
      message(project)
      load(FilePath[grep(project,FilePath)])#STARdata

      data <- STARdata[[dataType]]
      data <- filterGeneTypeExpr(expr =data,
                                fil_col = "gene_type",
                                filter = geneType)
      if(geneSymbol %in% rownames(data)){
        turexp <- splitTCGAmatrix(data = data,sample = "Tumor")
        norexp <- splitTCGAmatrix(data = data,sample = "Normal")
        if(!is.null(ncol(norexp))){
          if(ncol(norexp) >= nnorn){
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
df = getGeneExpData.fancancer(datafolder,
                              geneSymbol,
                              geneType = "protein_coding",
                              pattern = "STARdata.Rdata$",
                              nnorn = 10)

data = df
data$SampleType = gsub("[(].*?[)]","",data$Sample)

ggplotGeneFancancerExp <- function(data,folder = "."){
  p <- ggplot2::ggplot(data, aes(x= cancer, y = exp,color = SampleType)) +
    ggplot2::geom_boxplot(aes(color = SampleType),alpha =1,
                 lwd = 0.1, outlier.size = 1,
                 outlier.colour = "white")+ #color = c("red", "blue"),
    theme_bw()+
    labs(y= expression(log[2](TPM + 1)))+
    ggpubr::stat_compare_means(label = "p.signif") +
    theme(legend.position = "top",
          axis.text.x = element_text(face = "bold",colour = "#1A1A1A"),
          axis.text.y = element_text(face = "bold",colour = "#1A1A1A"),
          axis.title = element_text(size = 12,face = "bold", colour = "#1A1A1A"),
          legend.title = element_blank(),
          axis.title.x = element_blank(),
          legend.text = element_text(size = 12, face = "bold",colour = "#1A1A1A")
    )
  ggsave(filename = paste0(folder,"/all_cancer.pdf"),plot = p,height=5,width= 8)
}
