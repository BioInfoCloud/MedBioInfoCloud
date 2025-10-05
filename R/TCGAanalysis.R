#' TCGA_Clin_Analysis
#'
#' @param data For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param clin_feature For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param feature For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param title For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param width For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param height For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param save For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param folder For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#'
#' @return For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @export TCGA_Clin_Analysis
#'
TCGA_Clin_Analysis <- function(data,clin_feature,feature,title="",
                               width = 5,height = 5,save = TRUE,folder = "."){
  if((clin_feature %in% colnames(data)) & (feature %in% colnames(data))){
    filterdata <- data[!is.na(data[,clin_feature]),c(clin_feature,feature)]
    if(nrow(filterdata) >= 6){
      colnames(filterdata) <- c("clin_feature","feature")
      # 使用combn函数生成所有可能的两个元素的组合
      # m参数设置为2，表示每次组合两个元素
      # simplify=FALSE使得结果是一个列表，而不是矩阵
      compaired <- combn(unique(filterdata[,"clin_feature"]), m = 2, simplify = FALSE)


      p1 = ggplot(filterdata,aes(clin_feature,feature,fill= clin_feature))+
        geom_boxplot(aes(fill = clin_feature),notch = FALSE,
                     position = position_dodge(width =0.01,preserve = "single"),
                     outlier.alpha  =1,width=0.4) +
        scale_fill_manual(values=c(brewer.pal(length(unique(filterdata[,"clin_feature"])),"Set1")))+
        geom_signif(comparisons = compaired,
                    step_increase = 0.1,
                    map_signif_level = T,
                    margin_top = 0.2,
                    test = "wilcox.test")+
        labs(y= paste0("The expression of ",feature,"\nlog2(TPM + 1)"),title= title)+
        theme_classic()+
        theme(panel.background=element_rect(fill="white",colour="black",size=0.25),
              plot.title = element_text(hjust = 0.5),
              axis.line=element_line(colour="black",size=0.25),
              axis.title=element_text(size=10,face="plain",color="black"),
              axis.text.x = element_text(angle = 45,face = "plain",colour = "black",hjust=1,vjust=1),
              axis.title.x=element_blank(),
              axis.title.y=element_text(size=12,face="plain",color="black"),
              axis.text = element_text(size=12,color="black"),
              legend.position="none"
        )
      p2 = ggplot(filterdata,aes(clin_feature,feature,fill= clin_feature))+
        geom_violin(aes(fill = clin_feature),trim = FALSE)+
        geom_signif(comparisons = compaired,
                    step_increase = 0.1,
                    map_signif_level = T,
                    margin_top=0.2,
                    tip_length =0.02,
                    test = "wilcox.test")+
        geom_boxplot(width = 0.1,fill = "white")+
        scale_fill_manual(values=c(brewer.pal(length(unique(filterdata[,"clin_feature"])),"Dark2")))+
        theme_classic()+
        labs(y= paste0("The expression of ",feature,"\nlog2(TPM + 1)"),title= title)+
        theme(panel.background=element_rect(fill="white",colour="black",size=0.25),
              plot.title = element_text(hjust = 0.5),
              axis.line=element_line(colour="black",size=0.25),
              axis.title.x = element_blank(),
              axis.text.x = element_text(angle = 45,face = "plain",colour = "black",hjust=1,vjust=1),
              axis.text = element_text(size=12,face="plain",color="black"),
              legend.position="none"
        )
      if(save == TRUE){
        if(!dir.exists(folder)){
          dir.create(folder,recursive = TRUE)
        }
        ggsave(filename = paste0(folder,"/",feature,"-",clin_feature,"-boxplot.pdf"),
               plot = p1,width = width,height = height)
        ggsave(filename = paste0(folder,"/",feature,"-",clin_feature,"-violin.pdf"),
               plot = p2,width = width,height = height)

      }
      if(clin_feature == "ajcc_stage"){
        filterdata$stage = ""
        filterdata$stage[which(filterdata$clin_feature %in% c("Stage I", "Stage II"))] = "Stage I & II"
        filterdata$stage[which(filterdata$clin_feature %in% c("Stage III", "Stage IV"))] = "Stage III & IV"

        if(length(unique(filterdata[,"stage"])) >=2){
          compaired <- combn(unique(filterdata[,"stage"]), m = 2, simplify = FALSE)
        }else{stop()}

        p3 = ggplot(filterdata,aes(stage,feature,fill= stage))+
          geom_boxplot(aes(fill = stage),notch = FALSE,
                       position = position_dodge(width =0.01,preserve = "single"),
                       outlier.alpha  =1,width=0.4) +
          scale_fill_manual(values=c(brewer.pal(length(unique(filterdata[,"stage"])),"Set1")))+
          geom_signif(comparisons = compaired,
                      step_increase = 0.1,
                      map_signif_level = T,
                      margin_top = 0.2,
                      test = "wilcox.test")+
          labs(y= paste0("The expression of ",feature,"\nlog2(TPM + 1)"),title= title)+
          theme_classic()+
          theme(panel.background=element_rect(fill="white",colour="black",size=0.25),
                plot.title = element_text(hjust = 0.5),
                axis.line=element_line(colour="black",size=0.25),
                axis.title=element_text(size=10,face="plain",color="black"),
                axis.text.x = element_text(angle = 45,face = "plain",colour = "black",hjust=1,vjust=1),
                axis.title.x=element_blank(),
                axis.title.y=element_text(size=12,face="plain",color="black"),
                axis.text = element_text(size=12,color="black"),
                legend.position="none"
          )
        p4 = ggplot(filterdata,aes(stage,feature,fill= stage))+
          geom_violin(aes(fill = stage),trim = FALSE)+
          geom_signif(comparisons = compaired,
                      step_increase = 0.1,
                      map_signif_level = T,
                      margin_top=0.2,
                      tip_length =0.02,
                      test = "wilcox.test")+
          geom_boxplot(width = 0.1,fill = "white")+
          scale_fill_manual(values=c(brewer.pal(length(unique(filterdata[,"stage"])),"Dark2")))+
          theme_classic()+
          labs(y= paste0("The expression of ",feature,"\nlog2(TPM + 1)"),title= title)+
          theme(panel.background=element_rect(fill="white",colour="black",size=0.25),
                plot.title = element_text(hjust = 0.5),
                axis.line=element_line(colour="black",size=0.25),
                axis.title.x = element_blank(),
                axis.text.x = element_text(angle = 45,face = "plain",colour = "black",hjust=1,vjust=1),
                axis.text = element_text(size=12,face="plain",color="black"),
                legend.position="none")

        if(save == TRUE){
          ggsave(filename = paste0(folder,"/",feature,"-",clin_feature,"-bind-boxplot.pdf"),
                 plot = p3,width = width*0.6,height = height)
          ggsave(filename = paste0(folder,"/",feature,"-",clin_feature,"-bind-violin.pdf"),
                 plot = p4,width = width*0.6,height = height)

        }
      }

    }
  }

}
