

#' common_bar
#'
#' @param data For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param palette For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param axisTitle.x For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param axisTitle.y For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param title For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#'
#' @return list
#' @export common_bar
#'
common_bar <- function(data, palette, axisTitle.x, axisTitle.y, title) {
  ggplot(data = data, aes(x = Count, y = Term, fill = -log10(pvalue))) +
    geom_bar(stat = "identity", width = 0.8) +
    scale_fill_distiller(palette = palette, direction = 1) +
    labs(x = axisTitle.x, y = axisTitle.y, title = title) +
    theme_classic() +
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 11),
          plot.title = element_text(size = 14, hjust = 0.5),
          legend.title = element_text(size = 13),
          legend.text = element_text(size = 11))
}

#' Title
#'
#' @param enrichResult For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param showCategory For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param palette For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param axisTitle.x For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param axisTitle.y For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param title For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param save For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param folder For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param fileName For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param height For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param width For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#'
#' @return list
#' @export GO_KEGG.enrichVisual_barplot
#'
GO_KEGG.enrichVisual_barplot <- function(obj,
                                         showCategory = 6,
                                         palette = "RdPu",
                                         axisTitle.x = "Number of Gene",
                                         axisTitle.y = "Term name",
                                         title = "Enrichment barplot",
                                         save = FALSE,
                                         folder = "./",
                                         fileName = "EnrichBar",
                                         height = 6,
                                         width = 10){
  if(class(obj) == "enrichResult"){
    showCategory = min(showCategory,nrow(obj))
    if(showCategory > 0){
      enrich_result <- as.data.frame(obj)  %>%
        slice_head(n = showCategory) %>%
        arrange(desc(p.adjust)) %>%
        dplyr::select(Description, everything())
      enrich_result <- enrich_result[1:showCategory,]
      enrich_result$Term <- factor(stringr::str_to_title(enrich_result$Description),
                                   levels = rev(stringr::str_to_title(enrich_result$Description)))
      # 创建第一个图形
      p1 <- common_bar(enrich_result, palette, axisTitle.x, axisTitle.y, title)

      # 创建第二个图形，添加额外设置
      p2 <- common_bar(enrich_result, palette, axisTitle.x, axisTitle.y, title) +
        geom_bar(stat = "identity", width = 0.8, alpha = 0.7) +
        geom_text(aes(x = 0.03, label = Term), hjust = 0) +
        theme(plot.title = element_text(size = 12, hjust = 0.5),
              axis.text.y = element_blank())
      if(save == TRUE){
        ifelse(dir.exists(folder),FALSE,dir.create(folder,recursive = T))
        ggsave(filename = paste0(folder,fileName,"-style1.pdf"),plot = p1,
               height = height,width = width)
        ggsave(filename = paste0(folder,fileName,"-style2.pdf"),plot = p2,
               height = height,width = width*0.6)
      }
      return(list(style1 = p1,style2 = p2))
    }else{return(NULL)}
  }
}

