####--------------------------绘制火山图
#' plotDEGvolcanoFig
#'
#' @param data For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param x For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param y For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param cut_pvalue For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param cutFC For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param title For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param group For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param colour For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param label For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#'
#' @return For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @export plotDEGvolcanoFig
#'
plotDEGvolcanoFig <- function(data,x,y,cut_pvalue,cutFC,title,group,label,colour = c("#DC143C","#00008B", "#808080")){
  data <- data[,c(x,y,group,label)]
  colnames(data) <- c("log2FC","FDR","group","label")
  p <-  ggplot(data = data, aes(x = log2FC, y = -log10(FDR), colour = group)) + #数据映射
    geom_point(alpha = 1,size=2) +#散点图，alpha就是点的透明度
    scale_color_manual(values = colour) + #手动调颜色c("#DC143C","#00008B", "#808080")
    theme_bw() +#设定主题
    geom_text_repel(label = data$label,#便签就是基因名 geom_text_repel：adds text directly to the plot.
                    size = 4,
                    segment.color = "black", #连接线的颜色，就是名字和点之间的线
                    show.legend = FALSE) +
    theme(axis.title=element_text(size=15,face="plain",color="black"),
          axis.text = element_text(size=12,face="plain",color="black"),
          legend.title = element_blank(),
          legend.background = element_blank(),
          panel.background = element_rect(fill = "transparent",colour = "black"),
          plot.background = element_blank(),
          plot.title = element_text(size=15, lineheight=.8,hjust=0.5, face="plain"),
          legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))+
    labs(title = title)+ #"Tumor vs. Normal"
    ylab(ifelse(y == "PValue",expression(-log[10]("P Value")),expression(-log[10]("FDR")))) +#expression的作用就是让log10的10下标
    xlab(expression(log[2]("Fold Change"))) +
    geom_vline(xintercept = c(-cutFC, cutFC), #加垂直线
               lty = 2,
               col = "black",
               lwd = 0.3) +
    geom_hline(yintercept = -log10(cut_pvalue),
               lty = 2,
               col = "black",
               lwd = 0.3)
  return(p)
}
