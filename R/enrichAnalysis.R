#' exeGO_KEGG
#'
#' @param geneset For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param geneType For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param organism For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param showCategory For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param pvalueCutoff For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param OrgDb For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param KEGG For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param save For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param folder For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param Prefix For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param height For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param width For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#'
#' @return list
#' @export exeGO_KEGG
#'

exeGO_KEGG <- function(geneset,
                       geneType = "SYMBOL",
                       organism = "hsa",
                       showCategory = 6,
                       pvalueCutoff = 0.05,
                       OrgDb="org.Hs.eg.db",
                       KEGG = TRUE,
                       save = TRUE,
                       folder = "./",
                       Prefix = "enrich",
                       height = 6,
                       width = 10){
  if(geneType == "SYMBOL"){
    geneset <- clusterProfiler::bitr(geneset,
                                     fromType= "SYMBOL",
                                     toType="ENTREZID",
                                     OrgDb=OrgDb)
    geneset = unique(geneset$ENTREZID)
  }
  if(geneType == "ENTREZID"){geneset = unique(geneset)}
  ifelse(dir.exists(folder),FALSE,dir.create(folder,recursive = T))
  enrichFigBar <- list()
  for(gotype in c("BP","MF","CC")){
    message("GO enrichment analysis")
    go <- clusterProfiler::enrichGO(geneset, OrgDb, ont = gotype, pvalueCutoff=pvalueCutoff)
    GOfileName <- paste0(Prefix,"-",gotype)
    enrichFigBar[[gotype]] <- GO_KEGG.enrichVisual_barplot(go,showCategory = showCategory,
                                                           save = save,
                                                           folder = folder,
                                                           fileName = GOfileName,
                                                           height = height,
                                                           width = width)
  }
  ##======================kegg
  if(KEGG == TRUE){
    message("starting KEGG...........")
    R.utils::setOption( "clusterProfiler.download.method",'auto' )
    # R.utils::setOption( "download.file.method",'wininet' )
    kegg = clusterProfiler::enrichKEGG(gene = geneset,keyType = "kegg",
                                       organism = organism,pvalueCutoff = pvalueCutoff)
    KEGGfileName <- paste0(Prefix,"-","KEGG")
    enrichFigBar[["KEGG"]] <- GO_KEGG.enrichVisual_barplot(go,
                                                           showCategory = showCategory,
                                                           save = save,
                                                           folder = folder,
                                                           fileName = KEGGfileName,
                                                           height = height,
                                                           width = width)
  }
  return(enrichFigBar)
}


#' exeGO_KEGG2
#'
#' @param geneset For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param keyType For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param organism For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param showCategory For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param color For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param pvalueCutoff For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param OrgDb For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param pAdjustMethod For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param qvalueCutoff For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param minGSSize For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param maxGSSize For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param readable For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param KEGG For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param save For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param folder For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param Prefix For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param height For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param width For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#'
#' @return list
#' @export exeGO_KEGG2
#'
exeGO_KEGG2 <- function(geneset,
                        keyType = "SYMBOL",
                        organism = "hsa",
                        showCategory = 4,
                        color= NULL,
                        pvalueCutoff = 0.05,
                        OrgDb="org.Hs.eg.db",
                        pAdjustMethod = "BH",
                        qvalueCutoff = 0.2,
                        minGSSize = 10,
                        maxGSSize = 500,
                        readable = FALSE,
                        KEGG = TRUE,
                        save = TRUE,
                        folder = "./",
                        Prefix = "enrich",
                        height = 5,
                        width = 5.5){
  go <- clusterProfiler::enrichGO(geneset,
                                  OrgDb,
                                  keyType = keyType,
                                  ont = "ALL",
                                  pvalueCutoff=pvalueCutoff,
                                  pAdjustMethod = pAdjustMethod,
                                  qvalueCutoff = qvalueCutoff,
                                  minGSSize = minGSSize,
                                  maxGSSize = maxGSSize,
                                  readable = readable)
  go_top <- as.data.frame(go) %>%
    group_by(ONTOLOGY) %>%
    slice_head(n = 4) %>%  # 获取每个 GO 分类的前四个
    arrange(desc(p.adjust)) %>%
    ungroup() %>%  # 取消分组
    dplyr::select(ONTOLOGY, everything())
  enrich <- list( GO = go)
  ##======================kegg
  if(KEGG == TRUE){
    message("starting KEGG...........")
    R.utils::setOption( "clusterProfiler.download.method",'auto' )
    # R.utils::setOption( "download.file.method",'wininet' )
    kegg = clusterProfiler::enrichKEGG(gene = geneset,keyType = "kegg",
                                       organism = organism,pvalueCutoff = pvalueCutoff)
    enrich[["KEGG"]] <- kegg
    # 示例：KEGG 分析
    kegg_top <- as.data.frame(kegg) %>%
      dplyr::arrange(desc(p.adjust)) %>%
      dplyr::slice(1:4) %>%
      dplyr::select(ID:Count) %>%
      dplyr::mutate(ONTOLOGY = "KEGG") %>%
      dplyr::select(ONTOLOGY, everything())

    # 合并 GO 和 KEGG 数据
    enrich_result <- rbind(go_top, kegg_top) %>%
      dplyr::mutate(ONTOLOGY = factor(ONTOLOGY, levels = rev(c("BP", "CC", "MF", "KEGG")), ordered = TRUE)) %>%
      dplyr::arrange(desc(ONTOLOGY), -log10(p.adjust))
  }else{
    enrich_result <- go_top %>%
      dplyr::mutate(ONTOLOGY = factor(ONTOLOGY, levels = rev(c("BP", "CC", "MF")),
                                      ordered = TRUE)) %>%
      dplyr::arrange(desc(ONTOLOGY), -log10(p.adjust))
  }
  enrich_result$Term <- factor(stringr::str_to_title(enrich_result$Description),
                               levels = rev(stringr::str_to_title(enrich_result$Description)))
  # 清洗数据
  enrich_result$geneID <- gsub("/", ".", enrich_result$geneID)
  enrich_result$CountNumber <- enrich_result$Count / 1e3

  # 选择颜色
  if (is.null(color)) {
    color <- rev(c("#FFA39F", "#A7CB65", "#D0D8EA", "#16C2A5"))
  }
  plot <- ggplot(enrich_result) +
    geom_bar(aes(x = -log10(p.adjust), y = interaction(Term, ONTOLOGY),
                 fill = ONTOLOGY), stat = "identity") +
    scale_fill_manual(values = color, name = "ONTOLOGY") +
    geom_text(aes(x = 0.1, y = interaction(Term, ONTOLOGY),
                  label = Term), size = 3, hjust = 0, color = "black") +
    geom_text(aes(x = 0.1, y = interaction(Term, ONTOLOGY),
                  label = geneID), size = 2, hjust = 0, vjust = 2.5, color = "black") +
    geom_point(aes(x = -max(Count)/20, y = interaction(Term, ONTOLOGY),
                   size = Count, fill = ONTOLOGY), shape = 21) +
    geom_text(aes(x = -max(Count)/20, y = interaction(Term, ONTOLOGY), label = Count), size = 3) +
    scale_size(range = c(4, 8), guide = guide_legend(override.aes = list(fill = "black"))) +
    # scale_x_continuous(expand = expansion(mult = c(0, 0.2)), limits = c(-2, max(data$Count))) +
    guides(fill = guide_legend(reverse = TRUE)) +
    labs(x = expression(-log[10]("FDR")), y = "Trem name") +
    theme(
      legend.title = element_text(color = "#000080", size = 12),
      legend.text = element_text(color = "#000000", size = 10),
      axis.text.x = element_text(color = "#000000", size = 12),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_text(color = "#000080", size = 14),
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.background = element_blank()
    )
  enrich[["plot"]] <- plot
  if(save == TRUE){
    ifelse(dir.exists(folder),FALSE,dir.create(folder,recursive = T))
    ggsave(filename = paste0(folder,Prefix,".pdf"),plot = plot,height = height,width = width)
  }
  return(enrich)
}

#' exe.gseGO_GSEA
#'
#' @param data For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param gseGO.ont For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param keyType For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param OrgDb For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param species For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param TERM2GENE For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param category For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param subcategory For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param showCategory For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param pvalueCutoff For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param fileName For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param save For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param height For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param width For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param folder For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#'
#' @return list
#' @export exe.gseGO_GSEA
#'
exe.gseGO_GSEA <- function(data,
                           gseGO.ont = "BP",
                           keyType = "SYMBOL",
                           OrgDb = org.Hs.eg.db,
                           species = "Homo sapiens",
                           TERM2GENE = "msigdbr",
                           category = "C5",
                           subcategory = "GO:BP",
                           showCategory = 6,
                           pvalueCutoff=0.01,
                           fileName = "enrich",
                           save = TRUE,
                           height =4,
                           width = 7,
                           folder = "./"){
  ifelse(dir.exists(folder),"FALSE",dir.create(folder,recursive = T))
  sordeg <- arrange(data, desc(logFC))
  geneList = sordeg$logFC
  names(geneList) = sordeg$gene_name
  geneList <- geneList[!duplicated(names(geneList))]
  set.seed(42)
  message("starting gseGO...........")
  enrich <- gseGO(geneList,
                  keyType = keyType,
                  OrgDb = OrgDb,
                  ont = gseGO.ont)
  p <- enrichmentNetwork(enrich@result[1:min(100,nrow(enrich@result)),])
  if(save == TRUE){
    ifelse(dir.exists(folder),FALSE,dir.create(folder,recursive = T))
    ggsave(filename = paste0(folder,fileName,"-gseGO_enrichmentNetwork.pdf"),
           plot = p,
           height = 6,width = 6)
  }
  gseadata <- list(gseGO = list(),
                   GSEA = list())
  gseadata[["gseGO"]][[gseGO.ont]] <- enrich
  gseadata[["gseGO"]][["enrichmentNetwork"]] <- p
  if(TERM2GENE == "msigdbr"){
    gs = msigdbr::msigdbr(species = species, category = category,subcategory = subcategory)
    gs = gs[,c("gs_name","gene_symbol")]
    colnames(gs) = c("term","gene")
  }else{gs = TERM2GENE}
  egmt <- GSEA(geneList = geneList,
               TERM2GENE = gs,
               verbose = FALSE,
               pvalueCutoff = pvalueCutoff)
  gseadata[["GSEA"]][["GSEA"]] <- egmt
  if(nrow(egmt@result)>0){
    maxnum <- min(showCategory,nrow(egmt@result))
    fig <- gseaNb(egmt,
                  geneSetID = egmt$Description[1:maxnum],
                  curveCol = brewer.pal(maxnum,"Dark2"),
                  subPlot = 2)
    gseadata[["GSEA"]][["Figure"]] <- fig
    if(save == TRUE){
      ggsave(filename = paste0(folder,fileName,"-GSEA.pdf"),
             plot = fig,
             height = height,width = width)
    }
  }
  return(gseadata)
}
