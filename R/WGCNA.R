#' WGCNA.blockwiseModules
#'
#' @param expr TPM matrix.
#' @param recal TURE or FALSE.
#' @param folder Folder path
#'
#' @return list
#' @export WGCNA,ComplexHeatmap,grid
#'
WGCNA.blockwiseModules <- function(expr
                                   ,recal = FALSE
                                   ,folder="."
){
  library(WGCNA)
  ifelse(dir.exists(folder),"",dir.create(folder,recursive = T))
  filename <- paste0(folder,"/blockwiseModules_output.Rdata")
  if(!file.exists(filename) | recal== TRUE){
    WGCNA::enableWGCNAThreads() # 允许多线程运行
    wgcnadat <- log2(as.data.frame(t(expr)) + 1)
    powers = c(c(1:10), seq(from = 12, to=30, by=2))
    sft = WGCNA::pickSoftThreshold(wgcnadat, powerVector = powers, verbose = 5)
    if(!is.na(sft$powerEstimate)){
      # Plot the results:
      pdf(file = paste0(folder,"/Scale-free topology.power.pdf"),height = 5,width = 9)
      # WGCNA::sizeGrWindow(9, 5)
      par(mfrow = c(1,2))
      plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
           xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
           main = paste("Scale independence"));
      text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
           labels=powers,cex= 0.9,col="red");
      # this line corresponds to using an R^2 cut-off of h
      # 用红线标出R^2的参考值
      abline(h=0.90,col="red")
      # abline(h=0.80,col="red")
      # Mean connectivity as a function of the soft-thresholding power
      plot(sft$fitIndices[,1], sft$fitIndices[,5],
           xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
           main = paste("Mean connectivity"))
      text(sft$fitIndices[,1], sft$fitIndices[,5],
           labels=powers, cex=0.9,col="red")
      dev.off()
      ### 一步构建网络 ###

      message("Automatic network construction and module detection")
      net = WGCNA::blockwiseModules(wgcnadat, #处理好的表达矩阵
                                    power = sft$powerEstimate,#选择的软阈值 sft$powerEstimate
                                    #拓扑矩阵类型，none表示邻接矩阵聚类，unsigned最常用，构建无方向
                                    TOMType = "unsigned",
                                    minModuleSize = 30,#网络模块包含的最少基因数
                                    reassignThreshold = 0, #模块间基因重分类的阈值
                                    mergeCutHeight = 0.25,#合并相异度低于0.25的模块
                                    #true，返回模块的数字标签 false返回模块的颜色标签
                                    numericLabels = TRUE,
                                    #调用动态剪切树算法识别网络模块后，进行第二次的模块比较，合并相关性高的模块
                                    pamRespectsDendro = FALSE,
                                    saveTOMs = FALSE,#保存拓扑矩阵
                                    saveTOMFileBase = folder,
                                    verbose = 3)#0，不反回任何信息，＞0返回计算过程
      ### 样本聚类
      sampleTree = hclust(dist(wgcnadat), method = "average")

      netdata <- list(net = net,sft = sft,sampleTree =sampleTree)
      save(netdata,file = filename)
    }else{return(NULL)}
  }else(load(filename))
  if(exists("netdata")){
    return(netdata)
  }


}



#' WGCNA.ModulesPhenotype
#'
#' @param expr TPM matrix.
#' @param phenotype For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param recal For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param TCGA For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param filter For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param colors For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param heatmapH For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param heatmapW For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param moduleTraitCor For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param MM For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param GS For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param folder For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#'
#' @return finish
#' @export WGCNA,ComplexHeatmap,grid
#'
WGCNA.ModulesPhenotype <- function(expr
                                   ,phenotype
                                   ,TCGA = FALSE
                                   ,filter = "protein_coding"
                                   ,recal = FALSE
                                   ,outputNet = TRUE
                                   ,colors = blueWhiteRed(50)
                                   ,heatmapH=NULL
                                   ,heatmapW=NULL
                                   ,mtc = 0.7
                                   ,MM=0.6
                                   ,GS=0.6
                                   ,folder="."){
  if(TCGA == TRUE){
    # expr <- expr[["tpm"]]
    expr<- filterGeneTypeExpr(expr =expr,
                              fil_col = "gene_type",
                              filter = filter)
    expr <- splitTCGAmatrix(data = expr,sample = "Tumor")
    expr <- delTCGA_dup_sample(expr,col_rename = TRUE)
  }
  if(identical(colnames(phenotype),colnames(expr))){
    ifelse(dir.exists(folder),"",dir.create(folder,recursive = T))
    netdata <- WGCNA.blockwiseModules(expr,folder = folder,recal = recal)
    if(!is.null(netdata)){
      sft <- netdata[["sft"]]
      net <- netdata[["net"]]
      sampleTree <- netdata[["sampleTree"]]

      moduleLabels = net$colors
      moduleColors = labels2colors(net$colors) # 将标签转化为颜色标签
      MEs = net$MEs
      geneTree = net$dendrograms[[1]];
      #id = intersect(rownames(wgcnadat),rownames(Score))
      phenotype <- as.data.frame(t(phenotype))
      wgcnadat <- log2(as.data.frame(t(expr)) + 1)
      # identical(rownames(wgcnadat),rownames(phenotype))

      traitColors = numbers2colors(phenotype, signed = FALSE)
      # Plot the sample dendrogram and the colors underneath.


      pdf(file = paste0(folder,"/Sample dendrogram and trait heatmap.pdf"),
          width = 15,height = 9)
      plotDendroAndColors(sampleTree, traitColors,groupLabels = names(phenotype),
                          dendroLabels = F,main = "Sample dendrogram and trait heatmap")
      dev.off()
      # 定义基因和样本的数量
      nGenes = ncol(wgcnadat);
      nSamples = nrow(wgcnadat);
      # Recalculate MEs with color labels
      MEs0 = moduleEigengenes(wgcnadat, moduleColors)$eigengenes
      MEs = orderMEs(MEs0)
      # 计算模块与特征之间的相关性
      moduleTraitCor = cor(MEs, phenotype, use = "p");
      moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
      # Will display correlations and their p-values
      textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                          signif(moduleTraitPvalue, 1), ")", sep = "");
      dim(textMatrix) = dim(moduleTraitCor)
      # par(mar = c(6, 8.5, 3, 3));
      # Display the correlation values within a heatmap plot
      split = 1:length(names(MEs))#对应group分组的个数
      ha <- ComplexHeatmap::rowAnnotation(
        empty = ComplexHeatmap::anno_empty(border = FALSE),
        foo = ComplexHeatmap::anno_block(
          gp = grid::gpar(fill = substring(names(MEs), 3))
        )
      )
      if(is.null(heatmapH)){
        heatmapH <- max(length(split)*0.6,8)
      }
      if(is.null(heatmapW)){
        heatmapW <- max(ncol(phenotype),8)
      }
      # color.key<-c("#3300CC","#3399FF","white","#FF3333","#CC0000")
      pdf(file = paste0(folder,"/Module-trait relationships-Style-2.pdf"),
          width = heatmapW,height = heatmapH)
      #colorRampPalette(color.key)(50)
      ComplexHeatmap::Heatmap(
        matrix = as.matrix(moduleTraitCor),
        col= colors,# blueWhiteRed(50)
        name = " ",
        #column_split = split,
        row_split = split,
        left_annotation = ha,
        cluster_rows = F,
        column_title = NULL,
        row_title = NULL,
        column_names_rot = -45,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid::grid.text(textMatrix[i, j],
                          x, y,
                          gp = grid::gpar(fontsize = 10))}
      ) %>% print()
      dev.off()
      pdf(file = paste0(folder,"/Module-trait relationships-Style-1.pdf"),
          width = heatmapW,height = heatmapH)
      labeledHeatmap(Matrix = moduleTraitCor,
                     xLabels = names(phenotype),
                     yLabels = names(MEs),
                     ySymbols = names(MEs),
                     colorLabels = FALSE,
                     colors = colors,
                     textMatrix = textMatrix,
                     setStdMargins = FALSE,
                     cex.text = 0.7,
                     zlim = c(-1,1),
                     main = paste("Module-trait relationships"))
      dev.off()
      ##计算表达值与MEs之间的相关性
      modNames = substring(names(MEs), 3)#取出颜色名称字符串
      geneModuleMembership = as.data.frame(cor(wgcnadat, MEs, use = "p"));
      MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership),
                                                nSamples));
      names(geneModuleMembership) = paste("MM", modNames, sep="");
      names(MMPvalue) = paste("p.MM", modNames, sep="");
      dim(geneModuleMembership)
      ###计算表型与表达值之间的相关性
      geneTraitSignificance = as.data.frame(cor(wgcnadat, phenotype, use = "p"));
      GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance),
                                                nSamples));
      names(geneTraitSignificance) = paste("GS.", names(phenotype), sep="");
      names(GSPvalue) = paste("p.GS.", names(phenotype), sep="");

      # chooseOneHubInEachModule(wgcnadat,moduleColors)
      # chooseTopHubInEachModule(wgcnadat,moduleColors)
      ###绘制相关性散点图
      ###
      Trait = colnames(moduleTraitCor)
      module = rownames(moduleTraitCor)
      # mod = module[1]
      # sft = WGCNA::pickSoftThreshold(wgcnadat, powerVector = powers, verbose = 5)
      dirName <- paste0(folder,"/trait")
      ifelse(dir.exists(dirName),"",dir.create(dirName,recursive = TRUE))
      for(mod in module){
        mod_col = substring(mod, 3)
        inModule = (moduleColors== mod_col)
        moduleGenes = rownames(geneModuleMembership)[inModule]
        # Recalculate topological overlap
        # 是否输出网络
        if(outputNet == TRUE){
          TOM = TOMsimilarityFromExpr(wgcnadat, power = sft$powerEstimate);
          ## 也是提取指定模块的基因名
          # Select the corresponding Topological Overlap
          modTOM = TOM[inModule, inModule];
          dimnames(modTOM) = list(moduleGenes,moduleGenes)
          ## 模块对应的基因关系矩阵
          cyt = exportNetworkToCytoscape(
            modTOM,
            edgeFile = paste0(dirName,"/CytoscapeInput-edges-",mod, ".txt"),
            nodeFile = paste0(dirName,"/CytoscapeInput-nodes-",mod, ".txt"),
            weighted = TRUE,
            threshold = 0.02,
            nodeNames = moduleGenes,
            nodeAttr = moduleColors[inModule]
          )
        }
        #phe = Trait[5]
        for(phe in Trait){

          if(abs(moduleTraitCor[mod,phe]) > mtc){
            gs_names = colnames(geneTraitSignificance)
            df = data.frame(MM=geneModuleMembership[moduleGenes, paste0("MM",mod_col)],
                            GS=geneTraitSignificance[moduleGenes, grep(phe,gs_names)])
            rownames(df) = moduleGenes

            hubgene = moduleGenes[abs(df$MM)> MM & abs(df$GS)> GS]
            write(hubgene,file = paste0(dirName,"/",mod,"-",phe,"-hubgene.txt"))

            pdf(file = paste0(dirName,"/MM and GS",mod,"-",phe, "-relationships.pdf"),
                width = 6,height = 6.5)
            verboseScatterplot(abs(geneModuleMembership[moduleGenes, paste0("MM",mod_col)]),
                               abs(geneTraitSignificance[moduleGenes, grep(phe,gs_names)]),
                               xlab = paste("Module Membership in",substring(mod, 3), "module"),
                               ylab = paste("Gene significance for",phe),
                               main = paste("Module membership vs. gene significance\n"),
                               cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = mod_col)
            dev.off()
          }

        }
      }
    }

  }
}




#' WGCNA.ModulesScoresys
#'
#' @param expr For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param TCGA For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param geneset For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param method For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param study.type For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param heatmapH For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param heatmapW For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param moduleTraitCor For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param MM For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param GS For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param folder For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#'
#' @return For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @export WGCNA
#'
WGCNA.ModulesScoresys <- function(expr
                                  ,TCGA = FALSE
                                  ,geneset
                                  ,method = "ssgsea"
                                  ,study.type = 'bulk_RNAseq'
                                  ,heatmapH=6
                                  ,heatmapW=14
                                  ,mtc = 0.7
                                  ,MM=0.6
                                  ,GS=0.6
                                  ,folder="."){

  if(TCGA == TRUE){
    # expr <- expr[["tpm"]]
    expr<- filterGeneTypeExpr(expr =expr,
                              fil_col = "gene_type",
                              filter = "protein_coding")
    expr <- splitTCGAmatrix(data = expr,sample = "Tumor")
    expr <- delTCGA_dup_sample(expr,col_rename = TRUE)
  }
  ifelse(dir.exists(folder),"",dir.create(folder,recursive = T))
  phenotype <- scoringSys(expr
                          ,geneset = geneset
                          ,TCGA =  FALSE
                          ,method = method
                          ,study.type = study.type
                          ,save = TRUE
                          ,folder = paste0(folder,"/scoringSys"))
  WGCNA.ModulesPhenotype(expr
                         ,phenotype
                         ,recal = FALSE
                         ,colors = colors
                         ,heatmapH=heatmapH
                         ,heatmapW=heatmapW
                         ,mtc = mtc
                         ,MM=MM
                         ,GS=GS
                         ,folder=folder)
}


#' get.geneExpbasedonCorget.geneExpbasedonCor
#'
#' @param exprFor more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param gene For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param TCGA For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param filter For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param cor.method For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param cor.coef For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param cor.p.value For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param save For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param folder For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#'
#' @return For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @export WGCNA
#'
get.geneExpbasedonCor <- function(expr
                                  ,gene
                                  ,TCGA = FALSE
                                  ,filter = "protein_coding"
                                  ,cor.method = "pearson"
                                  ,cor.coef = 0.3
                                  ,cor.p.value = 0.05
                                  ,save = TRUE
                                  ,folder="."){
  if(TCGA == TRUE){
    # expr <- expr[["tpm"]]
    expr<- filterGeneTypeExpr(expr =expr,
                              fil_col = "gene_type",
                              filter = filter)
    expr <- splitTCGAmatrix(data = expr,sample = "Tumor")
    expr <- delTCGA_dup_sample(expr,col_rename = TRUE)
  }
  if(gene %in% rownames(expr)){
    ifelse(dir.exists(folder),"",dir.create(folder,recursive = T))
    # 计算相关性系数
    cordf = cor(as.data.frame(t(expr)),
                data.frame(gene = as.numeric(expr[gene,])),method = cor.method)
    # 计算p-value
    pv = WGCNA::corPvalueStudent(cordf,nSamples = ncol(expr))
    cordat = cbind(cordf,pv) %>% as.data.frame()
    colnames(cordat) = c("coefficient","p.value")
    ###==========相关性基因的筛选条件
    sigcordat = cordat[abs(cordat$coefficient) >cor.coef & cordat$p.value < cor.p.value,]
    expr = expr[intersect(rownames(sigcordat),rownames(expr)),]
    if(save == TRUE){
      write.csv(sigcordat,file = paste0(folder,"/",gene,"-",cor.method,"-coef.csv"))
      write.csv(expr,file = paste0(folder,"/",gene,"-",cor.method,"-expr.csv"))
    }
    return(expr)
  }else{stop("Error:gene")}

}


#' WGCNA.hubgeneBasedgeneCor
#'
#' @param expr For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param gene For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param TCGA For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param filter For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param cor.method For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param cor.coef For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param cor.p.value For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param phenotype For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param geneset For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param method For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param study.type For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param heatmapH For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param heatmapW For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param moduleTraitCor For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param MM For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param GS For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param save For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param folder For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#'
#' @return For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @export WGCNA
#'
WGCNA.hubgeneBasedgeneCor <- function(expr
                                      ,gene
                                      ,TCGA = FALSE
                                      ,filter = "protein_coding"
                                      ,cor.method = "pearson"
                                      ,cor.coef = 0.3
                                      ,cor.p.value = 0.05
                                      ,phenotype=NULL
                                      ,geneset=NULL
                                      ,method = "ssgsea"
                                      ,study.type = 'bulk_RNAseq'
                                      ,heatmapH=6
                                      ,heatmapW=14
                                      ,mtc = 0.7
                                      ,MM=0.6
                                      ,GS=0.6
                                      ,save= TRUE
                                      ,folder="."){
  expr <- get.geneExpbasedonCor(expr = expr
                                ,gene = gene
                                ,TCGA = TCGA
                                ,filter = filter
                                ,cor.method = cor.method
                                ,cor.coef = cor.coef
                                ,cor.p.value = cor.p.value
                                ,save = save
                                ,folder=paste0(folder,"/Cor"))
  if(is.null(phenotype) & !is.null(geneset)){
    ifelse(dir.exists(folder),"",dir.create(folder,recursive = T))
    phenotype <- scoringSys(expr
                            ,geneset = geneset
                            ,TCGA =  FALSE
                            ,method = method
                            ,study.type = study.type
                            ,save = TRUE
                            ,folder = paste0(folder,"/scoringSys"))
  }
  if(!is.null(phenotype)){
    WGCNA.ModulesPhenotype(expr
                           ,phenotype
                           ,recal = FALSE
                           ,colors = colors
                           ,heatmapH=heatmapH
                           ,heatmapW=heatmapW
                           ,mtc = mtc
                           ,MM=MM
                           ,GS=GS
                           ,folder=paste0(folder,"/WGCNA"))
  }

}

