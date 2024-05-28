
#' outputGmtFile
#'
#' @param input list or data.frame
#' @param description NA
#' @param folder For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param filename filepath
#'
#' @return
#' @export
#'
#' @examples
outputGmtFile <- function(input,description = NA,folder= ".",filename){
  ifelse(dir.exists(folder),"",dir.create(folder,recursive = T))
  if(is.data.frame(input)){
    nms <- unique(input[,1])
    ls <- lapply(nms, function(nm){
      input[input[,1]== nms[nm],2]
    })
    names(ls) <- nms
    input = ls
  }
  if(is.list(input)){
    output <- file(paste0(folder,"/",filename), open="wt")
    nms <- names(input)
    lapply(1:length(nms),function(x){
      if(!is.na(description)){
        dsc = description[x]
      }else{dsc = description}
      outlines = paste0(c(nms[x], dsc, input[[x]]),collapse='\t')
      writeLines(outlines, con=output)
    })
    close(output)
  }
}



#' tidy.gmt
#'
#' @param filepath For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param fun "stat" or "merge"
#' @param Source For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param termName For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param save For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param folder For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param filename For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#'
#' @return data.frame
#' @export clusterProfiler
#'
#' @examples
tidy.gmt <- function(filepath
                     ,fun = "stat"
                     ,Source = ""
                     ,termName=NULL
                     ,save = TRUE
                     ,folder = "."
                     ,filename = "geneset"
                     ){
  # fl <- dir(filepath,pattern = pattern)
  ifelse(dir.exists(folder),"",dir.create(folder,recursive = T))
  if(dir.exists(filepath)){
    fp <- dir(filepath,pattern = "gmt$",full.names = T)
  }else if(is.vector(filepath) & sum(grep("gmt$",filepath)) == length(filepath)){
    fp <- filepath
  }
  stat.gs <- data.frame()
  mergedf <- data.frame()
  gs <- c()
  for(x in fp){
    ldf <- clusterProfiler::read.gmt(x)
    mergedf <- rbind(mergedf,ldf)
    gs <- unique(c(gs,ldf$gene))
    df <- data.frame(term = gsub("_"," ",as.character(unique(ldf$term)))
                     ,Source = Source
                     ,`Gene count` = nrow(ldf)
    )
    stat.gs <- rbind(stat.gs,df)
  }
  totallgene <- data.frame(term = ""
                           ,Source = Source
                           ,`Gene count` = paste0(length(unique(gs)),"(unique)")
                           )
  stat.gs <- rbind(stat.gs,totallgene)
  if (!is.null(termName)) {
    colnames(stat.gs)[1] <- termName
  }
  if(fun == "stat"){
    if(save == TRUE){
      writeLines(unique(gs),con = paste0(folder,"/",filename,"-stat.unique.txt"))
      write.csv(stat.gs,file = paste0(folder,"/",filename,"-stat.csv"))
    }
    return(stat.gs)
  } else if(fun == "merge"){
    if(save == TRUE){
      writeLines(unique(gs),con = paste0(folder,"/",filename,"-merge.unique.txt"))
      outputGmtFile(input = mergedf
                    ,description = NA
                    ,folder= folder
                    ,filename = paste0(filename,"-merge.gmt"))
    }
    return(mergedf)
  }
}


#' tidyGene.fromeReactome
#'
#' @param filepath For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param fun For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param Source For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param termName For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param save For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param folder For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param filename For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#'
#' @return data.frame
#' @export stringr
#'
#' @examples
tidyGene.fromeReactome <- function(filepath
                                   ,fun = "stat"
                                   ,Source = "Reactome"
                                   ,termName=NULL
                                   ,save = TRUE
                                   ,folder = "."
                                   ,filename = "geneset"){
  ifelse(dir.exists(folder),"",dir.create(folder,recursive = T))
  if(dir.exists(filepath)){
    fp <- dir(filepath,pattern = "tsv$",full.names = T)
  }else if(is.vector(filepath) & sum(grep("tsv$",filepath)) == length(filepath)){
    fp <- filepath
  }
  stat.gs <- data.frame()
  mergedf <- data.frame()
  gs <- c()
  for(x in fp){
    if(Source == "Reactome"){
      rc <- read.csv(x,header = T,sep = "\t")[,c("MoleculeType","MoleculeName")]
      rc <- tidyr::separate(rc,MoleculeName,c("UniProt","gene")," ")
      rc <- rc[rc$MoleculeType == "Proteins",]
      rc <- rc[!duplicated(rc$gene),]
      filenm <- unlist(strsplit(x,"/"))
      filenm <- filenm[length(filenm)]
      # 提取[号前的内容
      # str_extract(filenm, "^[^\\[]+")
      term = stringr::str_extract(filenm,"^[^\\[]*")
      PathwayID = stringr::str_extract(filenm, "(?<=\\[)[^\\]]+")
    }else{
      message("开发中")
    }
    gs <- c(gs,rc$gene)
    statdf <- data.frame(term = term
                         ,Source = Source
                         ,PathwayID = PathwayID
                         ,`Gene count` = nrow(rc))
    stat.gs <- rbind(stat.gs,statdf)
    gmtdf <- data.frame(term = rep(term,nrow(rc))
                        ,gene = rc$gene)
    mergedf <- rbind(mergedf,gmtdf)
  }
  totallgene <- data.frame(term = ""
                           ,Source = ""
                           ,PathwayID =""
                           ,`Gene count` = paste0(length(unique(gs)),"(unique)")
  )
  stat.gs <- rbind(stat.gs,totallgene)
  if (!is.null(termName)) {
    colnames(stat.gs)[1] <- termName
  }

  if(fun == "stat"){
    if(save == TRUE){
      writeLines(unique(gs),con = paste0(folder,"/",filename,"-stat.unique.txt"))
      write.csv(stat.gs,file = paste0(folder,"/",filename,"-stat.csv"))
    }
    return(stat.gs)
  } else if(fun == "merge"){
    if(save == TRUE){
      outputGmtFile(input = mergedf
                    ,description = NA
                    ,folder= folder
                    ,filename = paste0(filename,"-merge.gmt"))
      writeLines(unique(gs),con = paste0(folder,"/",filename,"-merge.unique.txt"))
    }
    return(mergedf)
  }

}

