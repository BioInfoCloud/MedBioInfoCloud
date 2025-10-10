
#' getGeneLenFromeGTF
#'
#' @param gtf The gtf parameter represents the path to a GTF (Gene Transfer Format) format file.
#'
#' @return data.frame
#' @export getGeneLenFromeGTF
#'
getGeneLenFromeGTF <- function(gtf){
  # reference 1: https://mp.weixin.qq.com/s/dwbpJ0nhzyIp9fDv7fEWEQ
  # reference 2: https://mp.weixin.qq.com/s/lazavD3jzRVO4QkxysHQcg
  # 加载必要的包
  # library(GenomicFeatures)
  # library(GenomicRanges)
  txdb <- GenomicFeatures::makeTxDbFromGFF(gtf,format="gtf")
  exons.list.per.gene <- GenomicFeatures::exonsBy(txdb, by = "gene")
  exonic.gene.sizes <- lapply(exons.list.per.gene,
                              function(x){sum(GenomicRanges::width(GenomicRanges::reduce(x)))})
  eff_length <- do.call(rbind,lapply(exonic.gene.sizes, data.frame))
  eff_length <- data.frame(Ensembl = do.call(rbind,strsplit(rownames(eff_length),'\\.'))[,1],
                           effLen = eff_length[,1])
  return(eff_length)
}


#' getGeneTypeInfoFromeGTF
#'
#' @param gtf The gtf parameter represents the path to a GTF (Gene Transfer Format) format file.
#'
#' @return data.frame
#' @export getGeneTypeInfoFromeGTF
#'
getGeneTypeInfoFromeGTF <- function(gtf){
  if (is.character(gtf)) {
    if(!file.exists(gtf)) stop("Bad gtf file.")
    message("Treat gtf as file")
    gtf = data.table::fread(gtf, header = FALSE)
  } else {
    data.table::setDT(gtf)
  }
  gtf = gtf[gtf[[3]] == "gene", ]
  pattern_id = ".*gene_id \"([^;]+)\";.*"
  pattern_name = ".*gene_name \"([^;]+)\";.*"
  pattern_type = ".*gene_type \"([^;]+)\";.*"
  gene_id = sub(pattern_id, "\\1", gtf[[9]])
  gene_id = do.call(rbind,strsplit(gene_id,'\\.'))[,1]
  gene_name = sub(pattern_name, "\\1", gtf[[9]])
  gene_type = sub(pattern_type, "\\1", gtf[[9]])

  EnsemblTOGenename <- data.frame(Ensembl = gene_id,
                                  Symbol = gene_name,
                                  gene_type = gene_type,
                                  stringsAsFactors = FALSE)
  return(EnsemblTOGenename)
}

#' getGeneBaseInfo
#'
#' @param gtf The gtf parameter represents the path to a GTF (Gene Transfer Format) format file.
#'
#' @return data.frame
#' @export getGeneBaseInfo
#'
getGeneBaseInfo <- function(gtf){
  ens2symInfo <- getGeneInfoFromeGTF(gtf)
  eff_length <- getGeneLenFromeGTF(gtf)
  GeneInfo <- merge(ens2symInfo,eff_length,by = "Ensembl")
  return(GeneInfo)
}

#' RNAseqDataConversion
#'
#' @param data RNAseq expression data, either a matrix or a data frame
#' @param type must be one of "Counts2TPM", "Counts2FPKM", or "FPKM2TPM"
#' @param species must be one of "homo", "mus", or NULL.
#' @param gtf The gtf parameter represents the path to a GTF (Gene Transfer Format) format file.
#'
#' @return either a matrix or a data frame
#' @export RNAseqDataConversion
#'
RNAseqDataConversion <- function(data,type,species = "homo",gtf = NULL){
  if(type %in% c("Counts2TPM","Counts2FPKM","FPKM2TPM")){
    if(type == "FPKM2TPM"){
      return(FPKM2TPM(data))
    }else if(species == "homo"){
      ano <- MedBioInfoCloud::hsaGeneInfo
    }else if(species == "mus"){
      ano <- MedBioInfoCloud::musGeneInfo
    }else if(is.null(species) & !is.null(gtf)){
      if(is.character(gtf) & (length(gtf)==1) & grepl(".gtf$",gtf) & file.exists(gtf)){
        ano <- getGeneBaseInfo(gtf)
      }
    }
    ano <- dplyr::arrange(ano,Symbol,desc(effLen))
    ano <- ano[!duplicated(ano$Symbol),]
    rownames(ano) <- ano$Symbol
    congene <- intersect(rownames(data),rownames(ano))
    if(!is.numeric(congene)){
      ano <- ano[congene,]
      data <- data[congene,]
      if(type == "Counts2TPM"){
        return(Counts2TPM(counts = data,effLen = ano$effLen))
      }else if(type == "Counts2FPKM"){
        return(Counts2FPKM(counts = data,effLen = ano$effLen))
      }
    }else{message("Gene annotation information does not match the row names of the data.")}

  }
}




#' Counts2TPM
#'
#' @param counts a data.frame or matrix for raw count of RNAseq.
#' @param effLen The length of genes.
#'
#' @return a data.frame
#' @export Counts2TPM
#'
Counts2TPM <- function(counts, effLen){
  rate <- log(counts) - log(effLen)
  denom <- log(sum(exp(rate)))
  return(exp(rate - denom + log(1e6)))
}
#' Counts2FPKM
#'
#' @param counts a data.frame or matrix for raw count of RNAseq.
#' @param effLen The length of genes.
#'
#' @return a data.frame
#' @export Counts2FPKM
#'
Counts2FPKM <- function(counts, effLen){
  N <- sum(counts)
  exp( log(counts) + log(1e9) - log(effLen) - log(N) )
}
#' FPKM2TPM
#'
#' @param fpkm a data.frame or matrix for fpkm of RNAseq.
#'
#' @return a data.frame
#' @export FPKM2TPM
#'
FPKM2TPM <- function(fpkm){
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}




#' Generate a GMT file from gene set data
#'
#' This function converts gene set data (data frame or named list) into a GMT (Gene Matrix Transposed) file,
#' with validation and cleaning steps to ensure data integrity. It removes empty values and invalid entries
#' to produce a standard-compliant GMT file.
#'
#' @param input A data frame or named list containing gene set information:
#'   - If data frame: Must have at least two columns. First column = gene set names (terms),
#'     second column = genes in each term. Additional columns are ignored.
#'   - If named list: Each element is a character vector of genes; list names = gene set terms.
#'
#' @param description Character vector of descriptions for gene sets. If shorter than the number of gene sets,
#'   descriptions will be recycled. If `NA` (default), empty strings will be used.
#'
#' @param folder Path to directory where the GMT file will be saved. Creates directory recursively if it doesn't exist.
#'   Defaults to current working directory (".").
#'
#' @param filename Name of the output GMT file (e.g., "pathways.gmt"). Recommended to use ".gmt" extension.
#'
#' @return No return value; writes a GMT file to the specified path.
#'
#' @details The GMT format requires each line to follow: `term<tab>description<tab>gene1<tab>gene2<tab>...`
#'   This function ensures:
#'   - Duplicate genes within a set are removed
#'   - Empty gene sets (after cleaning) are skipped with warnings
#'   - NA values and empty strings in gene lists are filtered out
#'   - Proper handling of descriptions (recycling and empty value replacement)
#'
#' @examples
#' # Example 1: Data frame input
#' df <- data.frame(
#'   term = c("GO:0001", "GO:0001", "GO:0002"),
#'   gene = c("GeneA", NA, "GeneB"),  # Contains NA (will be filtered)
#'   stringsAsFactors = FALSE
#' )
#' outputGmtFile(
#'   input = df,
#'   description = c("Cell cycle", "Signal transduction"),
#'   folder = "output",
#'   filename = "test_genesets.gmt"
#' )
#'
#' # Example 2: Named list input
#' gene_sets <- list(
#'   "GO:0001" = c("GeneA", "GeneB", ""),  # Contains empty string (will be filtered)
#'   "GO:0002" = c("GeneC", NA)            # Contains NA (will be filtered)
#' )
#' outputGmtFile(
#'   input = gene_sets,
#'   filename = "test_genesets.gmt"
#' )
#'
#' @export
outputGmtFile <- function(input, description = NA, folder = ".", filename) {
  # 1. Check for required filename parameter
  if (missing(filename)) {
    stop("The 'filename' parameter must be provided")
  }
  # Warn about non-standard file extensions
  if (!grepl("\\.gmt$", filename, ignore.case = TRUE)) {
    warning("Filename is recommended to end with '.gmt' for compatibility with downstream tools")
  }

  # 2. Create output directory if it doesn't exist
  if (!dir.exists(folder)) {
    dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  }
  file_path <- file.path(folder, filename)

  # 3. Process data frame input: convert to named list of gene sets
  if (is.data.frame(input)) {
    if (ncol(input) < 2) {
      stop("Data frame input must contain at least two columns (terms and genes)")
    }
    # Rename first two columns to standardize access
    colnames(input)[1:2] <- c("term", "gene")

    # Extract unique terms BEFORE converting to list (critical for naming)
    terms <- unique(input$term)

    # Group genes by term with cleaning steps
    input <- lapply(terms, function(term) {
      # Extract genes for current term
      genes <- input$gene[input$term == term]
      # Remove NA values
      genes <- na.omit(genes)
      # Convert to character and remove empty strings
      genes <- as.character(genes)
      genes <- genes[genes != ""]
      # Remove duplicate genes
      unique(genes)
    })
    # Name list elements using pre-extracted terms
    names(input) <- terms
  }

  # 4. Validate list input and write GMT file
  if (is.list(input)) {
    if (is.null(names(input))) {
      stop("List input must have names (used as gene set identifiers)")
    }

    # Open file connection with automatic closure on exit
    con <- file(file_path, "wt")
    on.exit(close(con))

    # Process each gene set
    for (i in seq_along(input)) {
      gene_set_name <- names(input)[i]
      genes <- input[[i]]

      # Final cleaning of gene list
      genes <- na.omit(genes)          # Remove any remaining NA values
      genes <- as.character(genes)     # Ensure character type
      genes <- genes[genes != ""]      # Remove empty strings

      # Skip gene sets with no valid genes
      if (length(genes) == 0) {
        warning(paste("Skipping empty gene set:", gene_set_name))
        next
      }

      # Process description field
      if (is.na(description[1])) {
        gene_set_desc <- ""  # Use empty string if no description provided
      } else {
        # Recycle descriptions if needed
        desc_idx <- ((i - 1) %% length(description)) + 1
        gene_set_desc <- as.character(description[desc_idx])
        # Replace empty descriptions with placeholder
        gene_set_desc <- ifelse(gene_set_desc == "", "No description", gene_set_desc)
      }

      # Construct GMT line and write to file
      gmt_line <- paste(c(gene_set_name, gene_set_desc, genes), collapse = "\t")
      writeLines(gmt_line, con)
    }

    message("GMT file generated successfully: ", file_path)
  } else {
    stop("Input must be either a data frame (with 'term' and 'gene' columns) or a named list")
  }
}


#' geneset2gmt
#'
#' @param geneset A vector of gene names or a path to a txt file (file content should contain one gene per line)
#' @param genesetname A character string specifying the name of the gene set
#' @param description A character string describing the gene set, default is NA
#' @param return The type of object to return, either "data.frame" or "GeneSetCollection"
#' @param folder Path to the output folder, default is current working directory
#' @param filename Name of the output GMT file (must end with .gmt)
#'
#' @return An object of the type specified by the 'return' parameter
#' @export
#'
geneset2gmt <- function(geneset,
                        genesetname,
                        description = NA,
                        return = "data.frame",
                        folder = ".",
                        filename) {
  # Check required parameters
  if (missing(genesetname) || is.null(genesetname) || genesetname == "") {
    stop("Must specify a valid gene set name (genesetname)")
  }

  if (missing(filename)) {
    stop("Must provide an output filename (filename)")
  }

  # Check if filename ends with .gmt
  if (!grepl("\\.gmt$", filename, ignore.case = TRUE)) {
    stop("Filename must end with .gmt")
  }

  # Process input geneset: if it's a txt file path, read genes
  if (length(geneset) == 1 && is.character(geneset)) {
    if (grepl("\\.txt$", geneset, ignore.case = TRUE)) {
      if (!file.exists(geneset)) {
        stop("Specified txt file does not exist: ", geneset)
      }
      # Read genes from txt file (one gene per line)
      geneset <- readLines(geneset, warn = FALSE)
      # Remove empty lines and whitespace
      geneset <- trimws(geneset)
      geneset <- geneset[geneset != ""]
    }
  }

  # Validate gene set
  if (!is.vector(geneset) || length(geneset) == 0) {
    stop("Gene set must be a valid vector or txt file containing genes")
  }

  # Construct input list
  input <- list(geneset)
  names(input) <- genesetname

  # Generate GMT file
  outputGmtFile(
    input = input,
    description = description,
    folder = folder,
    filename = filename
  )

  # Build full file path
  file_path <- file.path(folder, filename)

  # Return appropriate result based on 'return' parameter
  if (return == "data.frame") {
    if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
      stop("clusterProfiler package is required. Please install it with: install.packages('clusterProfiler')")
    }
    gs <- clusterProfiler::read.gmt(file_path)
  } else if (return == "GeneSetCollection") {
    if (!requireNamespace("GSEABase", quietly = TRUE)) {
      stop("GSEABase package is required. Please install it with: BiocManager::install('GSEABase')")
    }
    gs <- GSEABase::getGmt(file_path)
  } else {
    stop("return parameter must be either 'data.frame' or 'GeneSetCollection'")
  }

  return(gs)
}







#' tidy.gmt
#'
#' @param filepath For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param fun "stat" or "merge"
#' @param Source For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param termName For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param addTotal For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param save For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param folder For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param filename For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#'
#' @return data.frame
#' @export tidy.gmt
#'
tidy.gmt <-  function (filepath, fun = "stat", Source = "", termName = NULL,
                       addTotal = FALSE, save = TRUE, folder = ".", filename = "geneset") {
  ifelse(dir.exists(folder), "", dir.create(folder, recursive = T))
  if (dir.exists(filepath)) {
    fp <- dir(filepath, pattern = "gmt$", full.names = T)
  }else if (is.vector(filepath) & sum(grep("gmt$", filepath)) ==
            length(filepath)) {
    fp <- filepath
  }
  stat.gs <- data.frame()
  mergedf <- data.frame()
  gs <- c()
  for (x in fp) {
    ldf <- clusterProfiler::read.gmt(x)
    if(nrow(ldf) != 0){
      mergedf <- rbind(mergedf, ldf)
      gs <- unique(c(gs, ldf$gene))
      df <- data.frame(term = gsub("_", " ", as.character(unique(ldf$term))),
                       Source = Source, `Gene count` = nrow(ldf))
      stat.gs <- rbind(stat.gs, df)
    }

  }
  totallgene <- data.frame(term = "", Source = Source, `Gene count` = paste0(length(unique(gs)),
                                                                             "(unique)"))
  stat.gs <- rbind(stat.gs, totallgene)
  if (!is.null(termName)) {
    colnames(stat.gs)[1] <- termName
  }
  if (fun == "stat") {
    if (save == TRUE) {
      writeLines(unique(gs), con = paste0(folder, "/",
                                          filename, "-stat.unique.txt"))
      write.csv(stat.gs, file = paste0(folder, "/", filename,
                                       "-stat.csv"))
    }
    return(stat.gs)
  }
  else if (fun == "merge") {
    if (save == TRUE) {
      writeLines(unique(gs), con = paste0(folder, "/",
                                          filename, "-merge.unique.txt"))
      outputGmtFile(input = mergedf, description = NA,
                    folder = folder, filename = paste0(filename,
                                                       "-merge.gmt"))
      if (addTotal == TRUE) {
        addmergedf <- data.frame(term = ifelse(is.null(termName),
                                               "Total gene", termName), gene = unique(gs))
        mergedf <- rbind(mergedf, addmergedf)
        outputGmtFile(input = mergedf, description = NA,
                      folder = folder, filename = paste0(filename,
                                                         "-add-unique.gmt"))
      }
    }
    return(mergedf)
  }
}

# tidy.gmt(filepath = "G:/publicData/base_files/GeneSet/Cytoskeleton",
#          fun = "merge",
#          folder = "G:/myProject/Cytoskeleton/data/geneset")




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
#' @export tidyGene.fromeReactome
#'
tidyGene.fromeReactome <- function(filepath
                                   ,fun = "stat"
                                   ,Source = "Reactome"
                                   ,addTotal = TRUE
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
      if(addTotal == TRUE){
        addmergedf <- data.frame(term = ifelse(is.null(termName),"Total gene",termName)
                                 ,gene = unique(gs))
        mergedf <- rbind(mergedf,addmergedf)
        outputGmtFile(input = mergedf
                      ,description = NA
                      ,folder= folder
                      ,filename = paste0(filename,"-add-unique.gmt"))
      }
    }
    return(mergedf)
  }
}

#' read.gmt.to.getGmt
#'
#' @param genesetdf The result of the read.gmt() function
#'
#' @return GeneSetCollection object.
#' @export read.gmt.to.getGmt
#'
read.gmt.to.getGmt <- function(genesetdf){
  # 首先，按照 term 分组
  gene_sets_list <- split(genesetdf$gene, genesetdf$term)
  # 创建 GeneSet 对象列表，过滤掉空字符串
  gene_sets <- lapply(names(gene_sets_list), function(term) {
    genes <- gene_sets_list[[term]]
    genes <- genes[genes != ""]  # 去掉空字符串
    GeneSet(
      geneIds = genes,
      setName = term,
      setIdentifier = term
    )
  })
  # 创建 GeneSetCollection 对象
  gsc <- GeneSetCollection(gene_sets)
  return(gsc)
}


#' GeneSetCollection.to.df
#'
#' @param GeneSetCollection GeneSetCollection object.
#'
#' @return data.frame
#' @export GeneSetCollection.to.df
#'
GeneSetCollection.to.df <- function(GeneSetCollection){
  genesetname <- names(GeneSetCollection)
  geneset = lapply(genesetname, function(x){
    data.frame(term = x,gene = GeneSetCollection[[x]]@geneIds)
  })
  geneset <- do.call(rbind,geneset)
  return(geneset)
}



