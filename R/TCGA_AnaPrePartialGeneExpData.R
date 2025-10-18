#' Retrieve sample number information for TCGA cancer types
#'
#' This function returns a data frame containing the number of normal and tumor samples
#' across 33 TCGA cancer types, with optional filtering by specific cancer type(s).
#'
#' @param cancer A character string specifying the cancer type(s) to retrieve. Use "all" (default)
#'   to get information for all cancer types. Partial matches are supported via regex.
#'
#' @return A data frame with columns:
#'   \item{cancer}{TCGA cancer type identifier (e.g., "TCGA-BRCA")}
#'   \item{nNormal}{Number of normal samples}
#'   \item{nTumor}{Number of tumor samples}
TCGASamoleNumInfo <- function(cancer = "all") {
  TCGAsampleNum <- data.frame(
    cancer = c(
      "TCGA-BRCA",
      "TCGA-UCEC",
      "TCGA-LAML",
      "TCGA-CHOL",
      "TCGA-ACC",
      "TCGA-DLBC",
      "TCGA-CESC",
      "TCGA-ESCA",
      "TCGA-KICH",
      "TCGA-HNSC",
      "TCGA-COAD",
      "TCGA-BLCA",
      "TCGA-KIRP",
      "TCGA-GBM",
      "TCGA-MESO",
      "TCGA-SARC",
      "TCGA-PCPG",
      "TCGA-READ",
      "TCGA-LUAD",
      "TCGA-LUSC",
      "TCGA-PAAD",
      "TCGA-LIHC",
      "TCGA-KIRC",
      "TCGA-LGG",
      "TCGA-SKCM",
      "TCGA-PRAD",
      "TCGA-STAD",
      "TCGA-OV",
      "TCGA-THYM",
      "TCGA-UCS",
      "TCGA-UVM",
      "TCGA-TGCT",
      "TCGA-THCA"
    ),
    nNormal = as.integer(c(
      113,
      35,
      0,
      9,
      0,
      0,
      3,
      11,
      25,
      44,
      41,
      19,
      32,
      5,
      0,
      2,
      3,
      10,
      59,
      49,
      4,
      50,
      72,
      0,
      1,
      52,
      32,
      0,
      2,
      0,
      0,
      0,
      59
    )),
    nTumor = as.integer(c(
      1113,
      554,
      150,
      35,
      79,
      48,
      306,
      163,
      65,
      504,
      480,
      412,
      291,
      169,
      87,
      263,
      184,
      167,
      539,
      502,
      179,
      374,
      541,
      532,
      472,
      501,
      375,
      381,
      120,
      57,
      80,
      156,
      512
    )),
    check.names = FALSE
  )
  TCGAsampleNum = arrange(TCGAsampleNum, nNormal)
  if (cancer == "all") {
    return(TCGAsampleNum)
  } else {
    allcancer = TCGAsampleNum$cancer
    mapacncer = allcancer[grepl(paste(cancer, collapse = "|"), allcancer)]
    return(TCGAsampleNum[TCGAsampleNum$cancer %in% mapacncer, ])
  }
}

#' Select cancer types based on sample number thresholds
#'
#' This function filters TCGA cancer types based on the number of normal or tumor samples
#' using a specified threshold and direction.
#'
#' @param sample A character string specifying the sample type to filter by, either "nNormal" (default) or "nTumor"
#' @param n An integer specifying the threshold sample number
#' @param direction An integer indicating the filtering direction: 1 for greater than or equal to n,
#'   -1 for less than n
#'
#' @return A character vector of cancer type identifiers meeting the criteria
selectCancerBaseSamNum = function(sample = "nNormal", n = 3, direction = 1) {
  numinfo = TCGASamoleNumInfo()
  if (direction == 1) {
    return(numinfo[numinfo[, sample] >= n, ][, 1])
  }
  if (direction == -1) return(numinfo[numinfo[, sample] < n, ][, 1])
}


#' Prepare gene expression data for a single cancer type
#'
#' This function processes gene expression data for a specific cancer type, including filtering,
#' log2 transformation, and organization into unpaired and paired (if available) datasets.
#'
#' @param dataFolder A character string specifying the path to the folder containing expression data files
#' @param cancer A character string specifying the TCGA cancer type to process
#' @param genes A character vector of gene names to include in the analysis
#' @param returnPaired A logical value indicating whether to return paired sample data (default: TRUE)
#' @param pattern A character string specifying the regular expression pattern to match data files (default: "STARdata.Rdata$")
#' @param nNormal An integer specifying the minimum number of normal samples required (default: 10)
#' @param nPaired An integer specifying the minimum number of paired samples required (only used if returnPaired = TRUE, default: 3)
#'
#' @return A list containing:
#'   \item{unpaired}{A data frame with unpaired expression data}
#'   \item{paired}{A data frame with paired expression data (or NULL if unavailable)}
preGeneExpDatIn1Cancer <- function(
    dataFolder, # 数据所在的文件夹
    cancer, # 癌症类型
    genes, # 基因向量
    returnPaired = TRUE, # 是否返回配对样本的数据，默认返回
    pattern = "STARdata.Rdata$", # 需要匹配的文件正则表达式
    nNormal = 10, # 至少包含的样本数
    nPaired = 3 # 配对的样本数，仅仅当 returnPaired = TRUE是使用，应该小于等于nNormal。
) {
  sinfo = TCGASamoleNumInfo(cancer = cancer)
  if (nrow(sinfo) == 1 && sinfo$nNormal >= nNormal) {
    FilePath <- dir(dataFolder, pattern, full.names = T)
    load(FilePath[grep(cancer, FilePath)]) # STARdata
    tpm <- STARdata[["tpm"]]
    ##正常组织样本ID
    SamN <- TCGAbiolinks::TCGAquery_SampleTypes(
      barcode = colnames(tpm)[-c(1:3)],
      typesample = c("NT", "NB", "NBC", "NEBV", "NBM")
    )
    tpm <- filterGeneTypeExpr(expr = tpm, fil_col = "gene_type", filter = FALSE)
    ##判断基因是否在表达矩阵中
    gs <- intersect(genes, rownames(tpm))
    if (length(gs) != 0) {
      ##肿瘤组织样本ID
      SamT <- setdiff(colnames(tpm), SamN)
      ###去除重复样本
      nor_exp <- delTCGA_dup_sample(tpm[, SamN], col_rename = T)
      tur_exp <- delTCGA_dup_sample(tpm[, SamT], col_rename = T)
      ###============非配对样本==
      ##构建数据框
      nor_gsExp <- log2(data.frame(nor_exp[gs, ], check.names = T) + 1) %>%
        t() %>%
        as.data.frame() %>%
        mutate(
          Sample = rep(
            paste0("Normal(n=", ncol(nor_exp), ")"),
            ncol(nor_exp)
          ),
          .before = 1
        ) %>%
        melt(id.vars = "Sample")

      tur_gsExp <- log2(data.frame(tur_exp[gs, ], check.names = T) + 1) %>%
        t() %>%
        as.data.frame() %>%
        mutate(
          Sample = rep(paste0("Tumor(n=", ncol(tur_exp), ")"), ncol(tur_exp)),
          .before = 1
        ) %>%
        melt(id.vars = "Sample")
      data <- rbind(nor_gsExp, tur_gsExp)
      colnames(data) <- c("Sample", "gene", "exp")
      unPaired.data = dplyr::mutate(
        data,
        cancer = cancer,
        group = gsub("\\(.*?\\)", "", data$Sample),
        .before = 1
      )

      if (returnPaired) {
        ##===========配对样本=========
        id <- intersect(colnames(nor_exp), colnames(tur_exp))
        if (length(id) > nPaired) {
          paired_nor_gsExp <- nor_exp[gs, id]
          paired_tur_gsExp <- tur_exp[gs, id]

          paired_nor_gsExp <- log2(
            data.frame(paired_nor_gsExp, check.names = T) + 1
          ) %>%
            t() %>%
            as.data.frame() %>%
            mutate(
              Sample = rep("Normal", ncol(paired_nor_gsExp)),
              id = id,
              .before = 1
            ) %>%
            melt(id.vars = c("Sample", "id"))

          paired_tur_gsExp <- log2(
            data.frame(paired_tur_gsExp, check.names = T) + 1
          ) %>%
            t() %>%
            as.data.frame() %>%
            mutate(
              Sample = rep("Tumor", ncol(paired_tur_gsExp)),
              id = id,
              .before = 1
            ) %>%
            melt(id.vars = c("Sample", "id"))
          paired_data <- rbind(paired_nor_gsExp, paired_tur_gsExp)
          head(paired_data)
          colnames(paired_data) <- c("Sample", "id", "gene", "exp")
          paired_data = mutate(paired_data, cancer = cancer, .before = 1)
        } else {
          paired_data = NULL
        }
      } else {
        paired_data = NULL
      }
      datalist = list(unpaired = unPaired.data, paired = paired_data)
      return(datalist)
    } else {
      stop(
        "genes 中的基因不在表达数据中为找到，请检查基因名称，可考虑基因的其他别名"
      )
    }
  }
}


#' Prepare gene expression data across multiple cancer types (pan-cancer)
#'
#' This function processes gene expression data across multiple TCGA cancer types, aggregating
#' results from individual cancer type processing.
#'
#' @param dataFolder A character string specifying the path to the folder containing expression data files
#' @param genes A character vector of gene names to include in the analysis
#' @param returnPaired A logical value indicating whether to return paired sample data (default: TRUE)
#' @param pattern A character string specifying the regular expression pattern to match data files (default: "STARdata.Rdata$")
#' @param nNormal An integer specifying the minimum number of normal samples required for inclusion of a cancer type (default: 10)
#' @param nPaired An integer specifying the minimum number of paired samples required (only used if returnPaired = TRUE, default: 3)
#'
#' @return A list containing:
#'   \item{unpaired}{A data frame with aggregated unpaired expression data across cancer types}
#'   \item{paired}{A data frame with aggregated paired expression data across cancer types (or NULL if unavailable)}
preGeneExpDatInPanCancer <- function(
    dataFolder, # 数据所在的文件夹
    genes, # 基因向量
    returnPaired = TRUE, # 是否返回配对样本的数据，默认返回
    pattern = "STARdata.Rdata$", # 需要匹配的文件正则表达式
    nNormal = 10, # 至少包含的样本数
    nPaired = 3 # 配对的样本数，仅仅当 returnPaired = TRUE是使用，应该小于等于nNormal。
) {
  cancers = selectCancerBaseSamNum(n = nNormal)
  datalist = lapply(cancers, function(x) {
    preGeneExpDatIn1Cancer(
      dataFolder = dataFolder,
      genes = genes,
      returnPaired = returnPaired,
      pattern = pattern,
      cancer = x,
      nNormal = nNormal
    )
  })
  unpairedDat = do.call(
    rbind,
    lapply(datalist, function(x) {
      x[["unpaired"]]
    })
  )
  if (returnPaired) {
    pairedDat = do.call(
      rbind,
      lapply(datalist, function(x) {
        x[["paired"]]
      })
    )
  } else {
    pairedDat = NULL
  }

  data = list(unpaired = unpairedDat, paired = pairedDat)
  return(data)
}


#' Prepare gene expression data for a single cancer type with GTEx normal sample mapping
#'
#' This function processes gene expression data for a specific cancer type, using GTEx data as
#' normal tissue references when TCGA normal samples are limited.
#'
#' @param dataFolder A character string specifying the path to the folder containing TCGA-GTEx merged data files
#' @param cancer A character string specifying the TCGA cancer type to process
#' @param genes A character vector of gene names to include in the analysis
#' @param pattern A character string specifying the regular expression pattern to match data files (default: "GTEX_norm_tpm.Rdata$")
#'
#' @return A data frame containing merged expression data with GTEx normal samples and TCGA tumor samples
preGeneExpDatIn1CancerMapGTEx = function(
    dataFolder, # 数据所在的文件夹
    cancer,
    genes, # 基因向量
    pattern = "GTEX_norm_tpm.Rdata$" # 需要匹配的文件正则表达式
) {
  haveCancer = c(
    "TCGA-ACC",
    "TCGA-BLCA",
    "TCGA-BRCA",
    "TCGA-CESC",
    "TCGA-COAD",
    "TCGA-ESCA",
    "TCGA-GBM",
    "TCGA-KICH",
    "TCGA-KIRC",
    "TCGA-KIRP",
    "TCGA-LGG",
    "TCGA-LIHC",
    "TCGA-LUAD",
    "TCGA-LUSC",
    "TCGA-OV",
    "TCGA-PAAD",
    "TCGA-PRAD",
    "TCGA-SKCM",
    "TCGA-STAD",
    "TCGA-TGCT",
    "TCGA-UCEC",
    "TCGA-UCS"
  )
  if (length(cancer) == 1 & cancer %in% haveCancer) {
    #=======================
    TCGA_TARGET_GTEx_FilePath <- dir(dataFolder, pattern, full.names = T)
    load(TCGA_TARGET_GTEx_FilePath[grep(cancer, TCGA_TARGET_GTEx_FilePath)]) # GTEX_TCGA_norm_tpm
    nor_exp <- GTEX_TCGA_norm_tpm[, grep("GTEX", colnames(GTEX_TCGA_norm_tpm))]
    tur_exp <- GTEX_TCGA_norm_tpm[, grep("TCGA", colnames(GTEX_TCGA_norm_tpm))]
    colnames(tur_exp)
    ##构建数据框
    nor_gsExp <- log2(nor_exp[genes, ] + 1) %>%
      t() %>%
      as.data.frame() %>%
      mutate(
        Sample = rep(paste0("Normal(n=", ncol(nor_exp), ")"), ncol(nor_exp)),
        .before = 1
      ) %>%
      melt(id.vars = "Sample")
    nor_gsExp = mutate(nor_gsExp, group = "Normal", .before = 1)
    tur_gsExp <- log2(data.frame(tur_exp[genes, ])) %>%
      t() %>%
      as.data.frame() %>%
      mutate(
        Sample = rep(paste0("Tumor(n=", ncol(tur_exp), ")"), ncol(tur_exp)),
        .before = 1
      ) %>%
      melt(id.vars = "Sample")
    tur_gsExp = mutate(tur_gsExp, group = "Tumor", .before = 1)
    data <- rbind(nor_gsExp, tur_gsExp)
    head(data)
    colnames(data) <- c("group", "Sample", "gene", "exp")
    data = mutate(data, cancer = cancer, .before = 1)
    return(data)
  }
}


#' Prepare pan-cancer gene expression data with GTEx normal sample mapping
#'
#' This function processes gene expression data across multiple cancer types, using GTEx data as
#' normal tissue references for cancer types with limited TCGA normal samples.
#'
#' @param dataFolder A character string specifying the path to the folder containing TCGA-GTEx merged data files
#' @param genes A character vector of gene names to include in the analysis
#' @param pattern A character string specifying the regular expression pattern to match data files (default: "GTEX_norm_tpm.Rdata$")
#' @param maxnNormalinTCGA A numeric value specifying the maximum number of normal samples in TCGA to qualify for GTEx mapping (default: Inf)
#'
#' @return A data frame containing aggregated expression data across cancer types with GTEx normal samples
preGeneExpDatInPanCancerMapGTEx <- function(
    dataFolder, # 数据所在的文件夹
    genes, # 基因向量
    pattern = "GTEX_norm_tpm.Rdata$", # 需要匹配的文件正则表达式
    maxnNormalinTCGA = Inf # 在TCGA中最大的正常样本数
) {
  haveCancer = c(
    "TCGA-ACC",
    "TCGA-BLCA",
    "TCGA-BRCA",
    "TCGA-CESC",
    "TCGA-COAD",
    "TCGA-ESCA",
    "TCGA-GBM",
    "TCGA-KICH",
    "TCGA-KIRC",
    "TCGA-KIRP",
    "TCGA-LGG",
    "TCGA-LIHC",
    "TCGA-LUAD",
    "TCGA-LUSC",
    "TCGA-OV",
    "TCGA-PAAD",
    "TCGA-PRAD",
    "TCGA-SKCM",
    "TCGA-STAD",
    "TCGA-TGCT",
    "TCGA-UCEC",
    "TCGA-UCS"
  )
  projInTCGA = selectCancerBaseSamNum(n = maxnNormalinTCGA, direction = -1)
  cancers = intersect(haveCancer, projInTCGA)
  if (length(cancers) > 1) {
    datalist = lapply(cancers, function(cancer) {
      preGeneExpDatIn1CancerMapGTEx(
        dataFolder = dataFolder,
        cancer = cancer,
        genes = genes,
        pattern = pattern
      )
    })
    data = do.call(rbind, datalist)
    return(data)
  }
}
