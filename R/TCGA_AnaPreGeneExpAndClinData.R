#' Integrate TCGA expression data with clinical data
#'
#' Read TCGA expression data and clinical data from specified folders, perform preprocessing (filter tumor samples, remove duplicates, etc.) and merge them.
#'
#' @param expDataFolder Character string. Path to the folder containing TCGA expression data files (must exist).
#' @param ClinDataFolder Character string. Path to the folder containing clinical data files (must exist).
#' @param cancer Character string. Abbreviation of TCGA cancer type (e.g., "LUAD", "BRCA"), must be a non-empty string.
#' @param ClinDataFrom Character string. Source of clinical data, must be one of: 
#'   "GDCquery_clinic", "UCSCphenotype", "UCSCsurvial". Defaults to "GDCquery_clinic".
#' @param expDataType Character string. Type of expression data, can be "tpm" (transcripts per million) or "count" (read counts). 
#'   Defaults to "tpm".
#' @param genes Character vector. Genes of interest to filter. If NULL (default), no gene filtering is performed (returns all genes passing type filter).
#' @param filter Gene type(s) to retain:
#'   - NULL: retain all gene types (no filtering by type).
#'   - Character string or vector: retain genes where `gene_type` matches these values (e.g., "protein_coding", c("lncRNA", "miRNA")).
#'   Common values include "protein_coding", "lncRNA", "snRNA", "miRNA", "misc_RNA".
#' @param expDatapattern Character string. Regular expression pattern to match expression data files. Defaults to "STARdata.Rdata$".
#'
#' @return Data frame containing merged clinical information and expression data
#' @export
#' @importFrom dplyr arrange
#' @importFrom MedBioInfoCloud splitTCGAmatrix delTCGA_dup_sample
preGeneExpAndClinDate <- function(
  expDataFolder,
  ClinDataFolder,
  cancer,
  ClinDataFrom = "GDCquery_clinic",
  expDataType = "tpm",
  genes = NULL,
  filter = "protein_coding",
  expDatapattern = "STARdata.Rdata$"
) {
  # 参数校验：检查文件夹是否存在
  if (!dir.exists(expDataFolder)) {
    stop("表达数据文件夹不存在: ", expDataFolder)
  }
  if (!dir.exists(ClinDataFolder)) {
    stop("临床数据文件夹不存在: ", ClinDataFolder)
  }
  
  # 参数校验：检查cancer参数合法性
  if (length(cancer) != 1 || !is.character(cancer) || cancer == "") {
    stop("cancer必须为非空字符串（如'LUAD'）")
  }
  
  # 参数校验：检查临床数据来源合法性
  valid_clin_sources <- c("GDCquery_clinic", "UCSCphenotype", "UCSCsurvial")
  if (!ClinDataFrom %in% valid_clin_sources) {
    stop("ClinDataFrom必须为以下之一: ", paste(valid_clin_sources, collapse = ", "))
  }
  
  # 参数校验：检查表达数据类型合法性
  if (!expDataType %in% c("tpm", "count")) {
    stop("expDataType必须为'tpm'或'count'")
  }
  
  # 预处理表达数据：调用外部函数获取表达矩阵
  exp <- preGeneExpMatrixFromTCGA(
    dataFolder = expDataFolder,
    cancer = cancer,
    genes = genes,
    type = expDataType,
    filter = filter
  )
  
  # 筛选肿瘤样本并去重
  turexp <- MedBioInfoCloud::splitTCGAmatrix(data = exp, sample = "Tumor")  # 仅保留肿瘤样本
  turexp <- MedBioInfoCloud::delTCGA_dup_sample(turexp, col_rename = TRUE)  # 去除重复样本
  # 转置矩阵并添加样本名列（作为第一列）
  turexp <- t(turexp) %>% 
    as.data.frame() %>% 
    dplyr::mutate(sample = colnames(turexp), .before = 1)
  
  # 读取并合并临床数据
  flpaths <- dir(ClinDataFolder, full.names = TRUE)  # 获取临床文件夹下所有文件路径
  fl <- flpaths[grep(cancer, flpaths)]  # 筛选匹配当前癌症类型的文件
  
  # 检查临床文件是否存在
  if (length(fl) == 0) {
    stop("未在临床数据文件夹中找到匹配癌症类型 '", cancer, "' 的文件")
  }
  if (!file.exists(fl)) {
    stop("找到的临床文件路径不存在: ", fl)
  }
  
  # 根据不同临床数据来源进行处理
  if (ClinDataFrom == "UCSCphenotype") {
    clindata <- read.delim(fl, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    clindata <- dplyr::arrange(clindata, submitter_id)  # 按submitter_id排序
    clindata <- clindata[!duplicated(clindata$submitter_id), ]  # 去除重复的submitter_id
    mergdata <- merge(clindata, turexp, by.x = "submitter_id", by.y = "sample")  # 合并数据
  } else if (ClinDataFrom == "GDCquery_clinic") {
    load(fl)  # 加载文件（假设文件中包含processed_clindata对象）
    if (!exists("processed_clindata")) {  # 检查加载的对象是否存在
      stop("加载的GDC临床文件中未找到'processed_clindata'对象")
    }
    mergdata <- merge(processed_clindata, turexp, by.x = "submitter_id", by.y = "sample")  # 合并数据
  } else if (ClinDataFrom == "UCSCsurvial") {
    clindata <- read.delim(fl, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    clindata <- dplyr::arrange(clindata, sample)  # 按sample排序
    clindata <- clindata[!duplicated(clindata$X_PATIENT), ]  # 去除重复的X_PATIENT
    mergdata <- merge(clindata, turexp, by.x = "X_PATIENT", by.y = "sample")  # 合并数据
  }
  
  return(mergdata)
}
