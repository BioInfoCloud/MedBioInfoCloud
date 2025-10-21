#' Group samples by gene expression levels and prepare survival analysis data
#' 
#' This function processes input data to group samples into high/low expression groups 
#' based on specified methods (median or quartile), and prepares survival data by 
#' converting time units and removing missing values.
#' 
#' @param data A data frame containing clinical information, gene expression data, and survival data.
#' @param gene Character string. Name of the gene to analyze (must correspond to a column name in `data`).
#' @param surv_cols A character vector of length 2 specifying column names for survival data: first element is survival status (event indicator, 0/1), second is survival time. Defaults to `c("vitalStat", "surTime")`.
#' @param Timeunit A single positive numeric value specifying the conversion factor for survival time units (e.g., use 30 to convert days to months). Defaults to 1 (no conversion), and survival time will be divided by this value.
#' @param group Character string specifying the grouping method for gene expression. Options are "median" (split by median expression) or "quartile" (split by 25th and 75th percentiles, with middle samples excluded). Defaults to "quartile".
#' 
#' @return A data frame with columns: survival status, survival time (converted), gene expression, and a `group` column (factor with levels "Low" and "High").
#' @export
#' @importFrom survival Surv survfit coxph
#' @importFrom survminer ggsurvplot
#' @importFrom ggplot2 theme_classic element_text
#' @importFrom dplyr filter mutate select
#' @importFrom stats quantile median as.formula
#' @importFrom rlang sym
geneHighLowExpGroupData <- function(
  data,
  gene,
  surv_cols = c("vitalStat", "surTime"),
  Timeunit = 1,
  group = c("quartile", "median")
) {
  # 参数验证：基础类型检查
  if (!inherits(data, "data.frame")) {
    stop("'data' must be a data frame.")
  }
  if (!is.character(gene) || length(gene) != 1 || !gene %in% colnames(data)) {
    stop(
      "'gene' must be a single character string matching a column name in 'data'."
    )
  }
  if (!is.character(surv_cols) || length(surv_cols) != 2) {
    stop(
      "'surv_cols' must be a character vector of length 2 (survival status, survival time)."
    )
  }
  if (!is.numeric(Timeunit) || Timeunit <= 0 || length(Timeunit) != 1) {
    stop("'Timeunit' must be a single positive numeric value.")
  }
  
  # 限定分组方式为预设值（避免无效输入）
  group <- match.arg(group)
  
  # 解析生存数据列名（状态列和时间列）
  vital_stat_col <- surv_cols[1]  # 生存状态列名（如"vitalStat"）
  surv_time_col <- surv_cols[2]   # 生存时间列名（如"surTime"）
  
  # 检查必要列是否存在于输入数据中
  required_cols <- c(vital_stat_col, surv_time_col, gene)
  missing_cols <- setdiff(required_cols, colnames(data))
  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns in 'data': ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  # 验证生存状态和时间的类型（必须为数值型）
  if (!is.numeric(data[[vital_stat_col]])) {
    stop(paste0("'", vital_stat_col, "' must be numeric (0/1)."))
  }
  if (!is.numeric(data[[surv_time_col]])) {
    stop(paste0("'", surv_time_col, "' must be numeric (survival time)."))
  }
  
  # 数据预处理：提取必要列并移除含缺失值的样本
  analysis_data <- data %>%
    dplyr::select(all_of(required_cols)) %>%
    dplyr::filter(
      !is.na(!!sym(gene)),           # 排除基因表达缺失的样本
      !is.na(!!sym(vital_stat_col)), # 排除生存状态缺失的样本
      !is.na(!!sym(surv_time_col))   # 排除生存时间缺失的样本
    )
  
  # 检查有效样本量（至少需要2个样本进行分析）
  if (nrow(analysis_data) < 2) {
    stop("Insufficient valid samples (n < 2) after removing missing values.")
  }
  
  # 转换生存时间单位（如将天数转换为月数：除以30）
  analysis_data <- analysis_data %>%
    dplyr::mutate(!!surv_time_col := !!sym(surv_time_col) / Timeunit)
  
  # 按指定方法对基因表达进行高低分组
  if (group == "median") {
    # 按中位数分组：高于中位数为"High"，低于等于为"Low"
    analysis_data <- analysis_data %>%
      dplyr::mutate(
        median_val = median(!!sym(gene), na.rm = TRUE),  # 计算基因表达中位数
        group = ifelse(!!sym(gene) > median_val, "High", "Low"),  # 分组
        group = factor(group, levels = c("Low", "High"))  # 转换为因子（指定顺序）
      ) %>%
      dplyr::select(-median_val)  # 移除临时变量
  } else {
    # 按四分位数分组：仅保留低于25%分位数（Low）和高于75%分位数（High）的样本
    analysis_data <- analysis_data %>%
      dplyr::mutate(
        q25 = quantile(!!sym(gene), 0.25, na.rm = TRUE),  # 计算25%分位数
        q75 = quantile(!!sym(gene), 0.75, na.rm = TRUE)   # 计算75%分位数
      ) %>%
      dplyr::filter(!!sym(gene) <= q25 | !!sym(gene) >= q75) %>%  # 排除中间50%样本
      dplyr::mutate(
        group = ifelse(!!sym(gene) >= q75, "High", "Low"),  # 分组
        group = factor(group, levels = c("Low", "High"))    # 转换为因子（指定顺序）
      ) %>%
      dplyr::select(-q25, -q75)  # 移除临时变量
  }
  
  # 检查分组后的有效样本量
  if (nrow(analysis_data) < 2) {
    stop("Insufficient valid samples (n < 2) after quartile grouping.")
  } else {
    return(analysis_data)
  }
}

# #' @param title Character string. Main title of the survival curve plot.
# #' @param xtitle Character string. Title for the x-axis (should match the time unit after conversion).
# #' @param palette A vector of length 2 specifying colors for "Low" and "High" expression groups, respectively. Defaults to `c("#E7B800", "#2E9FDF")`.
# #' @param Timeunit Numeric value. Conversion factor for survival time (e.g., use 365 to convert days to years). Defaults to 1 (no conversion).
# SurvivalAnalysis = function(
#   analysis_data,
#   data_cols = c("vitalStat", "surTime", "group"), # 生存状态，生存时间，和分组的列名，要求顺序不要改变。
#   title,
#   xtitle,
#   palette = c("#E7B800", "#2E9FDF")
# ) {
#   adat.surv.obj <- Surv(
#     time = analysis_data[, data_cols[2]],
#     event = analysis_data[,data_cols[1]]
#   )
#   # #生存率的比较
#   adat.surv.treat <- survfit(adat.surv.obj ~ group, data = analysis_data)

#   # 关键修复：用字符串构建公式，避免非标准评估错误
#   # 公式格式：Surv(生存时间列, 生存状态列) ~ group
#   surv_formula_str <- paste0(
#     "Surv(",
#     data_cols[2],
#     ", ",
#     data_cols[1],
#     ") ~ group"
#   )
#   surv_formula <- as.formula(surv_formula_str) # 转换为公式对象

#   # Kaplan-Meier生存分析（使用构建的公式）
#   km_fit <- survival::survfit(
#     formula = surv_formula,
#     data = analysis_data
#   )
#   # Cox比例风险模型（使用相同公式）
#   cox_model <- survival::coxph(
#     formula = surv_formula,
#     data = analysis_data
#   )
#   cox_summary <- summary(cox_model)
#   # 提取HR统计量
#   coef_name <- grep("group", rownames(cox_summary$coefficients), value = TRUE)[
#     1
#   ]
#   if (is.na(coef_name)) {
#     stop("Could not extract coefficient for group comparison from Cox model.")
#   }
#   hr_stats <- data.frame(
#     HR = cox_summary$coefficients[coef_name, "exp(coef)"],
#     HR_low = cox_summary$conf.int[coef_name, "lower .95"],
#     HR_high = cox_summary$conf.int[coef_name, "upper .95"],
#     p_value = cox_summary$coefficients[coef_name, "Pr(>|z|)"],
#     row.names = NULL
#   )

#   # 绘制生存曲线（公式与数据框匹配）
#   surv_plot <- survminer::ggsurvplot(
#     fit = survfit(
#       Surv(
#         time = analysis_data[, data_cols[2]],
#         event = analysis_data[, data_cols[1]]
#       ) ~
#         group,
#       data = analysis_data
#     ),
#     data = analysis_data,
#     pval = TRUE,
#     title = cancer,
#     xlab = "Time(months)",
#     palette = c("#E7B800", "#2E9FDF"),
#     legend.title = gene,
#     legend.labs = c("Low", "High"),
#     ggtheme = theme_classic() +
#       theme(
#         plot.title = element_text(hjust = 0.5),
#         legend.text = element_text(size = 10, color = "black"),
#         legend.title = element_text(size = 10, color = "black"),
#         axis.text = element_text(size = 10, color = "black"),
#         axis.title = element_text(size = 10, color = "black")
#       )
#   )
#   # 返回结果
#   list(
#     sample_info = table(analysis_data$group),
#     surv_fit = km_fit,
#     cox_result = cox_summary,
#     hr_stats = hr_stats,
#     surv_plot = surv_plot
#   )
# }

