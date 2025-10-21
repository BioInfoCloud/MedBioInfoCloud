#' Plot Gene Expression Across Multiple Genes in One Cancer Type
#'
#' Generate boxplots for expression levels of multiple genes in a single cancer type, with statistical comparisons between groups.
#'
#' @param boxdata A data frame containing gene expression data. Must include columns: `gene` (gene names), `exp` (expression values, log2(TPM + 1)), `group` (comparison groups), and `cancer` (cancer type).
#' @param color A character vector of length 2 specifying fill colors for the two groups in the boxplots. Defaults to `c("#E7B800", "#2E9FDF")`.
#' @param facet A logical value indicating whether to facet the plot by cancer type. If `TRUE` (default), facets by `cancer`; automatically set to `FALSE` if only one cancer type is present.
#'
#' @return A ggplot object showing boxplots of gene expression across multiple genes, with group comparisons.
#' @importFrom ggplot2 ggplot aes geom_boxplot labs ylim scale_color_manual scale_fill_manual theme_bw theme element_blank element_text facet_wrap
#' @importFrom ggpubr stat_compare_means
#' @importFrom dplyr pull
#' @keywords visualization expression
multGenesBase1cancer = function(
    boxdata,
    color = c("#E7B800", "#2E9FDF"),
    facet = TRUE
) {
  # 核心：使用分位数计算y轴范围，排除极端离群值
  y_limits <- boxdata %>%
    pull(exp) %>% # 提取y值列
    quantile(probs = c(0.01, 0.99), na.rm = TRUE) # 取1%和99%分位数，过滤极端值

  # 计算适当的缓冲（基于分位数范围的5%）
  buffer <- diff(y_limits) * 0.05
  y_lower <- y_limits[1] - buffer # 下限
  y_upper <- y_limits[2] + buffer # 上限
  boxdata$cancer = gsub("TCGA-", "", boxdata$cancer)
  p = ggplot(boxdata, aes(x = gene, y = exp, color = group, fill = group)) +
    # 箱线图（调整线条粗细和离群点样式）
    geom_boxplot(
      alpha = 1,
      lwd = 0.5, # 适当增加线条粗细更清晰
      outlier.size = 1,
      outlier.colour = "white",
      width = 0.3 # 设置箱型图的宽度，可根据需要调整该数值
    ) +
    # 添加统计检验（组间比较）
    ggpubr::stat_compare_means(
      label = "p.signif", # 显示显著性标记（*/*/***）
      method = "wilcox.test", # 箱线图常用非参数检验，可根据数据类型调整
      vjust = 1,
      hide.ns = FALSE, # 显示不显著(ns)的标注，如需隐藏可设为TRUE
      show.legend = FALSE # 添加这一行，不显示统计检验的图例
    ) +
    # 坐标轴标签
    labs(y = expression(log[2](TPM + 1))) +
    ylim(c(y_lower, y_upper)) +
    # 颜色设置（确保与填充色区分开更易读）
    scale_color_manual(values = c("#1A1A1A", "#1A1A1A")) + # 修正原代码颜色重复问题
    scale_fill_manual(values = color) +
    # 主题设置
    theme_bw() + # 基于白色背景的主题，便于去除网格线
    theme(
      # 移除网格线
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      # 图例位置在顶部
      legend.position = "top",
      # 图例标题不显示，文本样式调整
      legend.title = element_blank(),
      legend.text = element_text(size = 12, face = "bold", colour = "#1A1A1A"),
      # 坐标轴文本样式
      axis.text.x = element_text(
        face = "bold",
        size = 10,
        angle = 45,
        hjust = 1,
        colour = "#1A1A1A"
      ),
      axis.text.y = element_text(face = "bold", size = 10, colour = "#1A1A1A"),
      # 坐标轴标题样式（x轴标题不显示）
      axis.title.y = element_text(size = 12, face = "bold", colour = "#1A1A1A"),
      axis.title.x = element_blank(),
      strip.background = element_rect(fill = "white", color = NA) # color=NA去除边框
    )
  if (length(unique(boxdata$cancer)) == 1) {
    facet = FALSE # 只有一个癌症时，没必要画分面图
  }
  if (facet) {
    p = p + facet_wrap(~cancer, nrow = 2)
  }
  return(p)
}


#' Plot Gene Expression Across Multiple Cancers for One Gene
#'
#' Generate boxplots for expression levels of a single gene across multiple cancer types, with statistical comparisons between groups.
#'
#' @param boxdata A data frame containing gene expression data. Must include columns: `cancer` (cancer types), `exp` (expression values, log2(TPM + 1)), `group` (comparison groups), and `gene` (target gene name).
#' @param color A character vector of length 2 specifying fill colors for the two groups in the boxplots. Defaults to `c("#E7B800", "#2E9FDF")`.
#'
#' @return A ggplot object showing boxplots of gene expression across multiple cancers, with group comparisons.
#' @importFrom ggplot2 ggplot aes geom_boxplot labs ylim scale_color_manual scale_fill_manual theme_bw theme element_blank element_text facet_wrap scale_x_discrete
#' @importFrom ggpubr stat_compare_means
#' @importFrom dplyr pull
#' @keywords visualization expression
multCancerBase1gene = function(
    boxdata,
    color = c("#E7B800", "#2E9FDF")
) {
  # 核心：使用分位数计算y轴范围，排除极端离群值
  y_limits <- boxdata %>%
    pull(exp) %>% # 提取y值列
    quantile(probs = c(0.01, 0.99), na.rm = TRUE) # 取1%和99%分位数，过滤极端值

  # 计算适当的缓冲（基于分位数范围的5%）
  buffer <- diff(y_limits) * 0.05
  y_lower <- y_limits[1] - buffer # 下限
  y_upper <- y_limits[2] + buffer # 上限
  boxdata$cancer = gsub("TCGA-", "", boxdata$cancer)
  p = ggplot(boxdata, aes(x = cancer, y = exp, color = group, fill = group)) +
    # 箱线图（调整线条粗细和离群点样式）
    geom_boxplot(
      alpha = 1,
      lwd = 0.5, # 适当增加线条粗细更清晰
      outlier.size = 1,
      outlier.colour = "white",
      width = 0.3, # 设置箱型图的宽度，可根据需要调整该数值
      position = position_dodge(0.5)
    ) +
    # 添加统计检验（组间比较）
    ggpubr::stat_compare_means(
      label = "p.signif", # 显示显著性标记（*/*/***）
      method = "wilcox.test", # 箱线图常用非参数检验，可根据数据类型调整
      vjust = 1,
      hide.ns = FALSE, # 显示不显著(ns)的标注，如需隐藏可设为TRUE
      show.legend = FALSE # 添加这一行，不显示统计检验的图例
    ) +
    scale_x_discrete(
      expand = expansion(mult = 0.05, add = 0.1) # 控制分组两侧的留白，值越大间距越宽
    ) +
    # 坐标轴标签
    labs(y = expression(log[2](TPM + 1))) +
    ylim(c(y_lower, y_upper)) +
    # 颜色设置（确保与填充色区分开更易读）
    scale_color_manual(values = c("#1A1A1A", "#1A1A1A")) + # 修正原代码颜色重复问题
    scale_fill_manual(values = color) +
    # 主题设置
    theme_bw() + # 基于白色背景的主题，便于去除网格线
    theme(
      # 移除网格线
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      # 图例位置在顶部
      legend.position = "top",
      # 图例标题不显示，文本样式调整
      legend.title = element_blank(),
      legend.text = element_text(size = 12, face = "bold", colour = "#1A1A1A"),
      # 坐标轴文本样式
      axis.text.x = element_text(
        face = "bold",
        size = 10,
        colour = "#1A1A1A"
      ),
      axis.text.y = element_text(face = "bold", size = 10, colour = "#1A1A1A"),
      # 坐标轴标题样式（x轴标题不显示）
      axis.title.y = element_text(size = 12, face = "bold", colour = "#1A1A1A"),
      axis.title.x = element_blank(),
      strip.background = element_rect(fill = "white", color = NA) # color=NA去除边框
    ) +
    facet_wrap(~gene, nrow = length(unique(boxdata$gene)))
  return(p)
}


#' Plot Gene Expression in a Single Cancer Type (Boxplot and Violin Plot)
#'
#' Generate both boxplot and violin plot for expression levels of a single gene in a specific cancer type, with statistical comparisons between two sample groups.
#'
#' @param data A data frame containing gene expression data. Must include columns: `gene` (gene names), `cancer` (cancer types), `Sample` (sample groups, e.g., "Normal" vs "Tumor"), `exp` (expression values, log2(TPM + 1)), and `id` (sample identifiers).
#' @param cancer A character string specifying the target cancer type to plot.
#' @param gene A character string specifying the target gene to plot.
#' @param color A character vector of length 2 specifying fill colors for the two sample groups. Defaults to `c("#E7B800", "#2E9FDF")`.
#'
#' @return A list containing two ggplot objects: `violin` (violin plot with scatter points and boxplot overlay) and `boxplot` (standard boxplot with significance annotations).
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_violin geom_point labs scale_fill_manual scale_color_manual theme_classic theme element_rect element_text
#' @keywords visualization expression
singleGeneExpIn1cancer = function(
    data,
    cancer,
    gene,
    color = c("#E7B800", "#2E9FDF")
) {
  ldat <- data[data$gene == gene, ]
  ldat <- ldat[ldat$cancer == cancer, ]
  compaired = list(unique(ldat$Sample))
  theme = theme_classic() +
    theme(
      panel.background = element_rect(
        fill = "white",
        colour = "black",
        linewidth = 0.25
      ),
      plot.title = element_text(hjust = 0.5),
      axis.line = element_line(colour = "black", linewidth = 0.25),
      axis.title = element_text(
        size = 10,
        face = "plain",
        color = "black"
      ),
      axis.text.x = element_text(
        face = "plain",
        colour = "black",
        vjust = 1
      ),
      axis.title.x = element_blank(),
      axis.title.y = element_text(
        size = 12,
        face = "plain",
        color = "black"
      ),
      axis.text = element_text(size = 12, color = "black"),
      legend.position = "none"
    )
  boxfig = ggplot(ldat, aes(Sample, exp, fill = Sample)) +
    geom_boxplot(
      aes(fill = Sample),
      notch = FALSE,
      position = position_dodge(width = 0.8), # 调整这个值控制两个箱子间距
      outlier.alpha = 1,
      width = 0.4
    ) +
    scale_x_discrete(expand = expansion(mult = 0.2, add = 0.2)) +
    scale_fill_manual(values = color) +
    ggpubr::geom_signif(
      comparisons = compaired,
      step_increase = 0.1,
      map_signif_level = T,
      margin_top = 0.05, # 减小这个值，使标注整体下移
      ## vjust = 1,
      test = "wilcox.test"
    ) +
    labs(
      y = paste0("The expression of ", gene, "\nlog2(TPM + 1)"),
      title = cancer
    ) +
    theme
  violinfig = ggplot(ldat, aes(x = Sample, y = exp)) +
    #violinplot本体
    geom_violin(
      aes(fill = Sample),
      color = NA,
      alpha = 0.6,
      width = 0.7,
      trim = TRUE,
      scale = "width"
    ) +
    #样本散点
    geom_point(
      aes(color = Sample, fill = Sample),
      show.legend = F,
      position = position_jitter(seed = 123456, width = 0.2),
      shape = 21,
      size = 2
    ) +
    #boxplot本体
    geom_boxplot(
      aes(fill = Sample),
      width = 0.5,
      size = 0.5,
      alpha = 0.6,
      outlier.shape = NA
    ) + #离群点不单独显示
    ggpubr::geom_signif(
      comparisons = compaired,
      step_increase = 0.1,
      map_signif_level = T,
      margin_top = 0.2,
      test = "wilcox.test"
    ) +
    scale_fill_manual(
      values = color
    ) +
    scale_color_manual(
      values = color
    ) +
    labs(
      y = paste0("The expression of ", gene, "\nlog2(TPM + 1)"),
      title = cancer
    ) +
    theme
  return(list(violin = violinfig, boxplot = boxfig))
}


#' Plot Pan-Cancer Gene Expression Across Multiple Cancers
#'
#' Generate boxplots with scatter points for expression levels of a single gene across multiple cancer types, with statistical comparisons between groups (paired or unpaired).
#'
#' @param data A data frame containing gene expression data. Must include columns: `gene` (gene names), `cancer` (cancer types), `group` (comparison groups), `exp` (expression values, log2(TPM + 1)).
#' @param gene A character string specifying the target gene to plot.
#' @param paired A logical value indicating whether the samples are paired (e.g., matched normal-tumor pairs). If `TRUE`, groups are treated as paired.
#' @param color A character vector of length 2 specifying colors for the two groups. Defaults to `c("#E7B800", "#2E9FDF")`.
#'
#' @return A ggplot object showing pan-cancer expression patterns with group comparisons.
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_point labs scale_color_manual scale_fill_manual theme_bw theme element_blank element_text position_jitterdodge
#' @importFrom ggpubr stat_compare_means
#' @importFrom dplyr pull
#' @keywords visualization expression pancancer
singlePancancer = function(data, gene, paired,color = c("#E7B800", "#2E9FDF")) {
  ldat <- data[data$gene == gene, ]
  ldat$cancer = gsub("TCGA-", "", ldat$cancer)
  # 核心：使用分位数计算y轴范围，排除极端离群值
  y_limits <- ldat %>%
    pull(exp) %>% # 提取y值列
    quantile(probs = c(0.02, 0.98), na.rm = TRUE) # 取2%和98%分位数，过滤极端值

  # 计算适当的缓冲（基于分位数范围的5%）
  buffer <- diff(y_limits) * 0.1
  y_lower <- y_limits[1] - buffer # 下限
  y_upper <- y_limits[2] + buffer # 上限

  if(paired){ldat$group = ldat$Sample}
  ggplot(ldat, aes(x = cancer, y = exp, color = group, fill = group)) +
    # 先绘制箱线图
    geom_boxplot(
      aes(color = group),
      alpha = 0.6,
      lwd = 0.1,
      outlier.size = 1,
      outlier.colour = "white"
    ) +
    # 然后绘制散点图，使用position_jitterdodge使其与箱线图对齐
    geom_point(
      aes(color = group), # 确保散点颜色与group映射
      size = 0.5, # 散点大小
      alpha = 0.8, # 散点透明度
      position = position_jitterdodge(
        jitter.width = 0.2, # 水平方向抖动幅度
        jitter.height = 0, # 垂直方向不抖动
        dodge.width = 0.75 # 与箱线图的dodge宽度保持一致
      )
    ) +
    theme_bw() +
    labs(y = paste0("The expression of ", gene, "\nlog2(TPM + 1)")) +
    # 添加统计检验（组间比较）
    ggpubr::stat_compare_means(
      label = "p.signif", # 显示显著性标记（*/*/***）
      method = "wilcox.test", # 箱线图常用非参数检验，可根据数据类型调整
      vjust = 0.5,
      hjust = 0.5,
      paired = paired,
      hide.ns = FALSE, # 显示不显著(ns)的标注，如需隐藏可设为TRUE
      show.legend = FALSE # 添加这一行，不显示统计检验的图例
    ) +
    # ylim(c(y_lower, y_upper)) +
    # 颜色设置 - 修改这里来设置不同组的颜色
    scale_color_manual(values = color) + # 替换为你想要的颜色
    scale_fill_manual(values = color) +
    # 主题设置
    theme_bw() + # 基于白色背景的主题，便于去除网格线
    theme(
      # 移除网格线
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      # 图例位置在顶部
      legend.position = "top",
      legend.direction = "horizontal",  # 水平排列
      # 图例标题不显示，文本样式调整
      legend.title = element_blank(),
      legend.background = element_rect(fill = "transparent"),
      legend.text = element_text(size = 10,  colour = "#1A1A1A"),
      # 坐标轴文本样式
      axis.text.x = element_text(
        size = 10,
        angle = 45,
        hjust = 1,
        colour = "#1A1A1A"
      ),
      axis.text.y = element_text(size = 10, colour = "#1A1A1A"),
      # 坐标轴标题样式（x轴标题不显示）
      axis.title.y = element_text(size = 10, colour = "#1A1A1A"),
      axis.title.x = element_blank(),
      strip.background = element_rect(fill = "white", color = NA) # color=NA去除边框
    )
}


#' Plot Paired Gene Expression with Connecting Lines
#'
#' Generate a point plot with connecting lines for paired samples (e.g., matched normal-tumor) to visualize expression changes of a single gene in a specific cancer type, with paired t-test results.
#'
#' @param data A data frame containing paired gene expression data. Must include columns: `gene` (gene names), `cancer` (cancer type), `Sample` (two groups, e.g., "Normal" and "Tumor"), `exp` (expression values, log2(TPM + 1)), and `id` (unique identifier for paired samples).
#' @param cancer A character string specifying the target cancer type to plot.
#' @param gene A character string specifying the target gene to plot.
#' @param color A character vector of length 2 specifying point colors for the two sample groups. Defaults to `c("#E7B800", "#2E9FDF")`.
#' @param linecolor A character string specifying the color of connecting lines between paired samples. Defaults to `"#1A1A1A"`.
#'
#' @return A ggplot object showing paired expression values with connecting lines and statistical annotations.
#' @importFrom ggplot2 ggplot aes geom_line geom_point labs scale_x_continuous scale_color_manual theme_classic theme element_rect element_text annotate
#' @importFrom dplyr group_by mutate ungroup
#' @importFrom reshape2 dcast
#' @keywords visualization expression paired
singleGenePointLinePolt = function(
    data,
    cancer,
    gene,
    color = c("#E7B800", "#2E9FDF"),
    linecolor = "#1A1A1A"
) {
  # 正确的数据筛选
  ldat <- data[data$gene == gene, ]
  ldat <- ldat[ldat$cancer == cancer, ]

  # 检查数据是否包含配对信息
  if (length(unique(ldat$Sample)) != 2) {
    stop("数据中必须包含两个样本组（Normal和Tumor）")
  }

  # 更安全的配对t检验数据准备
  dat1_paired <- reshape2::dcast(
    ldat,
    id ~ Sample, # id作为行，Sample作为列
    value.var = "exp"
  )

  # 检查数据完整性
  if (!all(c("Normal", "Tumor") %in% colnames(dat1_paired))) {
    stop("数据中必须包含Normal和Tumor两组")
  }

  # 移除缺失配对的数据
  dat1_paired <- na.omit(dat1_paired)

  # 配对t检验
  if (nrow(dat1_paired) > 1) {
    p.value <- t.test(
      dat1_paired$Normal,
      dat1_paired$Tumor,
      paired = TRUE,
      alternative = 'two.sided',
      conf.level = 0.95
    )
    pv <- p.value$p.value
  } else {
    pv <- NA
    warning("配对数据不足，无法进行t检验")
  }

  # 格式化p值显示
  stasig <- ifelse(
    is.na(pv),
    "Insufficient data",
    ifelse(
      pv >= 0.001,
      paste0('p = ', round(pv, 3)),
      "p < 0.001"
    )
  )

  # 设置随机种子保证可重复性
  set.seed(123)

  # 为每个样本对计算相同的抖动值
  ldat_jitter <- ldat %>%
    group_by(id) %>%
    mutate(
      jitter_amount = runif(1, -0.1, 0.1) # 为每个id生成唯一的抖动值
    ) %>%
    ungroup() %>%
    mutate(
      x_jitter = as.numeric(factor(Sample)) + jitter_amount # 应用相同的抖动
    )
  # 确定p值显示的坐标位置
  normal_mean <- mean(ldat$exp[ldat$Sample == "Normal"], na.rm = TRUE)
  tumor_mean <- mean(ldat$exp[ldat$Sample == "Tumor"], na.rm = TRUE)
  if(pv < 0.05){
    nmin = min(ldat$exp[ldat$Sample == "Normal"])
    tmin = min(ldat$exp[ldat$Sample == "Tumor"])
    y.pos = max(nmin,tmin) - abs((nmin - tmin))/2
    x.pos = ifelse(normal_mean > tumor_mean,1,2)
  }else{
    x.pos = 1.5
    y.pos = max(ldat$exp)
  }
  # 使用手动计算的抖动位置绘图
  p <-  ggplot(
    ldat_jitter,
    aes(x = x_jitter, y = exp, colour = Sample)
  ) +
    # 连线使用抖动后的位置
    geom_line(aes(group = id), color = linecolor, linewidth = 0.3, alpha = 0.3) +
    # 点使用相同的抖动位置
    geom_point(size = 2, alpha = 0.6) +
    # 手动设置x轴标签和位置 - 调整expand参数减少两侧空白
    scale_x_continuous(
      breaks = 1:2,
      labels = levels(factor(ldat$Sample)),
      limits = c(0.5, 2.5),
      expand = expansion(mult = 0.05)  # 减少这个值来缩小两侧空白
    ) +
    scale_color_manual(
      values = color
    ) +
    theme_classic() +
    labs(
      y = paste0("Expression of ", gene, "\nlog2(TPM + 1)"),
      title = cancer
    ) +
    theme(
      panel.background = element_rect(
        fill = "white",
        colour = "black",
        linewidth = 0.25
      ),
      plot.title = element_text(hjust = 0.5),
      axis.line = element_line(colour = "black", linewidth = 0.25),
      axis.title.x = element_blank(),
      axis.text.x = element_text(face = "plain", colour = "black"),
      axis.text = element_text(size = 12, face = "plain", color = "black"),
      legend.position = "none"
    ) +
    annotate(
      'text',
      x = x.pos,
      y = y.pos,
      label = stasig,
      size = 4
    )

  return(p)
}
