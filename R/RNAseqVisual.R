#' Plot Volcano Figure for Differentially Expressed Genes
#'
#' Generates a volcano plot to visualize the relationship between log2 fold change and statistical significance of differentially expressed genes, with optional gene labeling.
#'
#' @param data A data frame containing differential expression results. Must include columns for log2 fold change, statistical significance, group classification, and gene labels.
#' @param x Character string specifying the column name in `data` representing log2 fold change (e.g., "log2FC").
#' @param y Character string specifying the column name in `data` representing statistical significance. Must be one of "pValue", "FDR", or "pAdj".
#' @param cut_pvalue Numeric threshold for statistical significance (e.g., 0.05) used to draw the horizontal dashed line.
#' @param cutFC Numeric threshold for log2 fold change (e.g., 1) used to draw vertical dashed lines.
#' @param title Character string specifying the plot title.
#' @param group Character string specifying the column name in `data` representing gene groups (e.g., "direction" with values "Up", "Down", "Ns").
#' @param label Character string specifying the column name in `data` containing gene labels. Empty strings indicate genes not to be labeled.
#' @param colour A character vector of length 3 specifying colors for up-regulated, down-regulated, and non-significant genes, respectively. Defaults to c("#E7B800", "#2E9FDF", "#58ae9a").
#'
#' @return A ggplot2 object (volcano plot)
#' @export plotDEGvolcanoFig
#' @import ggplot2
#' @import ggrepel
#'
plotDEGvolcanoFig <- function(data, x, y, cut_pvalue, cutFC, title, group, label,
                              colour = c("#E7B800", "#2E9FDF", "#58ae9a")) {
  # Input validation
  required_cols <- c(x, y, group, label)
  missing_cols <- setdiff(required_cols, colnames(data))
  if (length(missing_cols) > 0) {
    stop("Data frame is missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  if (!y %in% c("pValue", "FDR", "pAdj")) {
    stop("'y' must be one of: 'pValue', 'FDR', 'pAdj'")
  }

  if (length(colour) != 3) {
    stop("'colour' must be a vector of length 3 (for up, down, and non-significant genes)")
  }

  # Data preprocessing
  plot_data <- data[, required_cols]
  colnames(plot_data) <- c("log2FC", "sig", "group", "label")

  # Volcano plot construction (使用segment.size替代segment.linewidth)
  p <- ggplot(plot_data, aes(x = log2FC, y = -log10(sig), colour = group)) +
    geom_point(aes(size = ifelse(label != "", 1.2, 1)), alpha = 0.8) +
    scale_color_manual(values = colour) +
    scale_size_identity() +
    geom_vline(xintercept = c(-cutFC, cutFC), linetype = "dashed",
               color = "#666666", linewidth = 0.5) +
    geom_hline(yintercept = -log10(cut_pvalue), linetype = "dashed",
               color = "#666666", linewidth = 0.5) +
    ggrepel::geom_text_repel(
      data = plot_data[plot_data$label != "", ],
      aes(label = label),
      size = 4,
      segment.color = "#999999",
      segment.size = 0.3,  # 替换为旧版本兼容的参数
      box.padding = 0.5,
      show.legend = FALSE
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
      axis.title = element_text(size = 14, color = "#333333"),
      axis.text = element_text(size = 12, color = "#333333"),
      legend.title = element_blank(),
      legend.position = "top",
      legend.key = element_blank(),
      legend.text = element_text(size = 12),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "#666666", linewidth = 0.6)
    ) +
    labs(
      x = expression(log[2]("Fold Change")),
      y = switch(y,
                 pValue = expression(-log[10]("p Value")),
                 FDR = expression(-log[10]("FDR")),
                 pAdj = expression(-log[10]("Adjusted p Value"))
      ),
      title = title
    )

  return(p)
}
