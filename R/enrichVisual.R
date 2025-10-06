
#' Create a common bar plot for enrichment analysis results
#'
#' Generates a horizontal bar plot to visualize enrichment analysis results,
#' with bars colored by statistical significance.
#'
#' @param data A data frame containing enrichment results. Must include columns:
#'   - 'Count': Number of genes in each term
#'   - 'Term': Enrichment terms (y-axis categories)
#'   - 'pvalue': P-values for coloring
#' @param fill_by A character string specifying the column name in `data` used for coloring bars.
#'   Typically "pvalue", "p.adjust", or "-log10(pvalue)". Defaults to "-log10(pvalue)".
#' @param color_palette A character string specifying the color palette name (from RColorBrewer)
#'   or a vector of custom colors for gradient filling. Defaults to "RdPu".
#' @param color_direction Numeric (1 or -1) indicating the direction of the color gradient.
#'   Use 1 for original order, -1 to reverse. Defaults to 1.
#' @param axisTitle.x Character string for x-axis title.
#' @param axisTitle.y Character string for y-axis title.
#' @param title Character string for plot title.
#'
#' @return A ggplot2 object
#' @export common_bar
#' @import ggplot2
#' @import RColorBrewer
common_bar <- function(data,
                       fill_by = "-log10(pvalue)",
                       color_palette = "RdPu",
                       color_direction = 1,
                       axisTitle.x,
                       axisTitle.y,
                       title) {
  # Validate required columns
  required_cols <- c("Count", "Term", ifelse(fill_by == "-log10(pvalue)", "pvalue", fill_by))
  missing_cols <- setdiff(required_cols, colnames(data))
  if (length(missing_cols) > 0) {
    stop("Data is missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # Calculate -log10(pvalue) if needed
  if (fill_by == "-log10(pvalue)" && !"-log10(pvalue)" %in% colnames(data)) {
    data <- data %>% mutate(`-log10(pvalue)` = -log10(pvalue))
  }

  # Create base plot
  ggplot(data, aes(x = Count, y = Term, fill = .data[[fill_by]])) +
    geom_bar(stat = "identity", width = 0.8) +
    scale_fill_distiller(palette = color_palette, direction = color_direction) +
    labs(x = axisTitle.x,
         y = axisTitle.y,
         title = title,
         fill = ifelse(fill_by == "-log10(pvalue)", "-log10(pvalue)", fill_by)) +
    theme_classic() +
    theme(
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 11),
      plot.title = element_text(size = 14, hjust = 0.5),
      legend.title = element_text(size = 13),
      legend.text = element_text(size = 11)
    )
}


#' Visualize enrichment results as bar plot
#'
#' Generates bar plots for GO/KEGG enrichment results from clusterProfiler's enrichResult object.
#' Supports multiple styles and optional saving to file.
#'
#' @param obj An enrichResult object from clusterProfiler package.
#' @param showCategory Integer specifying the number of top enrichment terms to display.
#'   Defaults to 6.
#' @param fill_by Column name used for coloring bars. Options include "pvalue", "p.adjust",
#'   or "-log10(pvalue)". Defaults to "-log10(pvalue)".
#' @param color_palette Color palette name (from RColorBrewer) or custom color vector.
#'   Defaults to "RdPu".
#' @param color_direction Direction of color gradient (1 or -1). Defaults to 1.
#' @param axisTitle.x X-axis title. Defaults to "Number of Genes".
#' @param axisTitle.y Y-axis title. Defaults to "Term Name".
#' @param title Plot title. Defaults to "Enrichment Barplot".
#' @param save Logical indicating whether to save the plot. Defaults to FALSE.
#' @param folder Directory path for saving the plot. Defaults to "./".
#' @param fileName File name prefix for saved plot. Defaults to "EnrichBar".
#' @param height Plot height in inches. Defaults to 6.
#' @param width Plot width in inches. Defaults to 10.
#' @param style Plot style. Options: "style1" (standard bar plot) or "style2" (with term labels
#'   inside bars). Defaults to "style1".
#'
#' @return A ggplot2 object if plot is generated; NULL otherwise.
#' @export GO_KEGG.enrichVisual_barplot
#' @import ggplot2
#' @import dplyr
#' @import stringr
GO_KEGG.enrichVisual_barplot <- function(obj,
                                         showCategory = 6,
                                         fill_by = "-log10(pvalue)",
                                         color_palette = "RdPu",
                                         color_direction = 1,
                                         axisTitle.x = "Number of Genes",
                                         axisTitle.y = "Term Name",
                                         title = "Enrichment Barplot",
                                         save = FALSE,
                                         folder = "./",
                                         fileName = "EnrichBar",
                                         height = 6,
                                         width = 10,
                                         style = "style1") {
  # Validate input object
  if (!inherits(obj, "enrichResult") || nrow(obj) == 0) {
    warning("Input is not a valid enrichResult object or contains no results")
    return(NULL)
  }

  # Limit to number of categories to show
  showCategory <- min(showCategory, nrow(obj))
  if (showCategory <= 0) {
    warning("No categories to display (showCategory <= 0)")
    return(NULL)
  }

  # Prepare data
  enrich_result <- obj %>%
    as.data.frame() %>%
    slice_head(n = showCategory) %>%
    arrange(desc(p.adjust)) %>%
    mutate(
      Term = factor(
        stringr::str_to_title(Description),
        levels = rev(stringr::str_to_title(Description))
      )
    )

  # Generate plot based on style
  p <- switch(style,
              "style1" = common_bar(
                data = enrich_result,
                fill_by = fill_by,
                color_palette = color_palette,
                color_direction = color_direction,
                axisTitle.x = axisTitle.x,
                axisTitle.y = axisTitle.y,
                title = title
              ),

              "style2" = common_bar(
                data = enrich_result,
                fill_by = fill_by,
                color_palette = color_palette,
                color_direction = color_direction,
                axisTitle.x = axisTitle.x,
                axisTitle.y = axisTitle.y,
                title = title
              ) +
                geom_bar(stat = "identity", width = 0.8, alpha = 0.7) +
                geom_text(aes(x = 0.03, label = Term), hjust = 0, size = 3.5) +
                theme(
                  axis.text.y = element_blank(),
                  axis.ticks.y = element_blank()
                ),

              stop("Invalid style. Choose 'style1' or 'style2'")
  )

  # Save plot if required
  if (save) {
    if (!dir.exists(folder)) {
      dir.create(folder, recursive = TRUE, showWarnings = FALSE)
    }
    file_path <- file.path(folder, paste0(fileName, "-", style, ".pdf"))
    ggsave(
      filename = file_path,
      plot = p,
      height = height,
      width = width,
      device = "pdf"
    )
    message("Plot saved to: ", file_path)
  }

  return(p)
}

