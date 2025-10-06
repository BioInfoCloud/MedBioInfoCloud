
#' Preprocess enrichment results for visualization
#'
#' Filters significant enrichment entries and extracts key information for plotting.
#' Works with both GO and KEGG enrichment results from clusterProfiler.
#'
#' @param obj An object of class 'enrichResult' (from clusterProfiler package)
#' @param type Type of enrichment analysis, must be either "GO" or "KEGG"
#' @param nTerm Maximum number of terms to show per group (for GO: BP/MF/CC), default 6
#' @param p.adjust Significance threshold for filtering, default 0.05
#'
#' @return Data frame with columns:
#'   - group: Functional category (BP/MF/CC for GO, "KEGG" for KEGG)
#'   - term: Description of the enrichment term
#'   - count: Number of genes in the term
#' @export
#' @importFrom dplyr group_by arrange slice_head ungroup filter select mutate
#' @importFrom stringr str_to_title
preVisEnrishResults <- function(obj,type,nTerm = 6,p.adjust = 0.05){
  if (!inherits(obj, "enrichResult") || nrow(obj) == 0) {
    stop("obj is not a valid enrichResult object or contains no results")
  }
  # filter.result = obj@result %>% filter(p.adjust < p.adjust)
  filter.result = obj@result
  filter.result = filter.result[filter.result$p.adjust < p.adjust,]

  if(nrow(filter.result)==0){
    stop("没有显著性的条目被富集，可将p.adjust的值设置大一些")
  }
  # 加载dplyr包用于数据处理
  library(dplyr)
  if(type == "GO"){
    data = filter.result[,c(1,3,13)]
    # 提取每个ONTOLOGY类别中Count前10的条目
    top_data <- data %>%
      # 按ONTOLOGY分组
      group_by(ONTOLOGY) %>%
      # 每组按Count降序排序（取Count最大的在前）
      arrange(desc(Count), .by_group = TRUE) %>%
      # 每组取前nTerm行，不足nTerm行则取全部
      slice_head(n = nTerm) %>%
      # 取消分组
      ungroup()
    top_data$Description = str_to_title(top_data$Description)
  }else if(type == "KEGG"){
    top_data = filter.result[1:min(nTerm,nrow(filter.result)),c(2,4,14)]
  }
  return(top_data)

}

#' Prepare data for circular bar chart visualization
#'
#' Processes enrichment results to create formatted data structures required for
#' generating polar bar charts, including IDs, label angles, and baseline positions.
#'
#' @param data Data frame with columns: group, term, count (output from preVisEnrishResults)
#'
#' @return List containing:
#'   - data: Original data with added ID column
#'   - label_data: Data with label angles and justification for text placement
#'   - base_data: Data for drawing group baseline segments
#'   - grid_data: Data for drawing circular grid lines
#' @export
#' @importFrom dplyr group_by summarize rowwise mutate arrange
preCirBarchartData <- function(data){
  colnames(data) = c("group","term","count")
  # Set a number of 'empty bar' to add at the end of each group
  empty_bar <- 0
  to_add <- data.frame( matrix(NA, empty_bar*nlevels(data$group), ncol(data)) )
  colnames(to_add) <- colnames(data)
  to_add$group <- rep(levels(data$group), each=empty_bar)
  data <- rbind(data, to_add)
  data <- data %>% arrange(group)
  data$id <- seq(1, nrow(data))

  # Get the name and the y position of each label
  label_data <- data
  number_of_bar <- nrow(label_data)
  angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
  label_data$hjust <- ifelse( angle < -90, 1, 0)
  label_data$angle <- ifelse(angle < -90, angle+180, angle)

  # prepare a data frame for base lines
  base_data <- data %>%
    group_by(group) %>%
    summarize(start=min(id), end=max(id) - empty_bar) %>%
    rowwise() %>%
    mutate(title=mean(c(start, end)))

  # prepare a data frame for grid (scales)
  grid_data <- base_data
  grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
  grid_data$start <- grid_data$start - 1
  grid_data <- grid_data[-1,]
  return(list(data = data,
              label_data = label_data,
              base_data = base_data,
              grid_data = grid_data))
}


#' Generate circular bar chart for enrichment results
#'
#' Creates a polar bar chart visualization of enrichment results with:
#' - Outer ring: Grouped bars (by GO category/KEGG) with height representing gene count
#' - Inner ring: Uniform height bars with color intensity showing gene count
#' - Value labels: Gene counts displayed in the inner ring
#'
#' @param data Data frame with columns: group, term, count (output from preVisEnrishResults).
#' @param fillColors Vector of colors for group distinction, default provides 5 colors.
#' @param ylim Range of coordinate axis y, a numerical vector of length 2.
#'
#' @return ggplot2 object with circular bar chart visualization.
#' @export
#' @import ggplot2
#' @importFrom RColorBrewer brewer.pal
enrichCirBarchar <- function(data,fillColors = c("#E7B800", "#2E9FDF","#58ae9a","#646e9a","#efa78e","#Afa78e","#08ae9a","#B7B800"),
                             ylim = c(-150,150)){
  predat  = preCirBarchartData(data)
  data = predat$data
  label_data = predat$label_data
  base_data = predat$base_data
  grid_data = predat$grid_data

  groupnum = length(unique(data$group))
  grouphjust = c(rep(1,groupnum%/%2),rep(0,groupnum - groupnum%/%2))
  # Make the plot
  p <- ggplot(data, aes(x=as.factor(id), y= count, fill=group)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar

    geom_bar(aes(x=as.factor(id), y=count, fill=group), stat="identity", alpha=0.1) +

    # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
    # geom_segment(data=grid_data, aes(x = end, y = 80, xend = start, yend = 80), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    # geom_segment(data=grid_data, aes(x = end, y = 60, xend = start, yend = 60), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    # geom_segment(data=grid_data, aes(x = end, y = 40, xend = start, yend = 40), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    # geom_segment(data=grid_data, aes(x = end, y = 20, xend = start, yend = 20), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +

    # Add text showing the value of each 100/75/50/25 lines
    # annotate("text", x = rep(max(data$id),4), y = c(20, 40, 60, 80),
    #          label = c("20", "40", "60", "80") ,
    #          color="grey", size=3 , angle=0, fontface="bold", hjust=-1) +
    scale_fill_manual(
      values = fillColors,  # 选项1：使用自定义颜色向量
      # values = brewer_colors  # 选项2：使用RColorBrewer预定义调色板
    ) +

    geom_bar(aes(x=as.factor(id), y=count, fill=group), stat="identity", alpha=0.5) +
    ylim(ylim[1],(min(max(data$count),ylim[2]) + 5)) +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.margin = unit(rep(-1,4), "cm")
    ) +
    coord_polar() +
    geom_text(data=label_data, aes(x=id, y=0, label= term, hjust=hjust), color="black",
              fontface="bold",alpha=0.6, size=2.5, angle= label_data$angle, inherit.aes = FALSE ) +

    # Add base line information
    geom_segment(data=base_data, aes(x = start, y = -5, xend = end, yend = -5),
                 colour = "black",
                 alpha=0.8, linewidth=0.6 , inherit.aes = FALSE )  +
    geom_text(data=base_data, aes(x = title, y = -18, label=group),
              hjust= grouphjust, colour = "black", alpha=0.8, size=4,
              fontface="bold", inherit.aes = FALSE)
  return(p)
}
