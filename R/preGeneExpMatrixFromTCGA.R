#' Extract gene expression matrix from TCGA data with gene type and list filtering
#' 
#' This function loads TCGA gene expression data (TPM or count), filters genes by type (e.g., protein-coding) and/or a specified gene list,
#' and returns the processed matrix with control over retaining or removing the specified genes.
#' 
#' @param dataFolder Character string. Path to the folder containing TCGA expression data files.
#' @param cancer Character string. TCGA cancer type abbreviation (e.g., "LUAD", "BRCA").
#' @param type Character string. Type of expression data to extract: "tpm" (transcripts per million) or "count" (read counts). Defaults to "tpm".
#' @param genes Character vector. Genes of interest to filter. If NULL (default), no gene list filtering is applied (returns all genes passing type filter).
#' @param keep_genes Logical. Only effective when `genes` is not NULL. If TRUE, retain genes in `genes`; if FALSE, remove genes in `genes` and keep others. Defaults to TRUE.
#' @param filter Gene type(s) to retain. Can be: 
#'   - NULL (default): retain all gene types (no filtering by type).
#'   - Character string or vector: retain genes where `gene_type` matches these values (e.g., "protein_coding", c("lncRNA", "miRNA")).
#'   Common values include "protein_coding", "lncRNA", "snRNA", "miRNA", "misc_RNA" (see details for full list).
#' @param pattern Character string. Regular expression pattern to match data files (default: "STARdata.Rdata$").
#' 
#' @return A matrix with genes as rows and samples as columns, containing the filtered expression values.
#' @details 
#' - Gene type filtering (via `filter`) is applied first: only genes with `gene_type` in `filter` (or all genes if `filter = FALSE`) are retained.
#' - Gene list filtering (via `genes` and `keep_genes`) is applied next: 
#'   - If `genes = NULL`, returns all genes passing the type filter.
#'   - If `genes` is provided, retains/removes these genes from the type-filtered results (based on `keep_genes`).
#' - Common gene types include: "protein_coding", "lncRNA", "snRNA", "miRNA", "misc_RNA", "snoRNA", "rRNA", etc.
#' @export
preGeneExpMatrixFromTCGA <- function(
  dataFolder,
  cancer,
  type = "tpm",
  genes = NULL,
  keep_genes = TRUE,
  filter =  NULL,
  pattern = "STARdata.Rdata$"
) {
  # Get full paths of files matching the pattern
  file_paths <- dir(dataFolder, pattern = pattern, full.names = TRUE)
  if (length(file_paths) == 0) {
    stop("No files matching the pattern found in 'dataFolder'.")
  }
  
  # Load the file corresponding to the specified cancer type
  target_file <- file_paths[grep(cancer, file_paths)]
  if (length(target_file) == 0) {
    stop(paste("No file found for cancer type:", cancer))
  }
  load(target_file)  # Loads object 'STARdata'
  
  # Validate 'type' parameter
  if (!type %in% c("tpm", "count")) {
    stop("'type' must be either 'tpm' or 'count'.")
  }
  expr_matrix <- STARdata[[type]]
  if (is.null(expr_matrix)) {
    stop(paste("Expression data of type '", type, "' not found in the file.", sep = ""))
  }
  
  # Filter by gene type using filterGeneTypeExpr
  # - If filter =  NULL: retain all gene types
  # - If filter is character(s): retain only those gene types
  expr_matrix <- filterGeneTypeExpr(expr = expr_matrix, fil_col = "gene_type", filter = filter)
  all_genes_in_data <- rownames(expr_matrix)
  
  # Check if any genes remain after type filtering
  if (nrow(expr_matrix) == 0) {
    stop("No genes retained after gene type filtering. Check 'filter' parameter (valid gene types may include 'protein_coding', 'lncRNA', etc.).")
  }
  
  # Determine genes to retain based on 'genes' and 'keep_genes' (applied after type filtering)
  retained_genes <- if (is.null(genes)) {
    # When genes is NULL: return all genes passing type filter, ignore keep_genes
    if (!keep_genes) {
      warning("'keep_genes' is ignored because 'genes' is NULL (no target genes to remove). Returning all genes passing type filter.")
    }
    all_genes_in_data
  } else {
    # When genes is provided: use keep_genes to decide retain/remove from type-filtered genes
    if (keep_genes) {
      intersect(genes, all_genes_in_data)  # Keep genes in both 'genes' and type-filtered data
    } else {
      setdiff(all_genes_in_data, genes)  # Remove 'genes' from type-filtered data
    }
  }
  
  # Validate retained genes after all filtering
  if (length(retained_genes) == 0) {
    stop("No genes retained after all filtering. Check 'genes', 'keep_genes', and 'filter' parameters.")
  }
  
  # Generate status messages for gene type filtering
  if (isFALSE(filter)) {
    message("No gene type filtering applied (retained all gene types).")
  } else {
    filter_types <- if (is.character(filter)) paste(filter, collapse = ", ") else as.character(filter)
    retained_type_count <- nrow(expr_matrix)
    message(paste0("Retained ", retained_type_count, " genes of type(s): ", filter_types))
  }
  
  # Generate status messages for gene list filtering
  if (is.null(genes)) {
    message(paste0("Returning all ", length(retained_genes), " genes (since 'genes' is NULL)."))
  } else {
    if (keep_genes) {
      mapped <- length(retained_genes)
      total_input <- length(genes)
      mapped_pct <- round(mapped * 100 / total_input, 2)
      message(paste0("Retained ", mapped, " (", mapped_pct, "%) of input genes (matched in type-filtered data)."))
      if (mapped < total_input) {
        unmapped <- setdiff(genes, all_genes_in_data)
        message(paste0("Unmatched genes (not in type-filtered data): ", paste(unmapped, collapse = ", ")))
      }
    } else {
      removed <- length(intersect(genes, all_genes_in_data))
      total_removed <- length(genes)
      removed_pct <- if (total_removed > 0) round(removed * 100 / total_removed, 2) else 0
      kept <- length(retained_genes)
      message(paste0("Removed ", removed, " (", removed_pct, "%) of input genes from type-filtered data."))
      message(paste0("Retained ", kept, " genes (not in input 'genes')."))
    }
  }
  
  # Return final filtered matrix (preserve matrix structure)
  return(expr_matrix[retained_genes, , drop = FALSE])
}
