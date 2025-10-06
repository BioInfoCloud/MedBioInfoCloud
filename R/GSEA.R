


#' Prepare Pre-Ranked Gene List for GSEA
#'
#' Converts input data (either a data frame or a named numeric vector) into a pre-ranked gene list
#' required for Gene Set Enrichment Analysis (GSEA). The output is a named numeric vector where:
#' - Names are gene identifiers (e.g., symbols);
#' - Values are ranking metrics (e.g., log2 fold change, t-statistic), sorted in descending order.
#'
#' @param data Input data for GSEA. Two valid formats:
#'   1. A data frame containing at least:
#'      - A column with gene identifiers (specified by `geneColName`);
#'      - A column with ranking values (specified by `rangeColName`).
#'   2. A pre-ranked numeric vector where names are gene identifiers and values are ranking metrics.
#' @param geneColName Character string. Column name in `data` (if data frame) corresponding to gene identifiers.
#'   Required only if `data` is a data frame; ignored if `data` is a numeric vector.
#' @param rangeColName Character string. Column name in `data` (if data frame) corresponding to ranking values.
#'   Required only if `data` is a data frame; ignored if `data` is a numeric vector.
#'
#' @return A named numeric vector (pre-ranked in descending order) where:
#'   - Names = gene identifiers (from `geneColName` if data frame, or original names if vector);
#'   - Values = ranking metrics (sorted in descending order).
#' @export
#' @importFrom dplyr select rename filter arrange distinct
#'
#' @examples
#' # Example 1: Data frame input
#' de_data <- data.frame(
#'   Symbol = c("TP53", "MYC", "EGFR", "TP53", "RB1"),
#'   Log2FC = c(-1.8, 2.5, 1.2, -0.5, -1.5),
#'   stringsAsFactors = FALSE
#' )
#' ranked_list1 <- preRankedGeneList(
#'   data = de_data,
#'   geneColName = "Symbol",
#'   rangeColName = "Log2FC"
#' )
#'
#' # Example 2: Named vector input
#' raw_vector <- c(-1.8, 2.5, 1.2, -1.5)
#' names(raw_vector) <- c("TP53", "MYC", "EGFR", "RB1")
#' ranked_list2 <- preRankedGeneList(data = raw_vector)
preRankedGeneList <- function(data,
                              geneColName = NULL,
                              rangeColName = NULL) {
  if (is.data.frame(data)) {
    # Validate required parameters for data frame input
    if (is.null(geneColName) || is.null(rangeColName)) {
      stop("For data frame input, both 'geneColName' and 'rangeColName' must be specified.", call. = FALSE)
    }
    if (!all(c(geneColName, rangeColName) %in% colnames(data))) {
      missing_cols <- setdiff(c(geneColName, rangeColName), colnames(data))
      stop(sprintf(
        "Data frame is missing required columns:\n%s",
        paste(missing_cols, collapse = ", ")
      ), call. = FALSE)
    }

    # Process data frame: filter, deduplicate, and sort
    ranked_data <- data %>%
      dplyr::select(!!geneColName, !!rangeColName) %>%
      dplyr::rename(gene = !!geneColName, rank_value = !!rangeColName) %>%  # Generic column names
      dplyr::filter(!is.na(gene), !is.na(rank_value)) %>%  # Remove missing values
      dplyr::arrange(desc(abs(rank_value))) %>%  # Prioritize largest absolute values for deduplication
      dplyr::distinct(gene, .keep_all = TRUE) %>%  # Remove duplicate genes
      dplyr::arrange(desc(rank_value))  # Sort by ranking value (descending)

    # Convert to named vector
    gene_list <- ranked_data$rank_value
    names(gene_list) <- ranked_data$gene

  } else if (is.numeric(data) && is.vector(data) && !is.null(names(data))) {
    # Process numeric vector input
    # Remove NA values and duplicate genes
    gene_list <- data[!is.na(data)]
    gene_list <- gene_list[!duplicated(names(gene_list))]
    # Sort in descending order
    gene_list <- sort(gene_list, decreasing = TRUE)

  } else {
    # Invalid input format
    stop(
      "Invalid 'data' format. Supported formats:\n",
      "1. Data frame with 'geneColName' and 'rangeColName'\n",
      "2. Named numeric vector (names = gene identifiers)",
      call. = FALSE
    )
  }

  # Validate output (non-empty)
  if (length(gene_list) == 0) {
    stop("No valid genes remaining after filtering (check for missing values or duplicates).", call. = FALSE)
  }

  return(gene_list)
}

#' Prepare TERM2GENE for GSEA based on MSigDB
#' @param species Character string. Species name for MSigDB gene sets (e.g., "Homo sapiens" for human, "Mus musculus" for mouse).
#'   Use `msigdbr::msigdbr_species()` to view all available species. Defaults to "Homo sapiens".
#' @param collection Character string. MSigDB collection abbreviation (e.g., "C5" for GO gene sets, "H" for hallmark gene sets).
#'   Use `msigdbr::msigdbr_collections()` to view all available collections. Defaults to "C5".
#' @param subcollection Character string. MSigDB subcollection abbreviation (e.g., "GO:BP" for GO Biological Process, "CGP" for chemical and genetic perturbations).
#'   Must match the specified `collection`; use `msigdbr::msigdbr_collections()` to view valid subcollections for each collection. Defaults to "GO:BP".
preTERM2GENEfromMSigDB <- function(species = "Homo sapiens",
                                   collection = "C5",
                                   subcollection = "GO:BP"){
  # -------------------------- 2. Validate MSigDB Parameters --------------------------
  # Validate species (must be in msigdbr's available species)
  valid_species <- msigdbr::msigdbr_species()$species_name
  if (!species %in% valid_species) {
    stop(sprintf(
      "Invalid 'species': '%s' is not supported.\nAvailable species:\n%s",
      species,
      paste(valid_species, collapse = ", ")
    ), call. = FALSE)
  }

  # Validate collection (must be in msigdbr's available collections)
  msig_collections <- msigdbr::msigdbr_collections()
  valid_collections <- unique(msig_collections$gs_collection)
  if (!collection %in% valid_collections) {
    stop(sprintf(
      "Invalid 'collection': '%s' is not supported.\nAvailable collections:\n%s",
      collection,
      paste(valid_collections, collapse = ", ")
    ), call. = FALSE)
  }

  # Validate subcollection (must match the specified collection)
  valid_subcols <- unique(msig_collections$gs_subcollection[msig_collections$gs_collection == collection])
  # Allow NULL subcollection for collections with no subcols (e.g., Hallmark "H")
  if (!is.null(subcollection) && !subcollection %in% valid_subcols) {
    stop(sprintf(
      "Invalid 'subcollection': '%s' is not available for collection '%s'.\nValid subcollections for '%s':\n%s",
      subcollection, collection, collection,
      paste(valid_subcols, collapse = ", ")
    ), call. = FALSE)
  }

  # -------------------------- 3. Prepare MSigDB Gene Set (TERM2GENE) --------------------------
  # Fetch MSigDB gene set (term-gene mapping)
  msig_gene_set <- msigdbr::msigdbr(
    species = species,
    collection = collection,
    subcollection = subcollection
  )

  # Check if gene set is empty (e.g., invalid subcollection for small species)
  if (nrow(msig_gene_set) == 0) {
    warning(sprintf(
      "No gene sets found for species='%s', collection='%s', subcollection='%s'.\nTry a different combination.",
      species, collection, subcollection
    ), call. = FALSE)
    return(NULL)
  }

  # Format term-gene mapping (required by clusterProfiler::GSEA)
  term2gene <- msig_gene_set %>%
    dplyr::select(gs_name, gene_symbol) %>%
    dplyr::rename(term = gs_name, gene = gene_symbol) %>%
    # Remove duplicate term-gene pairs (if any)
    unique()
  return(term2gene)
}


#' Perform GSEA Analysis Using MSigDB Gene Sets
#'
#' Conduct Gene Set Enrichment Analysis (GSEA) based on pre-ranked gene lists,
#' leveraging the Molecular Signatures Database (MSigDB) for functional gene set annotations.
#' Supports custom gene ranking data (either data frame with gene symbols and fold changes or pre-ranked vector).
#'
#' @param data Input data for GSEA. Two valid formats:
#'   1. A data frame containing at least a gene symbol column and a fold change/ranking column;
#'   2. A pre-ranked numeric vector where names are gene symbols and values are ranking metrics (e.g., log2FC).
#' @param geneColName Character string. Column name in `data` (if data frame) corresponding to gene symbols.
#'   Required only if `data` is a data frame; ignored if `data` is a numeric vector.
#' @param rangeColName Character string. Column name in `data` (if data frame) corresponding to ranking values
#'   (e.g., log2 fold change, t-statistic). Required only if `data` is a data frame; ignored if `data` is a numeric vector.
#' @param species Character string. Species name for MSigDB gene sets (e.g., "Homo sapiens" for human, "Mus musculus" for mouse).
#'   Use `msigdbr::msigdbr_species()` to view all available species. Defaults to "Homo sapiens".
#' @param collection Character string. MSigDB collection abbreviation (e.g., "C5" for GO gene sets, "H" for hallmark gene sets).
#'   Use `msigdbr::msigdbr_collections()` to view all available collections. Defaults to "C5".
#' @param subcollection Character string. MSigDB subcollection abbreviation (e.g., "GO:BP" for GO Biological Process, "CGP" for chemical and genetic perturbations).
#'   Must match the specified `collection`; use `msigdbr::msigdbr_collections()` to view valid subcollections for each collection. Defaults to "GO:BP".
#' @param pvalueCutoff Numeric. Adjusted p-value threshold for filtering significant enrichment results. Defaults to 0.05.
#'
#' @return An object of class `gseaResult` (from clusterProfiler package), containing GSEA results including:
#'   - Enriched gene sets (terms),
#'   - Enrichment scores (ES), normalized enrichment scores (NES),
#'   - Adjusted p-values,
#'   - Genes contributing to each enriched term.
#'   Returns `NULL` if no significant terms are found.
#' @export
#' @importFrom msigdbr msigdbr msigdbr_species msigdbr_collections
#' @importFrom dplyr select rename arrange
#' @importFrom clusterProfiler GSEA
#'
#' @examples
#' # Example 1: Use a data frame with gene symbols and log2FC
#' # Create mock differential expression data
#' de_data <- data.frame(
#'   GeneSymbol = c("TP53", "MYC", "EGFR", "CDKN1A", "RB1"),
#'   Log2FC = c(-1.8, 2.5, 1.2, -0.9, -1.5),
#'   stringsAsFactors = FALSE
#' )
#' # Run GSEA with GO Biological Process gene sets
#' gsea_result1 <- GSEA.baseMSIGDB(
#'   data = de_data,
#'   geneColName = "GeneSymbol",
#'   rangeColName = "Log2FC",
#'   species = "Homo sapiens",
#'   collection = "C5",
#'   subcollection = "GO:BP",
#'   pvalueCutoff = 0.05
#' )
#'
#' # Example 2: Use a pre-ranked numeric vector
#' # Create mock pre-ranked gene list (log2FC as ranking metric)
#' ranked_genes <- c(2.5, 1.2, -0.9, -1.5, -1.8)
#' names(ranked_genes) <- c("MYC", "EGFR", "CDKN1A", "RB1", "TP53")
#' # Run GSEA with Hallmark gene sets
#' gsea_result2 <- GSEA.baseMSIGDB(
#'   data = ranked_genes,
#'   species = "Homo sapiens",
#'   collection = "H",
#'   subcollection = NULL,  # Hallmark has no subcollections
#'   pvalueCutoff = 0.05
#' )
GSEA.baseMSIGDB <- function(data,
                            geneColName = NULL,
                            rangeColName = NULL,
                            species = "Homo sapiens",
                            collection = "C5",
                            subcollection = "GO:BP",
                            pvalueCutoff = 0.05) {
  # -------------------------- 1. Validate Dependencies --------------------------
  # Check if required packages are installed (clusterProfiler is critical for GSEA)
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    stop("Package 'clusterProfiler' is required for GSEA. Install it with:\nBiocManager::install('clusterProfiler')",
         call. = FALSE)
  }


  # -------------------------- 2. Prepare ranked gene list --------------------------
  if (is.numeric(data) && is.vector(data) && !is.null(names(data))){
    geneList = data
  }else {
    geneList = preRankedGeneList(data = data,geneColName = geneColName, rangeColName = rangeColName)
  }
  # -------------------------- 3. Prepare TERM2GENE --------------------------
  term2gene <- preTERM2GENEfromMSigDB (species = species,
                                       collection = collection,
                                       subcollection = subcollection)

  # -------------------------- 4. Run GSEA --------------------------
  gsea_result <- clusterProfiler::GSEA(
    geneList = geneList,
    TERM2GENE = term2gene,
    pvalueCutoff = pvalueCutoff,
    verbose = FALSE  # Disable verbose output to avoid clutter
  )

  # Check if any significant terms are found
  if (nrow(gsea_result@result) == 0) {
    warning(sprintf(
      "No significant GSEA results found with pvalueCutoff = %.2f.\nTry increasing 'pvalueCutoff' or using a different gene set.",
      pvalueCutoff
    ), call. = FALSE)
    return(NULL)
  }

  # -------------------------- 5. Return Result --------------------------
  return(gsea_result)
}

#' Perform gseGO and GSEA Analysis with Visualization
#'
#' Conducts two complementary enrichment analyses—gene set enrichment analysis for GO terms (gseGO)
#' and flexible gene set enrichment analysis (GSEA) using custom or MSigDB gene sets—with integrated
#' visualization (enrichment network for gseGO, enrichment curves for GSEA). Results are returned
#' in a structured list and can be optionally saved to files.
#'
#' @param data Input data for GSEA. Two valid formats:
#'   1. A data frame containing at least a gene symbol column and a fold change/ranking column;
#'   2. A pre-ranked numeric vector where names are gene symbols and values are ranking metrics (e.g., log2FC).
#' @param geneColName Character string. Column name in `data` (if data frame) corresponding to gene symbols.
#'   Required only if `data` is a data frame; ignored if `data` is a numeric vector.
#' @param rangeColName Character string. Column name in `data` (if data frame) corresponding to ranking values
#'   (e.g., log2 fold change, t-statistic). Required only if `data` is a data frame; ignored if `data` is a numeric vector.
#' @param TERM2GENE:
#'   - A data frame with two columns: "term" (gene set names) and "gene" (gene identifiers matching `geneColName` column):
#'     Custom gene set-term mapping for GSEA.

GSEA.baseCustomGeneSet <- function(data,
                                   TERM2GENE,
                                   geneColName = NULL,
                                   rangeColName = NULL) {

  # -------------------------- 1. Validate TERM2GENE --------------------------
  #
  if (!is.data.frame(TERM2GENE) || !all(c("term", "gene") %in% colnames(TERM2GENE))) {
    stop(
      "'TERM2GENE'must be a data frame with columns:\n- 'term': Gene set names\n- 'gene': Gene identifiers matching 'geneColName'",
      call. = FALSE
    )
  }
  # -------------------------- 2. Prepare ranked gene list --------------------------
  if (is.numeric(data) && is.vector(data) && !is.null(names(data))){
    geneList = data
  }else {
    geneList = preRankedGeneList(data = data,geneColName = geneColName, rangeColName = rangeColName)
  }


  # -------------------------- 3. Run GSEA --------------------------
  gsea_result <- clusterProfiler::GSEA(
    geneList = geneList,
    TERM2GENE = TERM2GENE,
    pvalueCutoff = 1,
    verbose = FALSE  # Disable verbose output to avoid clutter
  )
  return(gsea_result)
}


#' exe.gseGO_GSEA
#'
#' @param data Input data for GSEA. Two valid formats:
#'   1. A data frame containing at least a gene symbol column and a fold change/ranking column;
#'   2. A pre-ranked numeric vector where names are gene symbols and values are ranking metrics (e.g., log2FC).
#' @param geneColName Character string. Column name in `data` (if data frame) corresponding to gene symbols.
#'   Required only if `data` is a data frame; ignored if `data` is a numeric vector.
#' @param rangeColName Character string. Column name in `data` (if data frame) corresponding to ranking values
#' @param gseGO.ont For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param keyType For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param OrgDb For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param species For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param TERM2GENE For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param category For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param subcategory For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param showCategory For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param pvalueCutoff For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param fileName For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param save For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param height For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param width For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param folder For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#'
#' @return list
#' @export exe.gseGO_GSEA
#'
exe.gseGO_GSEA <- function(data,
                           geneColName = NULL,
                           rangeColName = NULL,
                           gseGO.ont = "BP",
                           keyType = "SYMBOL",
                           OrgDb = org.Hs.eg.db,
                           species = "Homo sapiens",
                           TERM2GENE = "msigdbr",
                           category = "C5",
                           subcategory = "GO:BP",
                           showCategory = 6,
                           pvalueCutoff=0.01,
                           fileName = "enrich",
                           save = TRUE,
                           height =4,
                           width = 7,
                           folder = "./"){
  ifelse(dir.exists(folder),"FALSE",dir.create(folder,recursive = T))
  # -------------------------- 2. Prepare ranked gene list --------------------------
  if (is.numeric(data) && is.vector(data) && !is.null(names(data))){
    geneList = data
  }else {
    geneList = preRankedGeneList(data = data,geneColName = geneColName, rangeColName = rangeColName)
  }

  set.seed(42)
  message("starting gseGO...........")
  enrich <- gseGO(geneList,
                  keyType = keyType,
                  OrgDb = OrgDb,
                  ont = gseGO.ont)
  p <- aPEAR::enrichmentNetwork(enrich@result[1:min(100,nrow(enrich@result)),])
  if(save == TRUE){
    ifelse(dir.exists(folder),FALSE,dir.create(folder,recursive = T))
    ggsave(filename = paste0(folder,fileName,"-gseGO_enrichmentNetwork.pdf"),
           plot = p,
           height = 6,width = 6)
  }
  gseadata <- list(gseGO = list(),
                   GSEA = list())
  gseadata[["gseGO"]][[gseGO.ont]] <- enrich
  gseadata[["gseGO"]][["enrichmentNetwork"]] <- p
  if(TERM2GENE == "msigdbr"){
    gs = msigdbr::msigdbr(species = species, category = category,subcategory = subcategory)
    gs = gs[,c("gs_name","gene_symbol")]
    colnames(gs) = c("term","gene")
  }else{gs = TERM2GENE}
  egmt <- GSEA(geneList = geneList,
               TERM2GENE = gs,
               verbose = FALSE,
               pvalueCutoff = pvalueCutoff)
  gseadata[["GSEA"]][["GSEA"]] <- egmt
  if(nrow(egmt@result)>0){
    maxnum <- min(showCategory,nrow(egmt@result))
    fig <- gseaNb(egmt,
                  geneSetID = egmt$Description[1:maxnum],
                  curveCol = brewer.pal(maxnum,"Dark2"),
                  subPlot = 2)
    gseadata[["GSEA"]][["Figure"]] <- fig
    if(save == TRUE){
      ggsave(filename = paste0(folder,fileName,"-GSEA.pdf"),
             plot = fig,
             height = height,width = width)
    }
  }
  return(gseadata)
}


