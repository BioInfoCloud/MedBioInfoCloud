Functions and data relevant to the tutorials published on the WeChat official account  [MedBioInfoCloud](https://github.com/BioInfoCloud/MedBioInfoCloud) are encapsulated in this R package. 

## Install package

This guide provides step-by-step instructions to install all required and optional dependencies for the R package, including distinctions between Bioconductor packages (common in bioinformatics) and CRAN packages (general R packages). Follow the steps below to resolve dependency-related errors during installation.

### 1. Prerequisites

Before installing dependencies:

- Ensure your R version is **4.2.0 or higher** (compatible with the latest versions of Bioconductor packages).
- Install `BiocManager` (the official tool for Bioconductor package installation) if you haven’t already:

```R
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")  # Install from CRAN
}
library(BiocManager)  # Load the package
```

### 2. Install Dependencies by Type

Dependencies are split into two categories for clarity. Run the code blocks sequentially to avoid missing packages.

#### 2.1 Bioconductor Packages (Bioinformatics-Specific)

These packages are not available on CRAN and require installation via `BiocManager::install()`. Copy and run the entire code block:

```R
# List of Bioconductor dependencies
bioc_packages <- c(
  "DESeq2",          # Differential expression analysis (RNA-seq)
  "edgeR",           # Differential expression analysis (RNA-seq)
  "GenomicFeatures", # Retrieve genomic annotations (e.g., gene coordinates)
  "GseaVis",         # GSEA result visualization
  "limma",           # Differential expression analysis (microarray/RNA-seq)
  "BiocGenerics",    # Core utility functions for Bioconductor packages
  "GSVA",            # Gene set variation analysis
  "miRBaseVersions.db", # miRBase annotation database
  "TCGAbiolinks",    # Access TCGA (The Cancer Genome Atlas) data
  "aPEAR",           # Alternative polyadenylation analysis
  "clusterProfiler", # Functional enrichment analysis (GO/KEGG)
  "org.Hs.eg.db",    # Human gene annotation database (Entrez ID ↔ Symbol)
  "msigdbr"          # Access MSigDB gene sets
)

# Install Bioconductor packages (skip already installed ones)
BiocManager::install(
  packages = bioc_packages,
  update = FALSE,    # Set to TRUE if you want to update existing packages
  ask = FALSE        # Auto-confirm installation (no manual prompts)
)
```

#### 2.2 CRAN Packages (General R Packages)

```R
# List of CRAN dependencies
cran_packages <- c(
  "dplyr",           # Data manipulation (filter, group_by, etc.)
  "ggplot2",         # Core visualization package
  "ggpubr",          # Publication-ready ggplot2 plots
  "ggrepel",         # Avoid text overlap in plots
  "ggsignif",        # Add significance labels (e.g., * p<0.05) to plots
  "purrr",           # Functional programming (e.g., map() for iteration)
  "magrittr",        # Pipe operator (%>%) for readable code
  "RColorBrewer",    # Color palettes for visualization
  "reshape2",        # Data reshaping (melt(), dcast())
  "data.table"       # Fast data frame manipulation
)

# Install CRAN packages (skip already installed ones)
install.packages(
  pkgs = cran_packages,
  repos = "https://cloud.r-project.org/",  # Use CRAN mirror for fast installation
  dependencies = TRUE,                     # Install required sub-dependencies
  quiet = FALSE                            # Show installation progress
)
```

### 3. Verify Successful Installation

After installation, run the following code to confirm all dependencies are installed and loadable. If no errors appear, the dependencies are ready:

```R
# Function to check if a package is installed and loadable
check_package <- function(pkg) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    message(sprintf("✓ %s is installed and loadable", pkg))
  } else {
    warning(sprintf("✗ %s is missing or failed to load", pkg), call. = FALSE)
  }
}

# Check all dependencies
all_packages <- c(bioc_packages, cran_packages)
invisible(lapply(all_packages, check_package))
```

### 4. Troubleshooting Common Errors

If you encounter installation failures, try these fixes:

#### 4.1 "Package Not Available" Error

- For Bioconductor packages: Ensure `BiocManager` is up-to-date with `BiocManager::install(update = TRUE)`.
- For CRAN packages: Use a different CRAN mirror (e.g., `repos = "https://cran.r-project.org/"`).

#### 4.2 Compilation Errors (e.g., "gcc not found")

- **Windows**: Install [Rtools](https://cran.r-project.org/bin/windows/Rtools/) (required for compiling C/C++ code in packages like `data.table`).
- **macOS**: Install Xcode Command Line Tools with `xcode-select --install` (in Terminal).
- **Linux**: Install system libraries via `sudo apt-get install r-base-dev` (Debian/Ubuntu) or `sudo dnf install R-devel` (Fedora).

#### 4.3 Out-of-Memory Errors

- Close unused R sessions or applications to free memory.
- Install packages one at a time (e.g., `BiocManager::install("DESeq2")` followed by `BiocManager::install("edgeR")`) to reduce memory usage.

### 5. Proceed to Install the Target R Package

Once all dependencies are successfully installed, you can install the target R package (e.g., from a local file or GitHub):

#### Example 1: Install from a Local `.tar.gz` File

```R
install.packages(
  pkgs = "path/to/your/package.tar.gz",  # Replace with your package path
  repos = NULL,
  type = "source",
  dependencies = TRUE  # Ensure dependencies are linked
)
```

#### Example 2: Install from GitHub (if applicable)

```R
# Install devtools first (if needed)
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
devtools::install_github("username/package-name")  # Replace with GitHub repo
```

### 6. Install MedBioInfoCloud

```R
devtools::install_github("BioInfoCloud/MedBioInfoCloud")
```

Some functions with the same parameters are uniformly explained as follows:

**save** - A boolean value that determines whether the data should be preserved. If `save` is set to `TRUE`, the data will be  saved. Conversely, if it's set to `FALSE`, the data will unsaved.

**folder** - A string that designates the destination directory for the saved data. This parameter is applicable only when `save` is set to `TRUE`. If the specified directory does not exist, it will be automatically created. Kindly note that the folder path should not terminate with a slash ("/"). If the function does not have a `save` parameter but includes a `folder` parameter, it indicates that the function will definitely save data, and may not necessarily return valid data. The `folder` parameter defaults to the current working directory, and it is recommended to set a specific data output directory.

For more learning materials, please refer to [MedBioInfoCloud - 生物信息云](https://bioinfocloud.github.io/note/package/MedBioInfoCloud/).

=========================================================================

