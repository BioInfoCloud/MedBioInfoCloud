Functions and data relevant to the tutorials published on the WeChat official account  [MedBioInfoCloud](https://github.com/BioInfoCloud/MedBioInfoCloud) are encapsulated in this R package. 

## Install package

```R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
DependencyPackage <- c("edgeR","DESeq2","TCGAbiolinks")
BiocManager::install(DependencyPackage)
# install.packages("devtools")
devtools::install_github("BioInfoCloud/MedBioInfoCloud")
```

The common parameters for all functions in this package are:

**save** - A boolean value that determines whether the data should be preserved. If `save` is set to `TRUE`, the data will be  saved. Conversely, if it's set to `FALSE`, the data will unsaved.

**folder** - A string that designates the destination directory for the saved data. This parameter is applicable only when `save` is set to `TRUE`. If the specified directory does not exist, it will be automatically created. Kindly note that the folder path should not terminate with a slash ("/").

For more learning materials, please refer to [MedBioInfoCloud - 生物信息云](https://bioinfocloud.github.io/note/WeChatOfficialAccount/MedBioInfoCloud/).

======================================================
