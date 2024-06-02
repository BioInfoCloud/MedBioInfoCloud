Functions and data relevant to the tutorials published on the WeChat official account  [MedBioInfoCloud](https://github.com/BioInfoCloud/MedBioInfoCloud) are encapsulated in this R package. 

## Install package

This R package has dependencies on several other packages in order to function properly. If you encounter errors during the installation process, it is typically because the required dependencies have not been properly installed. Please install the dependencies individually before proceeding with the installation of the R package.

```R
if (!requireNamespace("BiocManager", quietly = TRUE))  
    install.packages("BiocManager")  
DependencyPackage <- c("edgeR", "DESeq2", "TCGAbiolinks")  
BiocManager::install(DependencyPackage)  
   
if (!requireNamespace("devtools", quietly = TRUE))  
    install.packages("devtools")  
devtools::install_github("BioInfoCloud/MedBioInfoCloud")
```

Some functions with the same parameters are uniformly explained as follows:

**save** - A boolean value that determines whether the data should be preserved. If `save` is set to `TRUE`, the data will be  saved. Conversely, if it's set to `FALSE`, the data will unsaved.

**folder** - A string that designates the destination directory for the saved data. This parameter is applicable only when `save` is set to `TRUE`. If the specified directory does not exist, it will be automatically created. Kindly note that the folder path should not terminate with a slash ("/"). If the function does not have a `save` parameter but includes a `folder` parameter, it indicates that the function will definitely save data, and may not necessarily return valid data. The `folder` parameter defaults to the current working directory, and it is recommended to set a specific data output directory.

For more learning materials, please refer to [MedBioInfoCloud - 生物信息云](https://bioinfocloud.github.io/note/WeChatOfficialAccount/MedBioInfoCloud/).

=========================================================================

