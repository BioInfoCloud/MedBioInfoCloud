Some functions or data involved in the relevant tutorials published on the WeChat official account [MedBioInfoCloud](https://github.com/BioInfoCloud/MedBioInfoCloud) are encapsulated in this R package.

## install package

```R
# install.packages("devtools")
devtools::install_github("BioInfoCloud/MedBioInfoCloud")
```

TCGA数据库数据下载

下载的数据是R对象：

RNAseq：[微信公众号生物信息云提供的链接](https://pan.baidu.com/s/1VWz8bIlgKaUKR0ncughBhg?pwd=e6wz )，提取码：e6wz 
蛋白组：[微信公众号生物信息云提供的链接](https://pan.baidu.com/s/1CrO2jIrXh-R1L9hfuO-ESQ?pwd=ogqx) ，提取码：ogqx 
TCGA-miRNA_Isoform：[微信公众号生物信息云提供的链接](https://pan.baidu.com/s/1k8-ZTwbsjQRE49EgORWUxQ?pwd=mx43 )

## 相关函数的使用

### geneDEAnalysis：

```R
geneDEAnalysis <- function (data, group, comparison,method = "DESeq2", filter = TRUE)
```

data为表达数据，行为基因名称，列为样本名称，group是一个数据框，只有一列为group的值，其值是二分类的字符串标签（如：Tumor，Normal），行名为样本名称，其顺序与counts的列名一致。comparison是一个由group中的二分类标签值用-链接，如"Tumor-Normal"，表示Tumor组与Normal进行差异表达分析。method是DESeq2, edgeR和 limma中的一种，RNAseq数据建议使用DESeq2或edgeR，芯片数据使用limma。filter是否过滤数据，默认为TRUE。

