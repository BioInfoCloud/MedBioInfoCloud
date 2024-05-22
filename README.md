Functions and data relevant to the tutorials published on the WeChat official account  [MedBioInfoCloud](https://github.com/BioInfoCloud/MedBioInfoCloud) are encapsulated in this R package. For more learning materials, please refer to [BioInfoCloud](https://bioinfocloud.github.io/note/).

## 一. Install package

```R
# install.packages("devtools")
devtools::install_github("BioInfoCloud/MedBioInfoCloud")
```

## 二. TCGA数据库数据挖掘相关函数

TCGA数据库首页：[GDC Data Portal Homepage (cancer.gov)](https://portal.gdc.cancer.gov/)

### 1.下载RNAseq数据

getTCGA_RNAseqData()返回一个list，包括count,tpm和fpkm 3个数据框。

```R
STARdata <- getTCGA_RNAseqData("TCGA-LUAD")
```

### 2.下载蛋白组数据

getProteinExp()返回一个数据框。

```R
Proteome_data <- getProteinExp("TCGA-LUAD")
```

### 3.下载SNV（simple nucleotide variation）数据

数据类型为：Masked Somatic Mutation。

```R
snv.dat <- getSNV_Masked_data("TCGA-LUAD")
```

### 4.下载miRNA数据

参考文章：

[https://mp.weixin.qq.com/s/__EjCrJFc08itoF3xqawNg](https://mp.weixin.qq.com/s/__EjCrJFc08itoF3xqawNg)

[https://mp.weixin.qq.com/s/-FH0Vi4PaCjhPbEq4-lxbg](https://mp.weixin.qq.com/s/-FH0Vi4PaCjhPbEq4-lxbg)

[https://mp.weixin.qq.com/s/WxgMhwpMAJy_CKTqNdFj0g](https://mp.weixin.qq.com/s/WxgMhwpMAJy_CKTqNdFj0g)

#### （1）Isoform Expression Quantification

```R
IsoformEQ <- get_miRNA_IsoformEQ("TCGA-LUAD")
```

#### （2）miRNA Expression Quantification

```R
miRNAEQ <- get_miRNAEQ("TCGA-LUAD")
```

### 5.下载甲基化数据

getMetData下载Methylation Beta Value数据。

```R
MetData <- getMetData("TCGA-LUAD")
```

### 6. 下载CNV（Copy Number Variation）数据

getCNV.data()函数还在优化中：

```R
cnv.gl <- getCNV.data("TCGA-LUAD",save = FALSE,folder = ".",data.type = "Gene Level Copy Number")
```

```R
cnv.gls <-getCNV.data("TCGA-LUAD",save = FALSE,folder = ".",data.type = "Gene Level Copy Number Scores")
```

### 7. 下载临床数据

```R
cldat <- getClinicalData(project = "TCGA-LUAD",save = FALSE,folder = ".",trim = TRUE)
```

针对的癌症类型：

```R
c("TCGA-READ","TCGA-COAD","TCGA-PAAD","TCGA-ESCA","TCGA-KIRP","TCGA-HNSC",
             "TCGA-BLCA","TCGA-STAD","TCGA-CHOL","TCGA-SKCM","TCGA-LUAD","TCGA-LIHC",
             "TCGA-KIRC","TCGA-KICH","TCGA-MESO","TCGA-LUSC","TCGA-GBM","TCGA-UVM",
             "TCGA-BRCA","TCGA-TGCT","TCGA-THCA")
```

由于每种癌症类型的临床信息有差异，其他癌症类型，获取临床数据可能会报错，可以通过指定getClinicalData()中的trim = FALSE，返回原始未整理过的数据。

```R
cldat <- getClinicalData(project = "TCGA-LUAD",save = FALSE,folder = ".",trim = FALSE)
```

### 8. 过滤表达数据

filterGeneTypeExpr()根据某列里面是数据进行过滤，保留filter值的数据。该函数仅适用于getTCGA_RNAseqData获取的count,tpm和fpkm 3个数据框。

```R
STARdata <- getTCGA_RNAseqData("TCGA-LUAD")
expr <- STARdata[["count"]]
table(expr$gene_type)
pc.expr <- filterGeneTypeExpr(expr = expr,fil_col = "gene_type",filter = "protein_coding")
```

### 9.分割RNAseq表达数据

splitTCGAmatrix()，data的列应该是TCGA病人样本的barcode，参数sample的值为"Tumor"或"Normal"，指定sample ="Normal"时，当样本中没有正常样本返回NULL。

```R
turexp <- splitTCGAmatrix(data = expr[,-c(1:3)],sample = "Tumor")
norexp <- splitTCGAmatrix(data = expr[,-c(1:3)],sample = "Normal")
```

## 10. 数据打包下载

下载的数据是R对象：

RNAseq：[微信公众号生物信息云提供的链接](https://pan.baidu.com/s/1VWz8bIlgKaUKR0ncughBhg?pwd=e6wz )

蛋白组：[微信公众号生物信息云提供的链接](https://pan.baidu.com/s/1CrO2jIrXh-R1L9hfuO-ESQ?pwd=ogqx) 

TCGA-miRNA_Isoform：[微信公众号生物信息云提供的链接](https://pan.baidu.com/s/1k8-ZTwbsjQRE49EgORWUxQ?pwd=mx43 )

Survival和Phenotype数据（fromUCSC）：[微信公众号生物信息云提供的链接](https://pan.baidu.com/s/1_VmOO_yyjiaEkLWlHxRYWg?pwd=04au)

临床数据：[微信公众号生物信息云提供的链接](https://pan.baidu.com/s/1KDO2gx-lnejeuInVZSEPFQ?pwd=0k83)

### 11 . 获取某个基因在泛癌中的表达数据

geneSymbol是要分析的基因名称；dataType是tpm,fpkm和count中的一种；datafolder来自 getTCGA_RNAseqData()函数下载数据，并存放在某个文件夹中，或者从这里下载（RNAseq：[微信公众号生物信息云提供的链接](https://pan.baidu.com/s/1VWz8bIlgKaUKR0ncughBhg?pwd=e6wz )），但这里下载的数据没有fpkm；geneType参照函数filterGeneTypeExpr()中的fil_col，pattern正则表达式匹配datafolder中的数据文件；nnorm表示至少包含几个正常样本；得到的数据进行了log2转换。

```R
geneSymbol = "ATG7"
datafolder = "G:/DatabaseData/TCGA/new/processedTCGAdata/TCGA-STAR_Exp"
df = getGeneExpData.fancancer(datafolder,
                              geneSymbol,
                              geneType = "protein_coding",
                              dataType = "tpm",
                              pattern = "STARdata.Rdata$",
                              nnorm = 10)
```

得到的数据样式如下：

![](https://raw.githubusercontent.com/BioInfoCloud/ImageGo/main/20240521220111.png)

## 三.基础分析相关函数的使用

### 1.差异表达分析

geneDEAnalysis：

```R
geneDEAnalysis <- function (data, group, comparison,method = "DESeq2", filter = TRUE)
```

data为表达数据，行为基因名称，列为样本名称，group是一个数据框，只有一列为group的值，其值是二分类的字符串标签（如：Tumor，Normal），行名为样本名称，其顺序与data的列名一致。comparison是一个由group中的二分类标签值用-链接，如"Tumor-Normal"，表示Tumor组与Normal进行差异表达分析。method是DESeq2, edgeR和 limma中的一种，RNAseq数据建议使用DESeq2或edgeR，芯片数据使用limma。filter是否过滤数据，默认为TRUE。

