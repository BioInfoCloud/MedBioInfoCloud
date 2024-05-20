Some functions or data involved in the relevant tutorials published on the WeChat official account [MedBioInfoCloud](https://github.com/BioInfoCloud/MedBioInfoCloud) are encapsulated in this R package.

## Install package

```R
# install.packages("devtools")
devtools::install_github("BioInfoCloud/MedBioInfoCloud")
```

## TCGA数据库数据挖掘相关函数

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

### 4.miRNA数据下载

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



### 过滤表达数据

filterGeneTypeExpr()根据某列里面是数据进行过滤，保留filter值的数据。该函数仅使用于getTCGA_RNAseqData获取的count,tpm和fpkm 3个数据框。

```R
STARdata <- getTCGA_RNAseqData("TCGA-LUAD")
expr <- STARdata[["count"]]
table(expr$gene_type)
pc.expr <- filterGeneTypeExpr(expr = expr,fil_col = "gene_type",filter = "protein_coding")
```

## 数据下载

下载的数据是R对象：

RNAseq：[微信公众号生物信息云提供的链接](https://pan.baidu.com/s/1VWz8bIlgKaUKR0ncughBhg?pwd=e6wz )

蛋白组：[微信公众号生物信息云提供的链接](https://pan.baidu.com/s/1CrO2jIrXh-R1L9hfuO-ESQ?pwd=ogqx) 

TCGA-miRNA_Isoform：[微信公众号生物信息云提供的链接](https://pan.baidu.com/s/1k8-ZTwbsjQRE49EgORWUxQ?pwd=mx43 )

Survival和Phenotype数据（fromUCSC）：[微信公众号生物信息云提供的链接](https://pan.baidu.com/s/1_VmOO_yyjiaEkLWlHxRYWg?pwd=04au)


## 相关函数的使用

### geneDEAnalysis：

```R
geneDEAnalysis <- function (data, group, comparison,method = "DESeq2", filter = TRUE)
```

data为表达数据，行为基因名称，列为样本名称，group是一个数据框，只有一列为group的值，其值是二分类的字符串标签（如：Tumor，Normal），行名为样本名称，其顺序与counts的列名一致。comparison是一个由group中的二分类标签值用-链接，如"Tumor-Normal"，表示Tumor组与Normal进行差异表达分析。method是DESeq2, edgeR和 limma中的一种，RNAseq数据建议使用DESeq2或edgeR，芯片数据使用limma。filter是否过滤数据，默认为TRUE。

