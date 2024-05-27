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
Proteome_data <- getTCGA_ProteinExp("TCGA-LUAD")
```

### 3.下载SNV（simple nucleotide variation）数据

数据类型为：Masked Somatic Mutation。

```R
snv.dat <- getTCGA_SNV_Masked_data("TCGA-LUAD")
```

### 4.下载miRNA数据

参考文章：

[https://mp.weixin.qq.com/s/__EjCrJFc08itoF3xqawNg](https://mp.weixin.qq.com/s/__EjCrJFc08itoF3xqawNg)

[https://mp.weixin.qq.com/s/-FH0Vi4PaCjhPbEq4-lxbg](https://mp.weixin.qq.com/s/-FH0Vi4PaCjhPbEq4-lxbg)

[https://mp.weixin.qq.com/s/WxgMhwpMAJy_CKTqNdFj0g](https://mp.weixin.qq.com/s/WxgMhwpMAJy_CKTqNdFj0g)

#### （1）Isoform Expression Quantification

```R
IsoformEQ <- getTCGA_miRNA_IsoformEQ("TCGA-LUAD")
```

#### （2）miRNA Expression Quantification

```R
miRNAEQ <- getTCGA_miRNAEQ("TCGA-LUAD")
```

### 5.下载甲基化数据

getTCGA_MethylationData 下载Methylation Beta Value数据。

```R
MetData <- getTCGA_MethylationData("TCGA-LUAD")
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

### 9.分割数据

splitTCGAmatrix()，data的列应该是TCGA病人样本的barcode，参数sample的值为"Tumor"或"Normal"，指定sample ="Normal"时，当样本中没有正常样本返回NULL。

```R
turexp <- splitTCGAmatrix(data = expr[,-c(1:3)],sample = "Tumor")
norexp <- splitTCGAmatrix(data = expr[,-c(1:3)],sample = "Normal")
```

### 10. 删除重复病人样本

del_dup_sample()函数可以将列为barcode的数据，去除有重复的数据，TCGA数据库的病人有的可能做了几个重复。可以只需要一个。

```R
expr <- del_dup_sample(expr = pc.expr,col_rename = TRUE)
```

### 11. 数据打包下载

下载的数据是R对象：

RNAseq：[微信公众号生物信息云提供的链接](https://pan.baidu.com/s/1VWz8bIlgKaUKR0ncughBhg?pwd=e6wz )

蛋白组：[微信公众号生物信息云提供的链接](https://pan.baidu.com/s/1CrO2jIrXh-R1L9hfuO-ESQ?pwd=ogqx) 

TCGA-miRNA_Isoform：[微信公众号生物信息云提供的链接](https://pan.baidu.com/s/1k8-ZTwbsjQRE49EgORWUxQ?pwd=mx43 )

Survival和Phenotype数据（fromUCSC）：[微信公众号生物信息云提供的链接](https://pan.baidu.com/s/1_VmOO_yyjiaEkLWlHxRYWg?pwd=04au)

临床数据：[微信公众号生物信息云提供的链接](https://pan.baidu.com/s/1KDO2gx-lnejeuInVZSEPFQ?pwd=0k83)

### 12 . 获取某个基因在泛癌中的表达数据

geneSymbol是要分析的基因名称的向量；dataType是tpm,fpkm和count中的一种；datafolder来自 getTCGA_RNAseqData()函数下载数据，并存放在某个文件夹中，或者从这里下载（RNAseq：[微信公众号生物信息云提供的链接](https://pan.baidu.com/s/1VWz8bIlgKaUKR0ncughBhg?pwd=e6wz )），但这里下载的数据没有fpkm；geneType参照函数filterGeneTypeExpr()中的fil_col，pattern正则表达式匹配datafolder中的数据文件；paired指定是否只获取配对样本的数据；nnorm表示至少包含几个正常样本；得到的数据进行了log2转换。

```R
geneSymbol = c("ATG7","ATG12")
datafolder = "G:/DatabaseData/TCGA/new/processedTCGAdata/TCGA-STAR_Exp"
df = getGeneExpData.pancancer(datafolder,
                              geneSymbol,
                              geneType = "protein_coding",
                              dataType = "tpm",
                              pattern = "STARdata.Rdata$",
                              paired = FALSE,
                              nnorm = 10)
```

得到的数据样式如下：

![](https://raw.githubusercontent.com/BioInfoCloud/ImageGo/main/20240522130053.png)

### 13. 单基因在泛癌中表达的箱型图可视化

data是由 getGeneExpData.fancancer得到的数据，gene是一个基因，字符串类型；paired表示数据是否是配对样本。

```R
fig <- ggplotGenePancancerExp(data = df,gene= "ATG7",
                              save = FALSE,folder = ".",paired = FALSE)
```

![](https://raw.githubusercontent.com/BioInfoCloud/ImageGo/main/20240522131135.png)

## 三.一些数据处理和分析过程中的基础函数

### 1.输出gmt文件

outputGmtFile()函数中description默认为NA，如果指定，应该是一个长度与input相同，用于描述每个基因集的字符串向量。filename应该是一个.gmt结尾的文件名称，可包括路径。input是一个list或是一个data.frame，如果是list，list中每一个对象是一个向量（基因），每一个对象应该有一个合适的名称，相当于基因集的名称，下面是一个input接收list数据对象案例：

![](https://raw.githubusercontent.com/BioInfoCloud/ImageGo/main/20240523205319.png)

```R
outputGmtFile(input = genes,description = NA,filename = "./gs.gmt")
```

如果input是一个数据框，应该包括两列，第一列是关于基因集描述的术语，第二列是基因，这时，如果需要指定description的值，长度应该等于input第一列值作为集合的长度。下面是一个input输入作为数据框的案例：

![](https://raw.githubusercontent.com/BioInfoCloud/ImageGo/main/20240523211219.png)

```
outputGmtFile(input = cg,description = NA,filename = "./gs.gmt")
```

### 2.融合生存数据与特征数据

```R
se <- mergeSurExp(expr
                  ,survival
                  ,survivalFrome = NULL
                  ,Timeunit=1
                  ,TCGA = FALSE
                  ,TCGAfrome = "MedBioInfoCloud"
                  ,feature = NULL
                  ,save = FALSE
                  ,folder = "."
)
```

如果不是TCGA数据库的数据，只需要关注参数expr，survival，save，folder，Timeunit，其他参数不需要考虑，并且，expr行为特征（一般为基因），列为样本；Timeunit的值表示生存时间进行何种转换，Timeunit=1表示不进行任何转换，如果你的生存数据的时间是天，可设置Timeunit=365，转换为年；feature是一个特征子集向量，可以不指定，默认expr的所有行。如果处理TCGA的数据，TCGA应该指定为TRUE，expr应该是getTCGA_RNAseqData()返回结果中的表达数据，如下图：

![](https://raw.githubusercontent.com/BioInfoCloud/ImageGo/main/20240525165216.png)

这时，如果生存数据是自己整理的，行名应该和expr的样本一致，或者有交集，如果数据是来自Survival和Phenotype数据（fromUCSC）：[微信公众号生物信息云提供的链接](https://pan.baidu.com/s/1_VmOO_yyjiaEkLWlHxRYWg?pwd=04au)，可以直接指定survivalFrome = "UCSC2022"，数据样式如下：

![](https://raw.githubusercontent.com/BioInfoCloud/ImageGo/main/20240525171631.png)

 如果从这里下载临床数据：[微信公众号生物信息云提供的链接](https://pan.baidu.com/s/1KDO2gx-lnejeuInVZSEPFQ?pwd=0k83)，或是通过本包getClinicalData()函数【trim = TRUE】获取的数据，可以直接指定survivalFrome = "GDCquery_clinic"，数据样式如下（至少包含"submitter_id","vitalStat","surTime")：

![](https://raw.githubusercontent.com/BioInfoCloud/ImageGo/main/20240525171820.png)

如果指定save = TRUE，会在指定的folder文件夹下保存csv和Rdata的格式文件。

```R
STARdata <- getTCGA_RNAseqData("TCGA-LUAD")
cldat <- getClinicalData(project = "TCGA-LUAD",save = FALSE,folder = ".",trim = TRUE)
se <- mergeSurExp(expr = STARdata[["tpm"]],
                  survival = cldat,
                  survivalFrome = "GDCquery_clinic",
                  TCGA = TRUE)
```



## 四.基础分析相关函数的使用

### 1.差异表达分析

data为表达数据，行为基因名称，列为样本名称，group是一个数据框，只有一列为group的值，其值是二分类的字符串标签（如：Tumor，Normal），行名为样本名称，其顺序与data的列名一致。comparison是一个由group中的二分类标签值用-链接，如"Tumor-Normal"，表示Tumor组与Normal进行差异表达分析。method是DESeq2, edgeR和 limma中的一种，RNAseq数据建议使用DESeq2或edgeR，芯片数据使用limma。filter是否过滤数据，默认为TRUE。

```R
DEG <- geneDEAnalysis(data, group, comparison,method = "DESeq2", filter = TRUE)
```

## 五.预后模型构建相关函数

### 1.特征选择

featureSelect.baseSur()函数可以基于lasso回归，随机森林以及单因素COX回归进行特征选择。

```R
fs <- featureSelect.baseSur(data
                            ,dataFrom = NULL,
                            feature ="all"
                            ,method = "all"
                            ,cutoff = 0.05
                            ,save = TRUE
                            ,folder = ".")
```

data是一个数据框，列应该包括生存数据和特征，行为样本。如果数据来源mergeSurExp()函数，前3列的列名应该如下图所示，并且需要指定dataFrom = "mergeSurExp"，如果是自己整理的数据，第一列的列名可以随意，但第2列（生存状态）和第3列（生存时间）的列名必须与图中相同。如果没有第一列，将设置dataFrom = NULL。feature，表示要进行分析的特征，默认是所有特征（下图中除了前3列）；method的值有lasso、cox、randomForest，单独设置这3个值时，函数返回一个向量（即筛选出的特征），method="all"时，3种方法都执行，最后返回一个list，包括3种方法的结果。cutoff只有当method="cox"或method="all"时被使用。

![](https://raw.githubusercontent.com/BioInfoCloud/ImageGo/main/20240526115531.png)

### 2.多因素COX回归模型

MultivariateCOX()函数用于一键式构建多因素COX回归模型。

```R
MultivariateCOX(data
                ,dataFrom ="mergeSurExp",
                feature ="all"
                ,method = "all"
                ,train_prop = 0.8
                ,cutoff = 0.05
                ,save = TRUE
                ,folder = ".")
```

