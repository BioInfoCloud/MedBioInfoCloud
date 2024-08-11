


#' featureSelect.randomForest
#'
#' @param data For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param dataFrom For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param save For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param feature For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param folder For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#'
#' @return vector
#' @export featureSelect.randomForest
#'
featureSelect.randomForest <- function(data
                                       ,dataFrom ="mergeSurExp"
                                       ,save = TRUE
                                       ,feature = "all"
                                       ,folder = "."){

  if(!is.null(dataFrom) & dataFrom == "mergeSurExp"){
    data <- data[,-1]
  }
  if(!(feature == "all")){
    fs <- intersect(colnames(data),feature)
    data <- data[,c("vitalStat","surTime",fs)]
  }
  rsf <- randomForestSRC::rfsrc(Surv(surTime, vitalStat) ~ . ,
                                data = data,
                                ntree = 100, ## tree的个数
                                nsplit = 1) ## 随机拆分数（非负整数），可提高运算速度，默认值为0
  # # 设置树的个数为3
  # plot(randomForestSRC::get.tree(rsf, 3))
  # # plot(rsf)
  # randomForestSRC::print.rfsrc(rsf) ## 输出结果信息
  # pred <- predict(rsf,data,OOB=TRUE,type="response")
  ## two options : 'vh':速度慢，准确度高；'md':速度快，准确度相对较低。
  select.fs <- randomForestSRC::var.select(object = rsf,method="vh")
  rf.select.fs <- select.fs$topvars ## 提取分析结果中方差较大的变量

  if(save == TRUE){
    save(rf.select.fs,file = paste0(folder,"/randomForest.selectedFeature.Rdata"))
    writeLines(rf.select.fs,con = paste0(folder,"/randomForest.selectedFeature.txt"))
  }
  return(rf.select.fs)
}


#' UnivariateCOX
#'
#' @param data For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param variable For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#'
#' @return data.frame
#' @export survival
#'
UnivariateCOX <- function(data,variable){ ## 构建一个R function 便于后期调用
  FML <- as.formula(paste0('BaSurv~',variable)) ## 构建生存分析公式
  GCox <- coxph(FML, data = data) ## Cox分析
  GSum <- summary(GCox) ## 输出结果
  HR <- round(GSum$coefficients[,2],2) ## 输出HR值
  PValue <- round(GSum$coefficients[,5],3) ## 输出P值
  CI <- paste0(round(GSum$conf.int[,3:4],2),collapse = "-") ## 输出HR的执行区间
  Unicox <- data.frame("characteristics" = variable, ## 返回结果，并构建数据框
                       "Hazard Ratio" = HR,
                       "CI95" = CI,
                       HR.95L = GSum$conf.int[,"lower .95"],
                       HR.95H = GSum$conf.int[,"upper .95"],
                       "P Value" = PValue)
  return(Unicox)
}

#' featureSelect.UnivariateCox
#'
#' @param data For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param dataFrom For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param feature For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param cutoff For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param save For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param folder For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#'
#' @return vector
#' @export featureSelect.UnivariateCox
#'
featureSelect.UnivariateCox <- function(data
                              ,dataFrom ="mergeSurExp",
                              feature ="all"
                              ,cutoff = 0.05
                              ,save = TRUE
                              ,folder = "."){
  if(!is.null(dataFrom) & dataFrom == "mergeSurExp"){
    data <- data[,-1]
  }
  if(!(feature == "all")){
    VarNames <- intersect(colnames(data),feature)
    data <- data[,c("vitalStat","surTime",fs)]
  }else{
    VarNames <- colnames(data)[-c(1:2)]
  }
  BaSurv <- survival::Surv(time = data$surTime, ## 生存时间
                 event = data$vitalStat) ## 生存状态
  UniVar <- lapply(VarNames,function(x){
    UnivariateCOX(data,variable = x)
  }) ## 批量做Cox分析
  UniVar <- plyr::ldply(UniVar,data.frame) ## 将结果整理为一个数据框
  uniCox.select.fs <- UniVar$characteristics[which(UniVar$P.Value < cutoff)] %>% as.character() ## 筛选其中P值<0.2的变量纳入多因素cox分析。

  if(save == TRUE){
    save(uniCox.select.fs,file = paste0(folder,"/uniCox.selectedFeature.Rdata"))
    writeLines(uniCox.select.fs,con = paste0(folder,"/uniCox.selectedFeature.txt"))
    write.table(UniVar,file=paste0(folder,"/all.feature.coxResult.txt"),sep="\t",
                row.names=F,quote=F)
  }
  return(uniCox.select.fs)
}



#' cox.forestplot
#'
#' @param forest_data For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param ylab For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param save For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param folder For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#'
#' @return list
#' @export forestplot,viridis
#'
cox.forestplot <- function(forest_data,ylab,save = FALSE,folder = "."){
  #####==================森林图样式一
  HR <- forest_data$Hazard.Ratio
  CI_LL  <- forest_data$HR.95L
  CI_HL <-forest_data$HR.95H
  style1 <- ggplot2::ggplot(data=forest_data,
                            aes(x=Hazard.Ratio,y= characteristics,
                                color= P.Value))+
    geom_errorbarh(aes(xmax=CI_HL, xmin=CI_LL),
                   color="skyblue",height=0.2,size=1.2)+
    geom_point(aes(x=Hazard.Ratio,y= characteristics),size=4,shape=18)+
    geom_vline(xintercept = 1,linetype='dashed',linewidth=1.2, color = "skyblue")+
    labs(x = "Hazard Ratio",y = ylab)+
    viridis::scale_color_viridis()+
    theme_bw()


  ###=================样式二==========
  forest_table <- cbind(c(ylab, forest_data$characteristics),
                        c("HR (95% CI)", forest_data$Hazard.Ratio),
                        c("p-value", forest_data$P.Value))

  csize <- data.frame(mean=c(NA, as.numeric(forest_data$Hazard.Ratio)),
                      lower=c(NA, as.numeric(forest_data$HR.95L)),
                      upper=c(NA, as.numeric(forest_data$HR.95H)))
  head(csize)
  style2 <- forestplot::forestplot(labeltext = forest_table,
                                   mean=csize[,1],
                                   lower=csize[,2],
                                   upper=csize[,3],
                                   graph.pos = 3,
                                   graphwidth = unit(5, "cm"),
                                   zero = 1,
                                   cex = 2,
                                   lineheight = "auto",
                                   boxsize = 0.2,
                                   fn.ci_norm = fpDrawNormalCI,
                                   lwd.ci = 1,
                                   ci.vertices = TRUE,
                                   lwd.xaxis = 1,
                                   txt_gp=fpTxtGp(label=gpar(cex=1.25),#各种字体大小设置
                                                  ticks=gpar(cex=1.25),
                                                  xlab=gpar(cex = 1.25),
                                                  title=gpar(cex = 1.25)),
                                   xlab = "Hazard Ratio",
                                   ci.vertices.height = 0.08,
                                   font= 1.3,
                                   col = fpColors(box = "royalblue",
                                                  line = "darkblue"))
  if(save == TRUE){
    ggsave(filename =  paste0(folder,"/forestplot.style-1.Rdata"),plot = style1)
    ggsave(filename =  paste0(folder,"/forestplot.style-2.Rdata"),plot = style2)
  }
  return(list(style1 = style1,style2 = style2))
}

#' featureSelect.lasso
#'
#' @param data For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param dataFrom For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param feature For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param save For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param folder For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#'
#' @return vector
#' @export glmnet
#'
featureSelect.lasso <- function(data
                                ,dataFrom ="mergeSurExp",
                                feature ="all"
                                ,save = TRUE
                                ,folder = "."){
  if(!is.null(dataFrom) & dataFrom == "mergeSurExp"){
    data <- data[,-1]
  }
  if(!(feature == "all")){
    VarNames <- intersect(colnames(data),feature)
    data <- data[,c("vitalStat","surTime",fs)]
  }else{
    VarNames <- colnames(data)[-c(1:2)]
  }
  #### Step II 迭代特征选取 ####
  ## Chr I 随机森林（生存分析) ##
  fs_expr <- as.matrix(data[,-c(1:2)]) ## 根据R包的要求，将数据需要筛选的部分提取转换为矩阵,lasso不需要生存时间和生存状态，所以多删除两列

  response <- data.matrix(Surv(data$surTime,data$vitalStat))

  fit <- glmnet::glmnet(fs_expr
                        ,response
                        ,family = "cox"
                        ,maxit = 5000
                        ,alpha = 1
                        )

  cvfit = glmnet::cv.glmnet(x = fs_expr
                    ,Surv(data$surTime,data$vitalStat)
                    ,nfold=10#10倍交叉验证，非必须限定条件.
                    ,family = "cox"
                    ,alpha = 1
                    )

  coef.min = coef(cvfit, s = "lambda.min")  ## lambda.min & lambda.1se 取一个
  active.min = which(coef.min != 0 ) ## 找出那些回归系数没有被惩罚为0的
  lasso.select.fs <- colnames(fs_expr)[active.min]## 提取基因名称

  if(save == TRUE){
    save(lasso.select.fs,file = paste0(folder,"/lasso.selectedFeature.Rdata"))
    writeLines(lasso.select.fs,con = paste0(folder,"/lasso.selectedFeature.txt"))
    pdf(file = paste0(folder,"/lasso.cvfit.pdf"),height = 3,width = 4.5)
    plot(cvfit) ## 画图
    plot(fit, label = TRUE)
    dev.off()
  }
  return(lasso.select.fs)
}

#' featureSelect.baseSur
#'
#' @param data For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param dataFrom For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param feature For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param method For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param cutoff For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param save For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param folder For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#'
#' @return list or vector
#' @export featureSelect.baseSur
#'
featureSelect.baseSur <- function(data
                                  ,dataFrom ="mergeSurExp",
                                  feature ="all"
                                  ,method = "all"
                                  ,cutoff = 0.05
                                  ,save = TRUE
                                  ,folder = "."){
  if(method == "all"){
    rf <- featureSelect.randomForest(data = data
                                     ,dataFrom=dataFrom
                                     ,save = save
                                     ,feature = feature
                                     ,folder = folder
                                     )
    cox <- featureSelect.UnivariateCox(data = data
                                       ,dataFrom = dataFrom,
                                       feature = feature
                                       ,cutoff = cutoff
                                       ,save = save
                                       ,folder = folder
                                       )
    lasso <- featureSelect.lasso(data = data
                                 ,dataFrom =dataFrom,
                                 feature = feature
                                 ,save = save
                                 ,folder = folder)
    return(list(rf = rf,cox = cox,lasso = lasso))
  }else if(method == "lasso"){
    lasso <- featureSelect.lasso(data = data
                                 ,dataFrom =dataFrom,
                                 feature = feature
                                 ,save = save
                                 ,folder = folder)
    return(lasso)
  }else if(method == "cox"){
    cox <- featureSelect.UnivariateCox(data = data
                                       ,dataFrom = dataFrom,
                                       feature = feature
                                       ,cutoff = cutoff
                                       ,save = save
                                       ,folder = folder
                                       )
    return(cox)
  }else if(method == "randomForest"){
    rf <- featureSelect.randomForest(data = data
                                     ,dataFrom=dataFrom
                                     ,save = save
                                     ,feature = feature
                                     ,folder = folder
    )
    return(rf)
  }else{
    message("method")
  }
}

#' MultivariateCOX.verify
#'
#' @param dataset For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param MulCox.sigFactors For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param coef For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#'
#' @return list
#' @export MultivariateCOX.verify
#'
MultivariateCOX.verify <- function(dataset,MulCox.sigFactors,coef){
  signature <- as.matrix(subset(dataset,select = MulCox.sigFactors)) %*% as.matrix(exp(coef))
  dataset <- mutate(dataset, signature = as.numeric(signature),.before = 1)
  ## KM分析进行验证 ##
  Group <- ifelse(dataset$signature > median(dataset$signature),'High-Score','Low-Score') ## 这里的命名要注意！！！
  dataset <- mutate(dataset, signatureGruop = Group,.before = 2)
  fit <- survfit(Surv(surTime, vitalStat)~signatureGruop,data = dataset)
  p <- ggsurvplot(fit, conf.int=F, pval=T,
                  risk.table=T,
                  legend.labs = c('High-Score','Low-Score'),
                  legend.title='Risk Score',
                  palette = c("dodgerblue2","orchid2"),
                  risk.table.height = 0.3)
  return(list(dataset = dataset
              ,survfit = fit
              ,ggsurvplot = p
              ))
}



#' MultivariateCOX
#'
#' @param data For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param dataFrom For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param feature For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param train_prop For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param cutoff For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param save For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @param folder For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#'
#' @return For more learning materials, please refer to https://github.com/BioInfoCloud/MedBioInfoCloud.
#' @export MultivariateCOX
#'
MultivariateCOX <- function(data
                            ,dataFrom ="mergeSurExp",
                            feature ="all"
                            ,train_prop = 0.8
                            ,cutoff = 0.05
                            ,save = TRUE
                            ,folder = "."){
  if(!is.null(dataFrom) & dataFrom == "mergeSurExp"){
    data <- data[,-1]
  }
  if(!(feature == "all")){
    VarNames <- intersect(colnames(data),feature)
    data <- data[,c("vitalStat","surTime",VarNames)]
  }else{
    VarNames <- colnames(data)[-c(1:2)]
  }
  set.seed(2024)
  if(train_prop > 0 & train_prop < 1){
    data <- dplyr::arrange(data,vitalStat)
    sur0 <- sample(rownames(data)[data$vitalStat == 0]
                   ,round(sum(data$vitalStat == 0) * train_prop)
    )
    sur1 <-  sample(rownames(data)[data$vitalStat == 1]
                    ,round(sum(data$vitalStat == 1) * train_prop)
    )
    trainSet <- data[c(sur0,sur1),]
    testSet <- data[-match(c(sur0,sur1),rownames(data)),]
  }else if(train_prop == 1){
    trainSet = data
    testSet = data
  }

  BaSurv <- survival::Surv(time = trainSet$surTime,event = trainSet$vitalStat)
  fml <-  as.formula(paste0('BaSurv~',paste0(VarNames,collapse = '+')))
  MultiCox <- coxph(fml, data = trainSet) ## 多因素Cox回归
  MultiSum <- summary(MultiCox) ## 提取分析结果

  MultiName <- as.character(VarNames)
  MHR <- round(MultiSum$coefficients[,2],2)
  MPV <- round(MultiSum$coefficients[,5],3)
  MCIL <- round(MultiSum$conf.int[,3],2)
  MCIU <- round(MultiSum$conf.int[,4],2)
  MCI <- paste0(MCIL,'-',MCIU)
  MulCox <- data.frame('characteristics' = MultiName, ## 构建多因素cox分析结果数据框
                       'Hazard Ratio' = MHR,
                       'CI95' = MCI,
                       "HR.95L" = MultiSum$conf.int[,"lower .95"],
                       "HR.95H" = MultiSum$conf.int[,"upper .95"],
                       'P Value' = MPV)
  MulCox.sigFactors <- MulCox$characteristics[which(MulCox$P.Value < cutoff)] %>% as.character() ## 确定最后有生存意义的因子
  if(length(MulCox.sigFactors) >=2){
    fml <-  as.formula(paste0('BaSurv~',paste0(MulCox.sigFactors,collapse = '+'))) ## 通过多因素Cox回归构建一个评分
    MultiCox <- coxph(fml, data = trainSet)
    MultiSum <- summary(MultiCox)

    coef <- MultiSum$coefficients[,1] %>% as.numeric() %>% exp()
    trainSet.result = MultivariateCOX.assist(dataset = trainSet
                                ,MulCox.sigFactors = MulCox.sigFactors
                                ,coef = coef
                                )
    testSet.result <- MultivariateCOX.assist(dataset = testSet
                                    ,MulCox.sigFactors = MulCox.sigFactors
                                    ,coef = coef
    )
    if(save == TRUE){
      save(MulCox,trainSet.result,testSet.result,file = paste0(folder,"/MulCox.result.Rdata"))
      writeLines(MulCox.sigFactors,con = paste0(folder,"/MulCox.selectedFeature.txt"))
      write.csv(MultiSum$coefficients,paste0(folder,'/trainSetFinal_coefficients.csv'))
      write.table(MulCox,file=paste0(folder,"/all.feature.MulCoxResult.txt"),sep="\t",
                  row.names=F,quote=F)
    }
    return(list(MulCox = MulCox
                ,trainSet.result = trainSet.result
                ,testSet.result = testSet.result))
  }
}
