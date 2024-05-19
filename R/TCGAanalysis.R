
#' In the expression matrix of TCGA tumor samples, the data of duplicate patients were deleted.
#'
#' @param data  Gene expression matrix of TCGA tumor samples.
#' @param col_rename TRUE or FALSE
#'
#' @return data.frame
#' @export
#'
#' @examples
del_dup_sample <- function(data,col_rename =T){
  data <- data[,sort(colnames(data))]
  pid = gsub("(.*?)-(.*?)-(.*?)-.*","\\1-\\2-\\3",colnames(data))
  reps = pid[duplicated(pid)]##重复样本
  if(length(reps)== 0 ){
    message("There were no duplicate patient samples for this data")
  }else{
    data <- data[,!duplicated(pid)]
  }
  if(col_rename == T){
    colnames(data) <- gsub("(.*?)-(.*?)-(.*?)-.*","\\1-\\2-\\3",colnames(data))
  }
  return(data)
}
