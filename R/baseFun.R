
#' outputGmtFile
#'
#' @param input list or data.frame
#' @param description NA
#' @param filename filepath
#'
#' @return
#' @export
#'
#' @examples
outputGmtFile <- function(input,description = NA,filename){
  if(is.data.frame(input)){
    nms <- unique(input[,1])
    ls <- lapply(nms, function(nm){
      input[input[,1]== nms[nm],2]
    })
    names(ls) <- nms
    input = ls
  }
  if(is.list(input)){
    output <- file(filename, open="wt")
    nms <- names(input)
    lapply(1:length(nms),function(x){
      if(!is.na(description)){
        dsc = description[x]
      }else{dsc = description}
      outlines = paste0(c(nms[x], dsc, input[[x]]),collapse='\t')
      writeLines(outlines, con=output)
    })
    close(output)
  }
}

