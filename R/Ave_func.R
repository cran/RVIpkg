#' @title Optimizing data from UK Biobank
#' @description
#' The Ave_func() can optimize data from UK Biobank(UKB). It will rename field IDs of regional neuroimaging traits to abbreviation names, and then average data of left and right
#' hemispheres of the same field.
#' @param resp.range a numeric vector specifying column range of regional neuroimaging traits.
#' @param type a character string specifying data types of regional neuroimaging traits(i.e. All traits(type='all'), White matter(type='WM'),Gray matter(type='GM') or Subcortical(type='Subcortical'))
#' @param data a data frame contains regional neuroimaging traits with field IDs from UKBB. Default(type='all')
#' @return a dataframe of regional neuroimaging traits with abbreviated field names.
#' @note
#' The Ave_func() function is developed at the Maryland Psychiatric Research Center, Department of Psychiatry,
#' University of Maryland School of Medicine. This project is supported by NIH R01 EB015611 grant. Please cite our funding if
#' you use this software.
#' @references
#' Kochunov P, Fan F, Ryan MC, et al. Translating ENIGMA schizophrenia findings using the regional vulnerability index: Association
#' with cognition, symptoms, and disease trajectory (2020). Hum Brain Mapp. 2020;10.1002/hbm.25045. \doi{10.1002/hbm.25045}
#'
#' R Core Team (2020). R: A language and environment for statistical computing. R
#' Foundation for Statistical Computing, Vienna, Austria. URL
#' https://www.R-project.org/.
#' @export

Ave_func <- function(resp.range, type='all', data) {
  other.data <- data[,-resp.range,drop=F]
  resp.data <- data[,resp.range,drop=F]
  out.data <- resp.data
  for (i in 1:ncol(resp.data)) {
    num_colname <- substr(as.numeric(unlist(regmatches(colnames(resp.data)[i], gregexpr("[[:digit:]]+\\.*[[:digit:]]*", colnames(resp.data)[i])))), start = 1, stop = 5)[1]
    for (j in 1:nrow(FID_list)) {
      FID <- FID_list$ID[j]
      if(is.na(num_colname)==F){
        if(num_colname==FID){
          colnames(out.data)[i] <- FID_list$name[j]
        }
      }
    }
  }

  average.data <- as.data.frame(sapply(split.default(out.data, names(out.data)), rowMeans))

  if(type=='all'){
    average.data$CC <- rowMeans(average.data[,c('GCC','BCC','SCC')], na.rm = FALSE, dims = 1)
    average.data$CR <- rowMeans(average.data[,c('ACR','SCR','PCR')], na.rm = FALSE, dims = 1)
    average.data$IC <- rowMeans(average.data[,c('ALIC','PLIC','RLIC')], na.rm = FALSE, dims = 1)
    average.data1 <- average.data[,EP.WM$WM]
    average.data2 <- average.data[,EP.GM$GM]
    average.data3 <- average.data[,EP.Subcortical$Subcortical]
    average.order <- cbind(average.data1,average.data2,average.data3)
  }
  else if(type=='WM'){
    average.data$CC <- rowMeans(average.data[,c('GCC','BCC','SCC')], na.rm = FALSE, dims = 1)
    average.data$CR <- rowMeans(average.data[,c('ACR','SCR','PCR')], na.rm = FALSE, dims = 1)
    average.data$IC <- rowMeans(average.data[,c('ALIC','PLIC','RLIC')], na.rm = FALSE, dims = 1)
    average.order <- average.data[,EP.WM$WM]
    }
  else if(type=='GM'){average.order <- average.data[,EP.GM$GM]}
  else if(type=='Subcortical'){average.order <- average.data[,EP.Subcortical$Subcortical]}

  out.combine <- cbind(other.data,average.order)
  return (out.combine)
}
