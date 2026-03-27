#' @title Regional Vulnerability Index by Multi-site
#' @description
#' The Regional Vulnerability Index (RVI), a statistical measure of brain structural abnormality, quantifies an individual’s similarity to the expected pattern (effect size)
#' of deficits seen in schizophrenia derived from large-scale meta-analyses by the ENIGMA consortium. This version of RVI function will specifically calculate RVI for data
#' collected from different scanners/sites. The function outputs the inverse-normal transformed (INT) residuals, z-normalized INT residuals, RVI and Alignment Vulnerability Index (AVI).
#' @param ID a column name of subject IDs in data.
#' @param DXcontrol a character string specifying control subset(i.e. DXcontrol='DX==0'or DXcontrol='DX=="CN"'). Mean and standard deviation of z-normalization
#' should be calculated in healthy controls.
#' @param covariates a character vector specifying column names of covariates (i.e. Age, Sex). If covariates=NULL (the default), residuals will not be
#' adjusted for any covariate. If covariates are specified (i.e. covariates=c('Age','Sex')), residuals will be adjusted for covariates.
#' @param batch a character string of a column name indicating the scanning site or scanner or etc.
#' @param resp.range a numeric vector specifying column indices of regional neuroimaging traits.
#' @param EP a numeric vector specifying an expected pattern of regional neuroimaging traits. The expected patterns(EP.WM, EP.GM and EP.Subcortical) for white matter
#' fractional anisotropy (FA), cortical matter thickness and subcortical volume are included in the package (Note: If you use an expected pattern, you need to make
#' sure the order of regional neuroimaging traits in your data match up the corresponding order of the expected pattern). The patterns can be extract in the package
#' (i.e. RVIpkg::EP.WM$SSD, RVIpkg::EP.WM$MDD, RVIpkg::EP.WM$AD, RVIpkg::EP.WM$BD ,RVIpkg::EP.WM$PD .etc.). They were developed using neuroimaging data of UK Biobank (UKBB).
#' @param sign a logical value indicating whether the AVI use signs from RVI.
#' @param fisherZ a logical value indicating whether the result should generate fisher-z transformed RVI.
#' @param data a data frame contains a column of subject IDs, a column of controls, columns of covariates, columns of responses.
#' @details
#' The RVI is developed as a simple measure of agreement between an individual's pattern of regional neuroimaging traits and the expected pattern of schizophrenia.
#' First, all observations of each regional neuroimaging trait are regressed out optional covariates using linear regression, and then residuals are extracted from
#' the model after removing effects of the optional covariates. The optional covariates could be age, sex, intracranial brain volume and/or .etc within each site and then
#' the residuals are inverse-normal transformed (INT) based on residuals' ranks. All INT residuals data from all sites are combined and z-normalized/standardized using mean
#' and standard deviation from the healthy controls to get z-normalized residuals. For each subject, the RVI is then calculated as a Pearson correlation coefficient between the
#' z-normalized INT residuals of the traits and corresponding expected pattern of the traits and the AVI is the dot product of the z-normalized INT residuals of the
#' traits and corresponding expected pattern of the traits. These expected patterns include cortical thickness, subcortical volume, and white matter FA for mental
#' illnesses and metabolic diseases.
#' @return A list with the following elements:
#'   \item{i.norm.resid}{inverse-normal transformed residuals}
#'   \item{z.norm.resid}{z-normalized INT residuals}
#'   \item{RVI}{RVI: the Pearson correlation coefficient between the z-normalized INT residuals and corresponding expected pattern;
#'   AVI: the dot product of the z-normalized INT residuals and corresponding expected pattern;
#'   RVI.fisherz: Fisher z-transformed RVI
#'   }
#' @note
#' The RVI_func() function is developed at the the University of Texas Health Science Center at Houston, Department of Psychiatry and Behavioral Sciences,
#' McGovern Medical School. This project is supported by NIH R01 EB015611 grant. Please cite our funding if you use this software.
#' @references
#' Kochunov P, Fan F, Ryan MC, et al. Translating ENIGMA schizophrenia findings using the regional vulnerability index: Association
#' with cognition, symptoms, and disease trajectory (2020). Hum Brain Mapp. 2020;10.1002/hbm.25045. \doi{10.1002/hbm.25045}
#'
#' R Core Team (2020). R: A language and environment for statistical computing. R
#' Foundation for Statistical Computing, Vienna, Austria. URL
#' https://www.R-project.org/.
#' @examples
#' EP1 <- c(-0.37,0.31,-0.02,-0.08,-0.21,0.46,0.31,0.25)
#' RVI1 <- RVI_batch_func(ID='ID', DXcontrol='DX==0',batch='Site', covariates=c('Age','Sex'),
#' resp.range=c(6:13), EP=EP1, data=RVIpkg::data)
#' @export

RVI_batch_func <- function(ID, DXcontrol, covariates=NULL, batch,resp.range, EP, sign=FALSE, fisherZ=FALSE, data) {
  message("Note: If the old atlas from ENIGMA DTI pipeline was used (e.g. label TAP is missing), then the following columns in your data need to be renamed: IFO->UNC and UNC->TAP.")
  resp.names <- names(data)[resp.range]
  resp.data <- data.frame(data[,resp.range])
  colnames(resp.data) <- resp.names
  control.idx <-  rownames(subset(data,eval(parse(text=DXcontrol))))
  unibatch <- unique(data[,which(names(data)==batch)])
  i.norm.data <- data.frame()
  for (j in 1:length(unibatch)) {
    batch.j <- unibatch[j]
    batch.condition <- paste(batch, "==", batch.j)
    data.batch <- subset(data, eval(parse(text = batch.condition)))
    i.norm.data.batch <- data.frame(ID=data.batch[,which(names(data.batch)==ID)])
    if(is.null(covariates)==T){
      for (i in 1:ncol(resp.data)) {
        formula <- stats::as.formula(paste(resp.names[i], ' ~ ', 1))
        model <- stats::lm(formula, na.action=stats::na.exclude, data=data.batch)
        model.residuals <- stats::resid(model)
        i.norm.residuals <- stats::qnorm((rank(model.residuals, na.last = "keep")) / (sum(!is.na(model.residuals))+1))
        i.norm.data.batch <- cbind(i.norm.data.batch,i.norm.residuals)
      }
      colnames(i.norm.data.batch) <- c(ID,resp.names)
    }
    else{
      for (i in 1:ncol(resp.data)) {
        formula <- stats::as.formula(paste(resp.names[i], ' ~ ', paste(covariates, collapse = "+")))
        model <- stats::lm(formula,na.action=stats::na.exclude, data=data.batch)
        model.residuals <- stats::resid(model)
        i.norm.residuals <- stats::qnorm((rank(model.residuals, na.last = "keep")) / (sum(!is.na(model.residuals))+1))
        i.norm.data.batch <- cbind(i.norm.data.batch,i.norm.residuals)
      }
      colnames(i.norm.data.batch) <- c(ID,resp.names)
    }
    i.norm.data <- rbind(i.norm.data,i.norm.data.batch)
  }
  i.norm.data <- merge(data[,which(names(data)==ID),drop=F],i.norm.data,by=1,all.x = T)

  z.norm.data <- data.frame(ID=data[,which(names(data)==ID)])
  for (i in 2:ncol(i.norm.data)) {
    inorm.control.mean <- mean(i.norm.data[control.idx,i], na.rm = T)
    inorm.control.sd <- stats::sd(i.norm.data[control.idx,i], na.rm = T)
    z.norm.residuals <- (i.norm.data[,i] - inorm.control.mean)/inorm.control.sd
    z.norm.data <- cbind(z.norm.data,z.norm.residuals)
  }
  colnames(z.norm.data) <- c(ID,resp.names)

  RVI.data <- data.frame(ID=data[,which(names(data)==ID)], RVI=NA, AVI=NA)
  sumna <- function(x) {if(all(is.na(x))==T) NA else sum(x, na.rm = TRUE)}
  for (i in 1:nrow(z.norm.data)) {
    RVI.data$RVI[i] <- stats::cor(unlist(z.norm.data[i,2:(length(resp.range)+1)]), EP, use = "na.or.complete")
  }

  if(sign==TRUE){
    for (i in 1:nrow(z.norm.data)) {RVI.data$AVI[i] <- abs(sumna(unlist(z.norm.data[i,2:(length(resp.range)+1)])*EP))*sign(RVI.data$RVI[i])}
  }
  else{
    for (i in 1:nrow(z.norm.data)) {RVI.data$AVI[i] <- sumna(unlist(z.norm.data[i,2:(length(resp.range)+1)])*EP)}
  }

  if(fisherZ==TRUE){
    RVI.data$RVI.fisherZ <- NA
    for (i in 1:nrow(z.norm.data)) {RVI.data$RVI.fisherZ[i] <- 0.5*log((1+RVI.data$RVI[i])/(1-RVI.data$RVI[i]))}
  }

  out <- list()
  out$i.norm.resid <- i.norm.data
  out$z.norm.resid <- z.norm.data
  out$RVI <- RVI.data

  return (out)
}
