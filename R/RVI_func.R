#' @title Regional Vulnerability Index
#' @description
#' The Regional Vulnerability Index (RVI), a statistical measure of brain structural abnormality, quantifies an individual’s
#' similarity to the expected pattern of deficits seen in schizophrenia derived from large-scale meta-analyses by the ENIGMA
#' consortium. This package outputs the inverse-normal transformed residuals, z-normalized INT residuals, and the final RVI
#' coefficients.
#' @param ID a column name of subject IDs in data.
#' @param DXcontrol an character string specifying control subset(i.e. DXcontrol='DX==0'). Mean and standard deviation for
#' z-normalization should be calculated in healthy controls.
#' @param covariates an optional character vector specifying column names of covariates (i.e. Age, Sex). If covariates=NULL,
#' residuals will not be adjusted for covariates.
#' If covariates are specified(i.e. covariates=c('Age','Sex')), residuals will be adjusted for covariates.
#' @param resp.range a numeric vector specifying column range of responses
#' @param EP a numeric vector specifying an expected pattern of measurements.
#' @param data a data frame contains a column of subject IDs, a column of controls, columns of covariates, columns of responses.
#' @details
#' The RVI is developed as a simple measure of agreement between an individual's pattern of regional neuroimaging traits and the expected
#' pattern of schizophrenia. First, residuals are extracted from simple linear regression by regressing out effects of optional covariates
#' (age, sex) and intracranial brain volume using the full sample (patients and controls). The residuals are inverse-normal transformed (INT)
#' based on residuals' ranks, and then the transformed residuals are z-normalized using mean and standard deviation calculated in healthy
#' controls. The RVI is then calculated as a Pearson correlation coefficient between z-normalized values of subjects and their expected
#' pattern. These expected patterns include cortical, subcortical, and white matter intracranial brain volumes for Schizophrenia Spectrum
#' Disorders (SSD).
#' @return A list with the following elements:
#'   \item{i.norm.resid}{inverse-normal transformed(INT) residuals}
#'   \item{z.norm.resid}{z-normalized INT residuals}
#'   \item{RVI}{pearson correlation coefficient between the z-normalized INT residuals of subjects and their expected pattern}
#' @note
#' The RVI_func() function is developed at the Maryland Psychiatric Research Center, Department of Psychiatry,
#' University of Maryland School of Medicine. This project is supported by NIH R01 EB015611 grant. Please cite our funding if
#' you use this software.
#' @references
#' Kochunov P, Fan F, Ryan MC, et al. Translating ENIGMA schizophrenia findings using the regional vulnerability index: Association
#' with cognition, symptoms, and disease trajectory (2020). Hum Brain Mapp. 2020;10.1002/hbm.25045. \doi{10.1002/hbm.25045}
#'
#' R Core Team (2020). R: A language and environment for statistical computing. R
#' Foundation for Statistical Computing, Vienna, Austria. URL
#' https://www.R-project.org/.
#' @examples
#' E.P <- c(-0.37,0.31,-0.02,-0.08,-0.21,0.46,0.31,0.25)
#' RVI1 <- RVI_func(ID='ID', DXcontrol='DX==0', covariates=c('Age','Sex'), resp.range=c(5:12),
#' EP=E.P, data=RVIpkg::test)
#' RVI2 <- RVI_func(ID='ID', DXcontrol='DX==0', covariates=NULL, resp.range=c(5:12),
#' EP=E.P, data=RVIpkg::test)
#' @export

RVI_func <- function(ID, DXcontrol, covariates, resp.range,EP, data) {
  resp.names <- names(data)[resp.range]
  resp.data <- data.frame(data[,resp.range])
  colnames(resp.data) <- resp.names
  control.idx <-  rownames(subset(data,eval(parse(text=DXcontrol))))

  i.norm.data <- data.frame(ID=data[,which(names(data)==ID)])
  if(is.null(covariates)==T){
    for (i in 1:ncol(resp.data)) {
      formula <- stats::as.formula(paste(resp.names[i], ' ~ ', 1))
      model <- stats::lm(formula, data=data)
      model.residuals <- stats::resid(model)
      i.norm.residuals <- stats::qnorm((rank(model.residuals, na.last = "keep")) / (sum(!is.na(model.residuals))+1))
      i.norm.data <- cbind(i.norm.data,i.norm.residuals)
    }
    colnames(i.norm.data) <- c(ID,resp.names)
  }
  else{
    for (i in 1:ncol(resp.data)) {
      formula <- stats::as.formula(paste(resp.names[i], ' ~ ', paste(covariates, collapse = "+")))
      model <- stats::lm(formula, data=data)
      model.residuals <- stats::resid(model)
      i.norm.residuals <- stats::qnorm((rank(model.residuals, na.last = "keep")) / (sum(!is.na(model.residuals))+1))
      i.norm.data <- cbind(i.norm.data,i.norm.residuals)
    }
    colnames(i.norm.data) <- c(ID,resp.names)
  }

  z.norm.data <- data.frame(ID=data[,which(names(data)==ID)])
  for (i in 2:ncol(i.norm.data)) {
    inorm.control.mean <- mean(i.norm.data[control.idx,i])
    inorm.control.sd <- stats::sd(i.norm.data[control.idx,i])
    z.norm.residuals <- (i.norm.data[,i] - inorm.control.mean)/inorm.control.sd
    z.norm.data <- cbind(z.norm.data,z.norm.residuals)
  }
  colnames(z.norm.data) <- c(ID,resp.names)

  RVI.data <- data.frame(ID=data[,which(names(data)==ID)],RVI=NA)
  for (i in 1:nrow(z.norm.data)) {RVI.data$RVI[i] <- -stats::cor(unlist(z.norm.data[i,2:(length(resp.range)+1)]),EP)}

  out <- list()
  out$i.norm.resid <- i.norm.data
  out$z.norm.resid <- z.norm.data
  out$RVI <- RVI.data

  return (out)
}