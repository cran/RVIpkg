% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RVI_func.R
\name{RVI_func}
\alias{RVI_func}
\title{Regional Vulnerability Index}
\usage{
RVI_func(ID, DXcontrol, covariates = NULL, resp.range, EP, data)
}
\arguments{
\item{ID}{a column name of subject IDs in data.}

\item{DXcontrol}{a character string specifying control subset(i.e. DXcontrol='DX==0'). Mean and standard deviation for
z-normalization should be calculated in healthy controls.}

\item{covariates}{an optional character vector specifying column names of covariates (i.e. Age, Sex). If covariates=NULL (the default),
residuals will not be adjusted for covariates.
If covariates are specified(i.e. covariates=c('Age','Sex')), residuals will be adjusted for covariates.}

\item{resp.range}{a numeric vector specifying column range of regional neuroimaging traits.}

\item{EP}{a numeric vector specifying an expected pattern of measurements. Expected patterns(EP.WM, EP.GM and EP.Subcortical) for
WM, GM and Subcortical are included in the package(Note: If you use an expected pattern from the package, the order of regional neuroimaging
traits need to match the corresponding order of the expected pattern). Check patterns(i.e. RVIpkg::EP.WM$SSD, RVIpkg::EP.WM$MDD, RVIpkg::EP.WM$AD
EP.GM$SSD, RVIpkg::EP.GM$MDD, RVIpkg::EP.GM$AD, EP.Subcortical$SSD, RVIpkg::EP.Subcortical$MDD, RVIpkg::EP.Subcortical$AD .etc.)}

\item{data}{a data frame contains a column of subject IDs, a column of controls, columns of covariates, columns of responses.}
}
\value{
A list with the following elements:
\item{i.norm.resid}{inverse-normal transformed(INT) residuals}
\item{z.norm.resid}{z-normalized INT residuals}
\item{RVI}{RVI.r:pearson correlation coefficient between the z.norm.resid of subjects and their expected pattern. Fisher.Z:Fisher
transformation of RVI.r. Dot.product: dot or scalar product of z.norm.resid and a expected pattern. Z.projection: projection on a pattern.
cos: cosin between Z.projection and z.norm.resid. magnitude: magnitude of z.norm.resid.}
}
\description{
The Regional Vulnerability Index (RVI), a statistical measure of brain structural abnormality, quantifies an individual’s
similarity to the expected pattern of deficits seen in schizophrenia derived from large-scale meta-analyses by the ENIGMA
consortium. This package outputs the inverse-normal transformed residuals, z-normalized INT residuals, and the final RVI
coefficients.
}
\details{
The RVI is developed as a simple measure of agreement between an individual's pattern of regional neuroimaging traits and the expected
pattern of schizophrenia. First, residuals are extracted from simple linear regression by regressing out effects of optional covariates
(age, sex) and intracranial brain volume using the full sample (patients and controls). The residuals are inverse-normal transformed (INT)
based on residuals' ranks, and then the transformed residuals are z-normalized using mean and standard deviation calculated in healthy
controls. The RVI is then calculated as a Pearson correlation coefficient between z-normalized values of subjects and their expected
pattern. These expected patterns include cortical, subcortical, and white matter intracranial brain volumes for Schizophrenia Spectrum
Disorders (SSD).
}
\note{
The RVI_func() function is developed at the Maryland Psychiatric Research Center, Department of Psychiatry,
University of Maryland School of Medicine. This project is supported by NIH R01 EB015611 grant. Please cite our funding if
you use this software.
}
\examples{
E.P <- c(-0.37,0.31,-0.02,-0.08,-0.21,0.46,0.31,0.25)
RVI1 <- RVI_func(ID='ID', DXcontrol='DX==0', covariates=c('Age','Sex'), resp.range=c(5:12),
EP=E.P, data=RVIpkg::data)
RVI2 <- RVI_func(ID='ID', DXcontrol='DX==0', covariates=NULL, resp.range=c(5:12),
EP=E.P, data=RVIpkg::data)
RVI3 <- RVI_func(ID='ID', DXcontrol='DX==0', covariates=c('Age','Sex'), resp.range=c(5:12),
EP=RVIpkg::EP.Subcortical$SSD, data=RVIpkg::data)
}
\references{
Kochunov P, Fan F, Ryan MC, et al. Translating ENIGMA schizophrenia findings using the regional vulnerability index: Association
with cognition, symptoms, and disease trajectory (2020). Hum Brain Mapp. 2020;10.1002/hbm.25045. \doi{10.1002/hbm.25045}

R Core Team (2020). R: A language and environment for statistical computing. R
Foundation for Statistical Computing, Vienna, Austria. URL
https://www.R-project.org/.
}
