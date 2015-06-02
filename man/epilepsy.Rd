\name{epilepsy}
\alias{epilepsy}
\docType{data}
\title{
Epilitic Seizures Data
}
\description{
Longitudinal study on epilitic seizures (Thall and Vail, 1990; Diggle et al. 2002). The data consist into 59 individuals with five observations each: The baseline eight-week interval and measurements collected at subsequent visits every two-week.
}
\usage{data(epilepsy)}
\format{
\tabular{ll}{
\code{id} \tab patient's id . \cr
\code{age} \tab patient's age. \cr
\code{trt} \tab indicator if the patient is treated with progabide (\code{1}) or with placebo (\code{2}). \cr
\code{counts} \tab number of epileptic seizures. \cr
\code{time} \tab observation period in weeks (\code{8} for baseline and \code{2} for subsequent visits). \cr
\code{visit} \tab indicator if observation at baseline (\code{0}) or subsequent visit (\code{1}). \cr
}
}
\source{
Thall, P.F. and Vail S.C. (1990). Some covariance models for longitudinal count data with overdispersion. \emph{Biometrics} \bold{46}, 657--671.
}
\references{
Diggle, P.J., Heagerty, P., Liang, K.Y. and Zeger, S.L. (2002). \emph{Analysis of Longitudinal Data}. Oxford: Oxford University Press. Second edition.
}