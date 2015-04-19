\name{malaria}
\alias{malaria}
\docType{data}
\title{
Gambia Malaria Data
}
\description{
Malaria prevalence in children in Gambia. The data are constructed from \code{gambia} data frame in \code{geoR} package (Diggle and Ribeiro, 2007) by village aggregation. 
}
\usage{data(malaria)}
\format{
A data frame with the 65 observations with the following variables
\tabular{ll}{
\code{x} \tab x-coordinate of the village (UTM). \cr
\code{y} \tab y-coordinate of the village (UTM). \cr
\code{cases} \tab number of sampled children with malaria in each village. \cr
\code{size} \tab number of sampled children in each village. \cr
\code{age} \tab mean age of the sampled children in each village. \cr
\code{netuse} \tab frequency of sampled children who regularly sleep under a bed-net in each village. \cr
\code{treated} \tab frequency of sampled children whose bed-net is treated. \cr
\code{green} \tab satellite-derived measure of the green-ness of vegetation in the immediate vicinity of the village (arbitrary units). \cr
\code{phc} \tab indicator variable denoting the presence (1) or absence (0) of a health center in the village. \cr
\code{area} \tab indicator of the village area (Diggle et al., 2002).\cr
}
}
\source{
Diggle, P.J. and Ribeiro Jr, P.J. (2007). \emph{Model Based Geostatistics}. New York: Springer.
}
\references{
Thomson, M., Connor, S., D Alessandro, U., Rowlingson, B., Diggle, P., Cresswell, M. and Greenwood, B. (1999). Predicting malaria infection in Gambian children from satellite data and bednet use surveys: the importance of spatial correlation in the interpretation of results. \emph{American Journal of Tropical Medicine and Hygiene} \bold{61}, 2--8.

Diggle, P., Moyeed, R., Rowlingson, B. and Thomson, M. (2002). Childhood malaria in The Gambia: a case-study in model-based geostatistics, \emph{Applied Statistics} \bold{51}, 493--506.

}
\examples{
data(malaria)
}
\keyword{datasets}
