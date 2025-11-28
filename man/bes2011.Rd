\name{bes2011}
\alias{bes2011}
\docType{data}
\title{Georeferenced responses of Neighborhood Social Capital Scale from BES2011}
\description{
Responses of the Neighborhood Social Capital Scale (NSCS) from Baltimore
Ecosystem Study (BES) Household Telephone Survey Wave 5 (2011), which excluded
the samples from Queen Anne's County. The responses were geocoded and transformed
to a SpatialPointsDataFrame.
}
\usage{data("bes2011")}
\format{
a SpatialPolygonsDataFrame with 1,410 responses on the following 5 items:
\describe{
  \item{\code{willHelp}}{People in the neighborhood are willing to help one another}
  \item{\code{closeKnit}}{This is a close knit neighborhood}
  \item{\code{trust}}{People in this neighborhood can be trusted}
  \item{\code{solComProb}}{There are many opportunities to meet neighbors and work on community problems}
  \item{\code{relGroups}}{Churches, temples and other volunteer groups actively support the neighborhood}
  }
}
\source{
Cary Institute Of Ecosystem Studies, M. Grove, and D. Locke. 2018. BES Household
Telephone Survey ver 180. Environmental Data Initiative.
}
\examples{data("bes2011")}
\keyword{datasets}
