% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Utilities.R
\name{MaxMoles}
\alias{MaxMoles}
\title{Maximum limit of moles (program limit)}
\usage{
MaxMoles(x)
}
\arguments{
\item{x}{component; factor (\code{"C", "H", "O", "N", "S", "P", "M", "E", "Z"})}
}
\value{
numeric; value of maximum limit of moles
}
\description{
Returns the maximum limit of moles used by the
\code{\link{MFAssign}} function.
}
\details{
This is an internal fuction that will not work outside the \code{\link{MFAssign}} function
environment.
}
\examples{
MaxMoles("C")

MaxMoles("C")

}
