\name{bat_distM}
\alias{bat_distM}
\docType{data}
\title{
  Species pairwise distance matrix for bat species
}
\description{
  This dataset includes a bat species pairwise symmetric distance matrix for 16 species listed in the \code{bat_incidence_data} dataset. Each element in the matrix represents a Gower distance between two species computed from species traits. \cr
}
\usage{data("bat_distM")}
\format{
  \code{bat_distM} is a 16 x 16 (species by species) symmetric data.frame. Each element of the data.frame is between zero (for species pairs with identical traits) and one.  \cr
}

\keyword{datasets}
