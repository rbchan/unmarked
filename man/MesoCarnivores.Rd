\name{MesoCarnivores}
\alias{MesoCarnivores}
\docType{data}
\title{
  Occupancy data for coyote, red fox, and bobcat
}
\description{
  Occupancy data and site covariates for coyote, red fox, and bobcat from 1437 camera trap sites sampled 3 times. Each sampling period represents one week. This data is a simplified form of the dataset used by Rota et al. (2016).
}

\usage{data(MesoCarnivores)}

\format{
  A list with four elements:
  \describe{
    \item{\code{bobcat}}{A 1437x3 occupancy matrix for bobcat}
    \item{\code{coyote}}{A 1437x3 occupancy matrix for coyote}
    \item{\code{redfox}}{A 1437x3 occupancy matrix for red fox}
    \item{\code{sitecovs}}{A data frame containing covariates for the 1437 sites, with the following columns:
      \describe{
        \item{\code{Dist_5km}}{Proportion of disturbed land in 5 km radius}
        \item{\code{HDens_5km}}{Housing density in 5 km radius}
        \item{\code{Latitude}}{Latitude / 100}
        \item{\code{Longitude}}{Longitude / 100}
        \item{\code{People_site}}{Number of photos of people at site / 1000}
        \item{\code{Trail}}{1 if camera was on trail, 0 if not}
      }
    }
  }
}

\source{
  Used with permission of Roland Kays and Arielle Parsons at North Carolina State University and the North Carolina Museum of Natural Sciences.
}

\references{
Rota, C.T., et al. 2016. A multi-species occupancy model for two or more
  interacting species. Methods in Ecology and Evolution 7: 1164-1173.
}

\keyword{datasets}
