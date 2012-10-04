\name{formatDistData}
\alias{formatDistData}
\title{Bin distance data}
\usage{formatDistData(distData, distCol, transectNameCol, dist.breaks,
                      occasionCol)}
\description{Convert individual-level distance data to the transect-level
format required by \code{\link{distsamp}} or \code{\link{gdistsamp}}}
\arguments{
\item{distData}{data.frame where each row is a detected individual.
Must have at least 2 columns. One for distances and the other for
transect names.}
\item{distCol}{character, name of the column in distData that contains
  the distances. The distances should be numeric.}
\item{transectNameCol}{character, column name containing transect
  names. The transect column should be a factor.}
\item{dist.breaks}{numeric vector of distance interval cutpoints. Length
  must equal J+1.}
\item{occasionCol}{optional character. If transects were visited more
  than once, this can be used to format data for \code{gdistsamp}. It is
  the name of the column in distData that contains the occasion
  numbers. The occasion column should be a factor.}
}
\details{This function creates a site (M) by distance interval (J) response
matrix from a data.frame containing the detection distances for each
individual and the transect names. Alternatively, if each transect was
surveyed T times, the resulting matrix is M x JT, which is the format
required by \code{\link{gdistsamp}}, see\code{\link{unmarkedFrameGDS}}.}

\value{An M x J or M x JT matrix containing the binned distance
  data. Transect  names will become rownames and colnames will describe
  the distance intervals.}

\note{It is important that the factor containing transect names includes
levels for all the transects surveyed, not just those with >=1
detection. Likewise, if transects were visited more than once, the
factor containing the occasion numbers should include levels for all
occasions. See the example for how to add levels to a factor.}
\seealso{\code{\link{distsamp}}, \code{\link{unmarkedFrame}}}
\examples{
# Create a data.frame containing distances of animals detected
# along 4 transects.
dat <- data.frame(transect=gl(4,5, labels=letters[1:4]),
                  distance=rpois(20, 10))
dat

# Look at your transect names.
levels(dat$transect)

# Suppose that you also surveyed a transect named "e" where no animals were
# detected. You must add it to the levels of dat$transect
levels(dat$transect) <- c(levels(dat$transect), "e")
levels(dat$transect)

# Distance cut points defining distance intervals
cp <- c(0, 8, 10, 12, 14, 18)

# Create formated response matrix
yDat <- formatDistData(dat, "distance", "transect", cp)
yDat

# Now you could merge yDat with transect-level covariates and
# then use unmarkedFrameDS to prepare data for distsamp


## Example for data from multiple occasions

dat2 <- data.frame(distance=1:100, site=gl(5, 20),
                   visit=factor(rep(1:4, each=5)))
cutpt <- seq(0, 100, by=25)
y2 <- formatDistData(dat2, "distance", "site", cutpt, "visit")
umf <- unmarkedFrameGDS(y=y2, numPrimary=4, survey="point",
                        dist.breaks=cutpt, unitsIn="m")

}