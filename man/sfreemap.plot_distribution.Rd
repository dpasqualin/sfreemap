\name{sfreemap.plot_distribution}
\alias{sfreemap.plot_distribution}
\title{Plot the distribution of a state in a particular node}
\usage{
sfreemap.plot_distribution(node, states=NULL, conf_level=90, number_of_ticks=20
							, type='emr', ...)
}
\examples{
sm <- sfreemap.map(sfreemap.trees, sfreemap.tips, Q='empirical', parallel=FALSE)
map <- sfreemap.map_posterior_distribution(sm[[1]], sm)
node <- map$emr[,,120]
result <- sfreemap.plot_distribution(node)
}
\arguments{
	\item{node}{The data for the base node, extracted from \link{sfreemap.map_posterior_distribution}}
    \item{states}{The base state}
	\item{conf_level}{The confidence desired}
	\item{number_of_ticks}{you know it..}
	\item{type}{What are you plotting?}
	\item{...}{Other parameters..}
}
\description{
    This functions plots a distribution for a particular state in
	particular node, highlighting the confidence level desired.
}
%\details{
%    Details
%}
\value{
    Besides the plot, this function returns a list containing the following
    values:
    \itemize{
        \item \code{px} and \code{py}: The x and y axis of the plot
        \item \code{lx} and {ly}: The start and end value of the line discribing
        the confidence level
        \item \code{ticks}: The ticks (cuts) in the possible lengths of the
        ending in the node given as argument
        \item \code{probabilities}: A vector representing the probability of
        finding a dwelling time for the state/node in all trees correspondent
        with the tick in the same position.
    }
}

%\references{
%    \code{\link{sfreemap.map}}
%}

\author{Diego Pasqualin \email{dpasqualin@inf.ufpr.br}}
\seealso{
    \code{\link{sfreemap.map}}
}
\keyword{phylogenetics}
%\keyword{simulation}
