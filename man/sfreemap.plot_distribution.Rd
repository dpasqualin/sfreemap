\name{sfreemap.plot_distribution}
\alias{sfreemap.plot_distribution}
\title{Plot the distribution of a state in a particular node}
\usage{
sfreemap.plot_distribution(node, state, conf.level=90, ...)
}
\examples{
sm <- sfreemap.map(sfreemap.trees, sfreemap.tips, Q='empirical', parallel=FALSE)
analysis <- sfreemap.analyse(sm[[1]], sm)
node <- analysis$n120
result <- sfreemap.plot_distribution(node, 1, conf.level=90)
}
\arguments{
	\item{node}{The data for the base node, extracted from \link{sfreemap.analyse}}
    \item{state}{The base state}
	\item{conf.level}{The confidence desired}
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
