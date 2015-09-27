\name{sfreemap.plot_distribution}
\alias{sfreemap.plot_distribution}
\title{Plot the distribution of states in a particular given node}
\usage{
sfreemap.plot_distribution(node, states=NULL, conf_level=95, number_of_ticks=20
							, type='emr', ...)
}
\examples{
sm <- sfreemap.map(sfreemap.trees, sfreemap.tips, Q='empirical', parallel=FALSE)
map <- sfreemap.map_posterior_distribution(sm[[1]], sm)
node <- map$emr[,,120]
result <- sfreemap.plot_distribution(node)
}
\arguments{
	\item{node}{The data for the node to analyse, extracted from \link{sfreemap.map_posterior_distribution}}
    \item{states}{The states to be shown in the plot, defaults to all states.}
	\item{conf_level}{The confidence level desired}
	\item{number_of_ticks}{The number of intervals in the x axis.}
	\item{type}{\code{emr} for expected markov reward (dwelling times) or
	\code{lmt} for labelled makov transistions (number of transitions).}
	\item{...}{nothing special so far.}
}
\description{
    This functions plots the posterior distribution for a state (or all states)
	in particular node, highlighting the confidence level desired.
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

\author{Diego Pasqualin \email{dpasqualin@inf.ufpr.br}}
\seealso{
    \code{\link{sfreemap.map_posterior_distribution}}
}
\keyword{phylogenetics}
%\keyword{simulation}
