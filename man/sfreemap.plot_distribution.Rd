\name{sfreemap.plot_distribution}
\alias{sfreemap.plot_distribution}
\title{Plot the distribution of states in a given node}
\usage{
sfreemap.plot_distribution(map, nodes=NULL, trees=NULL, states=NULL, conf_level=95
	                        , number_of_ticks=20, type='emr')
}
\examples{
sm <- sfreemap.map(sfreemap.corals.trees, sfreemap.corals.tips, parallel=FALSE)
map <- map_posterior_distribution(sm[[1]], sm, parallel=FALSE)
p <- sfreemap.plot_distribution(map, 160)
print(p)
}
\arguments{
	\item{map}{A result of a \link{map_posterior_distribution} execution;}
    \item{nodes}{A vector containing the nodes to filter, or NULL to not filter;}
    \item{trees}{A vector containing the trees to filter, or NULL to not filter;}
    \item{states}{A vector containing the states to filter, or NULL to not filter;}
	\item{conf_level}{An integer representing the confidence level desired, ranging from 0 to 100. Defaults to 95 percent;}
	\item{number_of_ticks}{The number of intervals in the x axis. Defaults to 20;}
	\item{type}{\code{emr} for expected markov reward (dwelling times) or
	\code{lmt} for labelled markov transitions (number of transitions). Defaults to \code{emr};}
}
\description{
    This function plots the posterior distribution for a state (or all states)
	in particular node, highlighting the confidence interval desired calculated using the highest posterior density (HPD).
}
%\details{
%    Details
%}
\value{
    This function returns a \link[ggplot2]{ggplot} object that shows the plot when printed. The print itself can be changed as the pleased, using regular ggplot layers functions.
}

\author{Diego Pasqualin \email{dpasqualin@inf.ufpr.br}}
\seealso{
    \code{\link{map_posterior_distribution}},
    \code{\link{sfreemap.map}}
}
\keyword{phylogenetics}
%\keyword{simulation}
