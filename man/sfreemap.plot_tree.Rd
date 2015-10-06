\name{sfreemap.plot_tree}
\alias{sfreemap.plot_tree}
\title{Plot the distribution of states in a particular given node}
\usage{
sfreemap.plot_tree(base_tree, trees, state, type='emr'
					        , conf_level=95, number_of_ticks=20
                            , tip_states=NULL)
}
\examples{
sm <- sfreemap.map(sfreemap.trees, sfreemap.tips, Q='empirical', parallel=FALSE)
result <- sfreemap.plot_tree(sm[[1]], sm, state='1')
}
\arguments{
	\item{base_tree}{trololo}
    \item{trees}{base_tree}
	\item{state}{trololo}
	\item{type}{\code{emr} for expected markov reward (dwelling times) or
	\code{lmt} for labelled makov transistions (number of transitions).}
	\item{conf_level}{trololo}
	\item{number_of_ticks}{trololo}
    \item{tip_states}{tip states as provided to \code{\link{sfreemap.map}}}
}
\description{
    This functions plots something.
}
%\details{
%    Details
%}
\value{
    Besides the plot, this function returns a lot of things
}

\author{Diego Pasqualin \email{dpasqualin@inf.ufpr.br}}
\seealso{
    \code{\link{sfreemap.map_posterior_distribution}}
}
\keyword{phylogenetics}
%\keyword{simulation}
