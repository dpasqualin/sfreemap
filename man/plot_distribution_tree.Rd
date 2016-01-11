\name{plot_distribution_tree}
\alias{plot_distribution_tree}
\title{Plot the dwelling times distribution of a state in a given tree}
\usage{
plot_distribution_tree(map, state, type='emr', conf_level=95, number_of_ticks=20
                    , tip_states=NULL, fsize=0.7, ftype="i", lwd=3)
}
\examples{
sm <- sfreemap(sfreemap.corals.trees, sfreemap.corals.tips, parallel=FALSE)
map <- map_posterior_distribution(sm[[1]], sm, parallel=FALSE)
tree <- plot_distribution_tree(map, state='colonial', tip_states=sfreemap.corals.tips)
}
\arguments{
    \item{map}{A result of a \link{map_posterior_distribution} execution;}
	\item{state}{The state to be plotted}
	\item{type}{\code{emr} for expected markov reward (dwelling times) or
	\code{lmt} for labelled makov transistions (number of transitions). Defaults to \code{emr};}
	\item{conf_level}{An integer representing the confidence level desired, ranging from 0 to 100. Defaults to 95 percent;}
	\item{number_of_ticks}{The number of intervals in which the data will be divided. Defaults to 20;}
    \item{tip_states}{tip states as provided to \code{\link{sfreemap}};}
    \item{fsize}{relative font size for tip labels;}
    \item{ftype}{font type - options are "reg", "i" (italics), "b" (bold), or "bi" (bold-italics);}
    \item{lwd}{line width for three branches.}
}
\description{
    This functions plots the posterior distribution of a state in the tree used as base_tree on \link{map_posterior_distribution}.

    The legend shows the colors representing the incidence of \code{state} given as parameter in the trees analysize in \code{map}. The NA values means that there are no sufficient data to compute the values for the branch given the confidence interval (\code{conf_level} parameter). In other words, the branch is not present in a sufficient number of the analysed trees.
}
%\details{
%    Details
%}
\value{
    Returns the tree plotted, with the data used in the plot stored in \code{tree\$maps}.
}

\author{Diego Pasqualin \email{dpasqualin@inf.ufpr.br}}
\seealso{
    \code{\link{map_posterior_distribution}}
}
\keyword{phylogenetics}
%\keyword{simulation}
