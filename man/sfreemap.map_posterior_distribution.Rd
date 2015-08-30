\name{sfreemap.map_posterior_distribution}
\alias{sfreemap.map_posterior_distribution}
\title{Creates an object suitable analysing results over many phylogenies}
\usage{
sfreemap.map_posterior_distribution(base_tree, trees, scale=TRUE)
}
\arguments{
    \item{base_tree}{A tree result of a sfreemap.map execution}
    \item{trees}{A 'multiPhylo' object, result of sfreemap.map execution}
    \item{scale}{Whether the function should scale the dwelling times to the
    branch length, resulting in a proportion of the time instead of an absolute
    value}
}
\description{
    This functions returns an object that can be used with
    \link{sfreemap.plot_distribution} to plot the distribution dwelling time
    distribution of a state in a particular node.
}
%\details{
%    Details
%}
\value{
    Stuff..
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
