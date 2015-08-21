\name{sfreemap.analyse}
\alias{sfreemap.analyse}
\title{Creates an object suitable analysing results over many phylogenies}
\usage{
sfreemap.analyse(base_tree, trees)
}
\arguments{
    \item{base_tree}{A tree result of a sfreemap.map execution}
    \item{trees}{A 'multiPhylo' object, result of sfreemap.map execution}
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
