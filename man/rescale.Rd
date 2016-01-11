\name{rescale}
\alias{rescale}
\title{Rescale branches of a \code{phylo} or \code{multiPhylo} object to a given height}
\usage{
rescale(tree, height=NULL, parallel=FALSE)
}
\arguments{
    \item{tree}{One or more trees, result of a \link{sfreemap.map} execution;}
    \item{height}{The desired height to scale branches to.}
    \item{parallel}{Whether computation should run in parallel. Only useful for \code{multiPhylo}}
}
\description{
    Rescale a \code{phylo} or \code{multiPhylo} object to a given height.
}
%\details{
%    Details
%}
\value{
    A phylogenetic tree of class \code{phylo} or \code{multiPhylo} with its branches
    rescaled.
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
