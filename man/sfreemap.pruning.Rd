\name{sfreemap.pruning}
\alias{sfreemap.pruning}
\title{Prune a phylogetic tree}
\usage{
sfreemap.pruning(t1, t2, reroot=NULL)
}
\arguments{
    \item{t1}{The tree to be prunned (\code{phylo} object.}
    \item{t2}{The base tree (\code{phylo} object.)}
    \item{reroot}{Optional. The node that should represent the new root of
    \code{t2}}
}
\description{
Returns t1 with only the tips that are in t2 too. Optionally reroot t1.
}
\value{
A \code{phylo} object representing the prunned tree.
}

%\references{
%    References
%}

\author{Diego Pasqualin \email{dpasqualin@inf.ufpr.br}}
%\seealso{
%    See Also
%}
%\keyword{phylogenetics}
\keyword{utilities}
