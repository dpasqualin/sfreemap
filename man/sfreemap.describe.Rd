\name{sfreemap.describe}
\alias{sfreemap.describe}
\title{describe a \code{phylo} object modified by sfreemap.map}
\usage{
sfreemap.describe(tree)
}
\arguments{
    \item{tree}{A result of a sfreemap.map execution, can be of types
\code{phylo} or \code{multiPhylo}}

}
\description{
    This function returns the expected number of transitions and dwelling times
    of for the states of the tree passed as argument.

    If tree is a \code{multiPhylo} the result will be the mean of the trees.
}
%\details{
%    Details
%}
\value{
    A named list with values \code{transitions} and \code{dwelling_times} for
    each state.
}

%\references{
%    \code{\link{sfreemap.map}}
%}

\author{Diego Pasqualin \email{dpasqualin@inf.ufpr.br}}
\seealso{
    \code{\link{sfreemap.map}}
}
\keyword{utilities}
%\keyword{simulation}
