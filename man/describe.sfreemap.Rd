\name{describe.sfreemap}
\alias{describe.sfreemap}
\title{describe a \code{phylo} object modified by sfreemap}
\usage{
describe.sfreemap(tree, ...)
}
\arguments{
    \item{tree}{A result of a sfreemap execution, can be of types \code{phylo} or \code{multiPhylo}}
    \item{...}{Other arguments}
}
\description{
    This function returns the expected number of transitions and dwelling times
    of for the states of the tree passed as argument.

    If tree is a \code{multiPhylo} the result will be the mean value for the trees.
}
%\details{
%    Details
%}
\value{
    A named list with values \code{transitions} and \code{dwelling_times} for
    each state.
}

%\references{
%    \code{\link{sfreemap}}
%}

\author{Diego Pasqualin \email{dpasqualin@inf.ufpr.br}}
\seealso{
    \code{\link{sfreemap}}
}
\keyword{utilities}
%\keyword{simulation}
