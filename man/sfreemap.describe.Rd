\name{sfreemap.describe}
\alias{sfreemap.describe}
\title{describe a \code{phylo} object modified by sfreemap.map}
\usage{
sfreemap.describe(tree)
}
\arguments{
    \item{tree}{A result of a sfreemap.map execution (\code{phylo} or \code{multiPhylo})}

}
\description{
    This function returns the mean number of transitions and dwelling times of tree (which can be multiple trees as a multiPhylo object).
}
%\details{
%    Details
%}
\value{
    A named list with values \code{transitions} and \code{dwelling_times}.
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
