\name{add.data}
\alias{add.data}
\title{Adds data from \code{map} into \code{obj}.}
\usage{
add.data(obj, map, name)
}
\arguments{
    \item{obj}{A result of \link{sfreemap.correlation} function;}
    \item{map}{A single phylogenetic tree as a result of \link{sfreemap.map};}
    \item{name}{A name for this mapping;}

}
\description{
    Adds data from \code{map} into \code{obj}.
}
\details{
    This return of this function can later be plotted in a correlation matrix
    using the function \link{plot}.
}
\value{
    \code{obj} with new data added.
}

%\references{
%    References
%}

\author{Diego Pasqualin \email{dpasqualin@inf.ufpr.br}}
\seealso{
    \link{sfreemap.correlation}, \link{plot}
}
%\keyword{phylogenetics}
\keyword{utilities}
