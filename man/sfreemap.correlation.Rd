\name{sfreemap.correlation}
\alias{sfreemap.correlation}
\title{Creates an object that can be used to plot a correlation matrix}
\usage{
sfreemap.correlation(state)
}
\arguments{
    \item{state}{The character state that will be used to calculate The
    correlation;}
}
\description{
    Creates an object that can be used to plot a correlation matrix using
    \link{add.data} and \link{plot}.
}
\details{
    Correlation matrix is useful to check if different methods implies on
    similar results. To if you run \link{sfreemap.map} with different parameters
    you can easily check the correlation of the outcomes by creating an objects
    with this function and then adding data to it using \link{add.data}.

    Finally, plot the correlation matrix with \link{plot}.
}
\value{
    An object of class \code{sfreemap.correlation}
}

%\references{
%    References
%}

\author{Diego Pasqualin \email{dpasqualin@inf.ufpr.br}}
\seealso{
    \link{sfreemap.map}, \link{add.data}, \link{plot}
}
%\keyword{phylogenetics}
\keyword{utilities}
