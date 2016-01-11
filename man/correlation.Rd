\name{correlation}
\alias{correlation}
\title{Creates an object that can be used to plot a correlation matrix}
\usage{
correlation(map, state, name)
}
\examples{
r1 <- sfreemap(sfreemap.corals.trees[[1]], sfreemap.corals.tips, parallel=FALSE)
r2 <- sfreemap(sfreemap.corals.trees[[1]], sfreemap.corals.tips, method='mcmc',
                   n_simulation=1, parallel=FALSE)

cor <- correlation(r1, 'colonial', 'empirical') + correlation(r2, 'colonial', 'mcmc')
plot(cor)
}
\arguments{
    \item{map}{A \code{sfreemap} object, result of \link{sfreemap}
     execution;}
    \item{state}{The character state that will be used to calculate The
    correlation;}
    \item{name}{A string to uniquely identify this mapping.}
}
\description{
    Creates an object that can be used to plot a correlation matrix using
    \link{plot.correlation}.
}
\details{
    Correlation matrix is useful to check if different methods implies on
    similar results. To if you run \link{sfreemap} with different parameters
    you can easily check the correlation of the outcomes by summing up objects
    created with this function and plottig it using \link{plot.correlation}.

    This function returns an object that can handle the plus operator.
    See examples.
}
\value{
    An object of class \code{correlation}
}

%\references{
%    References
%}

\author{Diego Pasqualin \email{dpasqualin@inf.ufpr.br}}
\seealso{
    \link{sfreemap}, \link{plot.correlation}
}
%\keyword{phylogenetics}
\keyword{utilities}
