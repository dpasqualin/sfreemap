\name{sfreemap.pie_plot}
\alias{sfreemap.pie_plot}
\title{Plot a tree with a summary of time spent no each state by each node}
\usage{
sfreemap.pie_plot(tree, percent_only=FALSE)
}
\arguments{
    \item{tree}{A result of a simmap.map execution.}
    \item{percent_only}{(default FALSE) if TRUE just return the percentage of time spent by each state in the internal nodes}

}
\description{
    Plot a tree with a pizza chart on each node representing the percentage of time spent by each state on each node.

    When \code{"percent_only"} is set to TRUE, just return the calculate percentages, without actually plotting the tree.
}
%\details{
%    Details
%}
\value{
    Return nothing if \code{percentage_only=TRUE} or a matrix with each row being an internal node, columns representing the states and values representing the percentage of time that each state spent on each node.
}

%\references{
%    References
%}

\author{Diego Pasqualin \email{dpasqualin@inf.ufpr.br}}

\seealso{
    \code{\link{sfreemap.map}}
}
%\keyword{phylogenetics}
\keyword{utilities}
