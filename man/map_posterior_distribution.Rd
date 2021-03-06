\name{map_posterior_distribution}
\alias{map_posterior_distribution}
\title{Analyses a set of trees regarding to a specific tree and output an object that can be used by many other functions of this package}
\usage{
map_posterior_distribution(base_tree, trees, scale.branches=TRUE
                                    , scale.trees=FALSE, parallel=TRUE)
}
\examples{
    sm <- sfreemap(sfreemap.corals.trees, sfreemap.corals.tips, parallel=FALSE)
    map <- map_posterior_distribution(sm[[1]], sm, parallel=FALSE)
    p <- plot_distribution_chart(map, 160)
    print(p)
}
\arguments{
    \item{base_tree}{A single tree, result of the function \link{sfreemap};}
    \item{trees}{A \code{multiPhylo} object, result of \link{sfreemap} execution}
    \item{scale.branches}{Whether the function should scale the dwelling times to the
    branch length, resulting in a proportion of the time instead of an absolute
    value (the expected number of transitions is not scaled).}
    \item{scale.trees}{A value representing the maximum height to scale the tree, or FALSE if the trees should not be scaled.}
    \item{parallel}{Whether the function should run in parallel (defaults to TRUE)}
}
\description{
    This function will compared the \code{trees} against the \code{base_tree} and return an object that can be used with
    \link{plot_distribution_chart} to plot the posterior distribution for
    dwelling times and number of transistions on a particular node and \link{plot_distribution_tree}, which plots the posterior distribution in the entire tree.
}
%\details{
%    Details
%}
\value{
    A named list with three items:

    \itemize{
        \item \code{emr} for expected markov reward (the
        dwelling times for states) and \code{lmt} for the labelled markov
        transitions (the expected number of state transitions).

        Each criteria (emr and lmt) is composed by an array three dimentions, like this: $emr[trees, states, nodes]$.

        This structure represent the corresponding value for each node and state on each tree in \code{trees} when compared o \code{base_tree}. When a node from \code{base_tree} has no match, the corresponding state values will be set to NA.

        \item \code{base_tree} the original base_tree given as argument. This is useful for the functions that use the result of this function, cited in the description.
    }
}

\author{Diego Pasqualin \email{dpasqualin@inf.ufpr.br}}
\seealso{
    \code{\link{sfreemap}}
}
\keyword{phylogenetics}
%\keyword{simulation}
