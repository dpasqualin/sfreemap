\name{sfreemap.map_posterior_distribution}
\alias{sfreemap.map_posterior_distribution}
\title{Creates an object suitable analysing results over many phylogenies}
\usage{
sfreemap.map_posterior_distribution(base_tree, trees, scale=TRUE)
}
\arguments{
    \item{base_tree}{A single tree, result of the function \link{sfreemap.map}}
    \item{trees}{A \code{multiPhylo} object, result of \link{sfreemap.map} execution}
    \item{scale}{Whether the function should scale the dwelling times to the
    branch length, resulting in a proportion of the time instead of an absolute
    value (the expected number of transistions cannot be scaled).}
}
\description{
    This function will compared the \code{trees} with the \code{base_tree} given
    as argument and return an object that can be used with
    \link{sfreemap.plot_distribution} to plot the posterior distribution for
    dwelling times or number of transistions in a particular node.
}
%\details{
%    Details
%}
\value{
    A names list with two names: \code{emr} for expected markov reward (the
    dwelling times for states) and \code{lmt} for the labelled markov
    transistions (the expected number of transitions of states).

    Each criteria (emr and lmt) is composed by an array three dimentions, like
    this: $emr[number of trees, number of states, number of nodes]$.

    This structure represent the corresponding value for each node and state on
    each tree (from \code{trees}), when compared o \code{base_tree}. When a node
    from \code{base_tree} has no match, the corresponding state values will be
    NA (not a number).
}

\author{Diego Pasqualin \email{dpasqualin@inf.ufpr.br}}
\seealso{
    \code{\link{sfreemap.map}}
}
\keyword{phylogenetics}
%\keyword{simulation}
