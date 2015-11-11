\name{reorderSfreemap}
\alias{reorderSfreemap}
\title{Reorder a \code{phylo} object modified by sfreemap.map}
\usage{
reorderSfreemap(tree, order='cladewise', index.only = FALSE, ...)
}
\arguments{
    \item{tree}{A tree, result of a \link{sfreemap.map} execution;}
    \item{order}{The resulting order. Can be \code{cladewise} or \code{pruningwise}.}
    \item{index.only}{logical value indicating whether only an index should be returned.}
    \item{...}{other arguments, passed to \link[ape]{reorder.phylo}}
}
\description{
    This function reorders the edges (and all related objects) of an object \code{phylo} modified by \link{sfreemap.map}.

    From \code{\link[ape]{reorder.phylo}} (package ape):

    In the ``cladewise'' order each clade is formed by a series of contiguous rows. In the ``postorder'' order, the rows are arranged so that computations following pruning-like algorithm the tree (or postorder tree traversal) can be done by descending along these rows (conversely, a preorder tree traversal can be performed by moving from the last to the first row). The ``pruningwise'' order is an alternative ``pruning'' order which is actually a bottom-up traversal order (Valiente 2002). (This third choice might be removed in the future as it merely duplicates the second one which is more efficient.) The possible multichotomies and branch lengths are preserved.
}
%\details{
%    Details
%}
\value{
    A phylogenetic tree of class \code{phylo} with its edges reordered.
}

%\references{
%    \code{\link{sfreemap.map}}
%}

\author{Diego Pasqualin \email{dpasqualin@inf.ufpr.br}}
\seealso{
    \code{\link{sfreemap.map}}, \code{\link{reorder.phylo}}
}
\keyword{phylogenetics}
%\keyword{simulation}
