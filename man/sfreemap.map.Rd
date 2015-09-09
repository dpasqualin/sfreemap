\name{sfreemap.map}
\alias{sfreemap.map}
\title{Simulation free stochastic character mapping on a phylogenetic tree}
\usage{
sfreemap.map(tree, tip_states, Q='empirical', ...)

}
\arguments{
	\item{tree}{a phylogenetic tree as an object of class \code{"phylo"},
	or a list of trees as an object of class \code{"multiPhylo"}.}

	\item{tip_states}{Two formats are accepted:
		\itemize{
			\item A named vector containing the states of the nodes at the
			tips as values, and the taxa labels as names;
			\item A matrix with taxa labels as rownames, states as colnames and the
			numeric value 1 representing if the taxa is in the corresponding state, or
			0 otherwise.
		}
	}

	\item{Q}{
		The transition rate matrix. Can be given as a matrix or estimated by the
		program, according to the following options:

		\itemize{
			\item \code{"empirical"}: first it fits a continuous-time reversible Markov model for the evolution of x
and then simulates stochastic character histories using that model and the tip states on the tree. This
is the same procedure that is described in Bollback (2006), except that simulation is performed using
a fixed value of the transition matrix, Q, instead of by sampling Q from its posterior distribution (from Phytools).
			\item \code{"mcmc"}: samples Q nsim times from the posterior probability
		distribution of Q using MCMC, then it simulates nsim stochastic maps
		conditioned on each sampled value of Q.
			\item A transition rate matrix. Make sure to name rows and columns properly with the state labels.
			Also, this matrix \bold{must be symmetrical}.
		}

	}

	\item{...}{Optional parameters listed below:

		\code{"prior"}: the prior distribution on the root node of the tree. Options are:
		\itemize{
			\item \code{"equal"} (default): root node is sampled from the conditional
			scaled likelihood distribution at the root;
			\item \code{"estimated"}: the stationary distribution is estimated by numerically
			solving \code{pi*Q=0};
		}

		\code{"tol"} (default: 1e-8): the tolerance for zero elements in Q, elements less then
		tol will be set to zero;

		\code{"parallel"} (default: TRUE): when \code{tree} is of type \code{multiPhylo} we can run
		sfreemap in parallel. The number of processes will be the same as the cores available in your machine.

		When \code{Q="mcmc"} some other parameters might be set:
		\itemize{
			\item \code{"n_simulations"} (default: 100): The number of Q matrices that will be generated;
			\item \code{"burn_in"} (default: 1000): the burn in for the MCMC, when Q="mcmc";
			\item \code{"sample_freq"} (default: 100): sample frequency for the MCMC, when Q="mcmc".
		}
	}


}

\description{
  This function performs an analitic stochastic character mapping on a
	phylogenetic tree (algorithym proposed by Minin and Suchard).
}

%\details{
%	By now only symmetrical model is allowed.
%}

\value{
    Returns a modified object of class "phylo", adding the following data:

		\item{mapped.edge}{a matrix containing the total time spent in each state along each edge of the tree}

		\item{mapped.edge.lmt}{a matrix containing the number of labelled markov transitions for each state along each edge of the tree}

		\item{Q}{the given or estimated value of Q}

		\item{logL}{The likelihood of given or sampled Q}
}

\references{
Vladimir N Minin e Marc A Suchard. Fast, accurate and simulation-free
stochastic mapping. Philosophical Transactions of the Royal Society B:
Biological Sciences, 363 (1512): 3985-3995, 2008.
}

\author{Diego Pasqualin \email{dpasqualin@inf.ufpr.br}}
\seealso{
    See Also
}
\keyword{phylogenetics}
%\keyword{simulation}
