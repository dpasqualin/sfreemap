\name{sfreemap.map}
\alias{sfreemap.map}
\title{Simulation free stochastic character mapping on a phylogenetic tree}
\usage{
sfreemap.map(tree, tip_states, Q=NULL, type="standard", model="SYM"
                 , method="empirical", ...)

}
\arguments{
	\item{tree}{a phylogenetic tree as an object of class \code{"phylo"},
	or a list of trees as an object of class \code{"multiPhylo"}.}

	\item{tip_states}{Two formats are accepted:
		\itemize{
			\item \bold{When type='standard'}: A named vector containing the states of the nodes at the
			tips as values, and the taxa labels as names or a matrix with taxa
            labels as rownames, states as colnames and the numeric value 1
            representing if the taxa is in the corresponding state, or 0
            otherwise.

            \item \bold{When type='dna'}: A matrix with caracters as columns and
            taxa as rows.
		}
	}

	\item{Q}{
		The transition rate matrix. Can be given as a matrix with state names on dimensons or estimated by the
		program. Options for estimation depend on argumentos \code{method}, \code{model} and \code{type}.
    }

    \item{type}{
        The type of the \code{tip_states} being analysed. It can be "standard", usually used for morphological characters, or "dna", for nucleotideos. Default to "standard".
    }

    \item{model}{
    By choosing to estimate Q user can then select a model.

    When \code{type="standard"} the available methods are:
    \itemize{
        \item \code{"SYM" (default)}: symmetrical model, e.g, \code{matric(c(0,1,2,1,0,3,2,3,0))};
        \item \code{"ER"}: equal rates model, for example \code{matrix(c(0,1,1,0))};
        \item \code{"ARD"}: all rates different model, for example \code{matrix(c(0,1,2,0))};
    }

    When \code{type="dna"} the available methods are JC, F81, K80, HKY, TrNe, TrN, TPM1, K81, TPM1u, TPM2, TPM2u, TPM3, TPM3u, TIM1e, TIM1, TIM2e, TIM2, TIM3e, TIM3, TVMe, TVM, SYM and GTR.

    }

    \item{method}{
        The \code{method} argument is only used when \code{type='standard'} and the available options are:

		\itemize{
			\item \code{empirical (default)}: first it fits a continuous-time reversible Markov model for the evolution of x
and then simulates stochastic character histories using that model and the tip states on the tree. This
is the same procedure that is described in Bollback (2006), except that simulation is performed using
a fixed value of the transition matrix, Q, instead of by sampling Q from its posterior distribution (\link[phytools]{phytools}).
			\item \code{"mcmc"}: samples n_simulation Q matrices from the
            posterior probability distribution of Q using MCMC, then
            performs stochastic maps conditioned on each sampled value of Q.
		}
	}

	\item{...}{Optional parameters, listed below:

    \itemize{
		\item \code{"prior"}: the prior distribution on the root node of the tree. Options are:
		\itemize{
			\item \code{"equal"} (default): root node is sampled from the conditional
			scaled likelihood distribution at the root;
			\item \code{"estimated"}: the stationary distribution is estimated by numerically
			solving \code{pi*Q=0};
		}

		\item \code{"tol"} (default: 1e-8): the tolerance for zero elements in Q, elements less then
		tol will be set to tol;

		\item \code{"parallel"} (default: TRUE): when \code{tree} is of type \code{multiPhylo} we can run
		sfreemap in parallel. The number of processes created will be the same as the cores available in your machine.

		\item When \code{Q="mcmc"} some other parameters might be set:
		\itemize{
			\item \code{"n_simulations"} (default: 100): The number of Q matrices that will be generated;
			\item \code{"burn_in"} (default: 1000): the burn in for the MCMC;
			\item \code{"sample_freq"} (default: 100): number of generations for each sample taken;
		}
    }
    }

}

\description{
  This function performs an analitic stochastic character mapping on a
	phylogenetic tree (algorithym proposed by Minin and Suchard).

 It can be called with a combination of parameters, much like any vectorized
 function in R. In other words, calling it with \code{N} trees (\code{multiPhylo}
 object) and a single rate matrix \code{Q} will return \code{N} mapped trees.
 Calling \code{sfreemap.map} with a single tree (\code{phylo} object) and \code{M}
 \code{Q} matrices will result in \code{M} mapped trees, replicas of the single
 tree with the algorithm applied to it using all \code{Q} matrices. Same logic
 applies to \code{prior} as well. It is important to note though that it you
 pass on \code{N > 1} trees and \code{M > 1} rate matrices (or priors), \code{M}
  and \code{N} should be equal.
}

%\details{
%	By now only symmetrical model is allowed.
%}

\value{
    Returns a modified object of class "phylo", adding the following data:

		\item{mapped.edge}{a matrix containing the expected value for dwelling times for each state along each edge of the tree;}

		\item{mapped.edge.lmt}{a matrix containing the expected number of labelled markov transitions for each state along each edge of the tree;}

		\item{Q}{the given or estimated value of Q;}

		\item{logL}{The likelihood of calculated for the given or sampled Q;}

        \item{prior}{The prios given or calculated for the root node;}
}

\references{
Vladimir N Minin e Marc A Suchard. Fast, accurate and simulation-free
stochastic mapping. Philosophical Transactions of the Royal Society B:
Biological Sciences, 363 (1512): 3985-3995, 2008.
}

\author{Diego Pasqualin \email{dpasqualin@inf.ufpr.br}}

\keyword{phylogenetics}
%\keyword{simulation}
