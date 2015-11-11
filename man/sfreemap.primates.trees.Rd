\name{sfreemap.primates.trees.dna}
\docType{data}
\alias{sfreemap.primates.trees.dna}
\title{Taxonomic sampling of Old World monkeys}
\description{
Taxonomic sampling of Old World monkeys (parvorder Catarrhini, which includes humans, gorillas, chimpanzees, orangutans, gibbons, baboons, macaques and langurs) used by Pagel et al. (2006), plus one New World monkey (the tufted capuchin, Cebus apella), which was used as outgroup.  Complete and partial cytochrome b sequences were downloaded from GeneBank and aligned using the online version of MAFFT (Katoh et al., 2007). The alignment was converted to amino acid sequences and trimmed so that the reading frame starts at the first site. Data were partitioned by codon position and analyzed under the general-time reversible (GTR, Tavaré, 1986) substitution model with rate heterogeneity among sites also modeled by a 4-category gamma distribution and rate variation across branches accommodated by a relaxed (uncorrelated) lognormal molecular clock.  Four independent chains were run for 100 million generation and sampled every 100,000 generations. Stationarity and convergence among chains were assessed using Tracer (after a 10\% burnin) and the posterior distributions of trees were also "thinned". These posteriors were subsequently combined using LogCombiner (Drummond 2006) into a single distribution of 3,600 trees.
}
\usage{sfreemap.primates.trees.dna}
\format{A "multiPhylo" object containing 3596 phylogenetic trees}
\references{Barbeitos, Marcos S., Sandra L. Romano, and Howard R. Lasker. "Repeated Loss of Coloniality and Symbiosis in Scleractinian primates." Proceedings of the National Academy of Sciences 107, no. 26 (2010): 11877–82.

Pagel, Mark, and Andrew Meade. "Bayesian Analysis of Correlated Evolution of Discrete Characters by Reversible-Jump Markov Chain Monte Carlo." Am Nat 167, no. 6 (2006): 808–25.

Katoh, Kazutaka, Kei-ichi Kuma, Hiroyuki Toh, and Takashi Miyata. “MAFFT Version 5: Improvement in Accuracy of Multiple Sequence Alignment.” Nucleic Acids Research 33, no. 2 (2005): 511–18. doi:10.1093/nar/gki198.

Tavaré, S. "Some Probabilistic and Statistical Problems in the Analysis of DNA Sequences." American Mathematical Society: Lectures on Mathematics in the Life Sciences 17 (1986): 57–86.

Drummond, Alexei, and Andrew Rambaut. "BEAST: Bayesian Evolutionary Analysis by Sampling Trees." BMC Evolutionary Biology 7 (2007): 214.

}
\keyword{datasets}
