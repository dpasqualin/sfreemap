\name{sfreemap.corals.trees}
\docType{data}
\alias{sfreemap.corals.trees}
\title{Partial sequences of the 12S and 28S rRNA for 108 species of scleractinian corals and one sea anemone}
\description{
This dataset comprises partial sequences of the 12S and 28S rRNA for 108 species of scleractinian corals and one sea anemone (used as outgroup) and one binary character (0 = solitary and 1 = colonial species). The posterior distribution of trees was sampled via MCMC using the CIPRES Science Gateway ( http://www.phylo.org/sub_sections/portal/) implementation of BEAST (Drumond et al., 2007). Sequences were aligned with reference to a consensus secondary structure model (details in Barbeitos et al., 2010) and data were partitoned into loop and stems. Analyses were run under the Hasegwa Kishino Yano substitution model (HKY, Hasegawa et al., 1985) with rate heterogeneity among sites modeled by a gamma distribution discretized into 4 categories and a relaxed (uncorrelated) lognormal molecular clock. The chain was run for 10 million generations and sampled every 1,000 generations. Stationarity was assessed using Tracer (Rambaut et al., 2014) after the first 10\% chain was discarded. In order to maximize independence among samples, the ensuing distribution of trees was "thinned" to 901 trees, i.e, sub-sampled until the integrated autocorrelated time (IACT) dropped to 0 (see Pagel et al., 2006 for details).
}
\usage{sfreemap.corals.trees}
\format{A "multiPhylo" object containing 901 phylogenetic trees}
\references{Barbeitos, Marcos S., Sandra L. Romano, and Howard R. Lasker. "Repeated Loss of Coloniality and Symbiosis in Scleractinian Corals." Proceedings of the National Academy of Sciences 107, no. 26 (2010): 11877-82.

Drummond, Alexei, and Andrew Rambaut. "BEAST: Bayesian Evolutionary Analysis by Sampling Trees." BMC Evolutionary Biology 7 (2007): 214.

Hasegawa, M., H. Kishino, and T. Yano. "Dating the Human-ape Split by a Molecular Clock of Mitochondrial DNA." Journal of Molecular Evolution 22 (1985): 160-74.
}
\keyword{datasets}
