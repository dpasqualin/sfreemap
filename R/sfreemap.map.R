# Input
#   tree    a phylogenetic tree as an object of class "phylo" (from package
#           ape)
sfreemap.map <- function(tree, tip_states, Q='empirical', ...) {

    # Should this program run in parallel?
    parallel <- TRUE
    if (hasArg(parallel)) {
        parallel <- list(...)$parallel
    }

    # how many omp threads should be created?
    omp <- 1
    if (hasArg(omp)) {
        omp <- list(...)$omp
    }

    # tree sanity check
    if ('multiPhylo' %in% class(tree)) {
        # For Just call the same program multiple times...
        if (parallel == TRUE) {
            cores <- detectCores()
            mtrees <- mclapply(tree, sfreemap.map, tip_states, Q, ..., mc.cores=cores)
        } else {
            mtrees <- lapply(tree, sfreemap.map, tip_states, Q, ...)
        }

        # When Q=mcmc we will have length(trees)*n_simulation trees at the end
        # of the execution. Instead of lists of multPhylo objects we want to
        # return one multiPhylo object with all trees.
        if (Q == 'mcmc') {
            mtrees <- c(mapply(c, mtrees))
        }

        class(mtrees) <- "multiPhylo"
        return(mtrees)

    } else if (! "phylo" %in% class(tree)) {
        stop("'tree' should be an object of class 'phylo'")
    } else if (!is.rooted(tree)) {
        stop("'tree' must be rooted")
    }

    # Defining the prior distribution for the root node of the tree,
    # also known as "pi"
    prior <- "equal"
    if (hasArg(prior)) {
        prior <- list(...)$prior
    }

    # tol gives the tolerance for zero elements in Q.
    # (Elements less then tol will be reset to tol)
    tol <- 1e-8
    if (hasArg(tol)) {
        tol <- list(...)$tol
    }

    # prior a list containing alpha and beta parameters for the gamma
    # prior distribution on the transition rates in Q. Note that alpha
    # and beta can be single values or vectors, if different prior
    # are desired for each value in Q
    gamma_prior <- list(alpha=1, beta=1, use.empirical=FALSE)
    if (hasArg(gamma_prior)) {
        pr <- list(...)$gamma_prior
        gamma_prior[names(pr)] <- pr
    }

    # burn_in for the MCMC
    burn_in <- 1000
    if (hasArg(burn_in)) {
        burn_in <- list(...)$burn_in
    }

    # sample_freq for the MCMC
    sample_freq <- 100
    if (hasArg(sample_freq)) {
        sample_freq <- list(...)$sample_freq
    }

    # number os simulations for the MCMC
    n_simulations <- 100
    if (hasArg(n_simulations)) {
        n_simulations <- list(...)$n_simulations
    }

    # For now only "symmetrical" model is accepted
    model <- "SYM"
    #if (hasArg(model)) {
    #    model <- list(...)$model
    #}

    # A single numeric value or a vector containing the (normal)
    # sampling variances for the MCMC
    vQ <- 0.1
    if (hasArg(vQ)) {
        vQ <- list(...)$vQ
    }

    # Define the tip states as a matrix
    if (!is.matrix(tip_states)) {
        tip_states <- build_states_matrix(tree, tip_states)
    } else {
        tip_states <- tip_states[tree$tip.label,]
    }

    # Defining Q
    if (is.character(Q) && (Q == "empirical")) {
        # Phytools would replicate this result nsim times, but for now
        # we will return just one result.
        QP <- Q_empirical(tree, tip_states, prior, model, tol, omp)
    } else if (is.character(Q) && (Q == "mcmc")) {
        # TODO: This function will generate many Qs. We have to decide how to
        # deal with it. Maybe just run the program for every Q?
        QP <- Q_mcmc(tree, tip_states, model, prior, gamma_prior, tol, burn_in
                     , sample_freq, vQ, n_simulations)
        # Call sfreemap.map for each {Q,prior} returned by the mcmc simulation
        params <- list(...)
        params$tree <- tree
        params$tip_states <- tip_states
        res <- function(QP) {
            params$Q <- QP$Q
            params$prior <- QP$prior
            return (do.call(sfreemap.map, params))
        }

        if (parallel == TRUE) {
            mtrees <- mclapply(QP, res, mc.cores=detectCores())
        } else {
            mtrees <- lapply(QP, res)
        }
        class(mtrees) <- "multiPhylo"
        return(mtrees)

    } else if (is.matrix(Q)) {
        # Phytools would replicate this result nsim times, but for now
        # we will return just one result.
        QP <- Q_matrix(tree, tip_states, Q, model, prior, tol)
    } else {
        stop("Unrecognized format for 'Q'")
    }

    # Set the final value
    Q <- QP$Q
    prior <- QP$prior
    logL <- QP$logL

    # Vector with rewards
    rewards <- rep(1,nrow(Q))
    # FIXME: according to Minin & Suchar article this was suppose to be a
    # parameters and user could provide different rewards for the states. But
    # we tried with different values and the result just doesn't make sense.
    # We've tried to reach the authors but got no answer on this matter.
    #if (hasArg(rewards)) {
    #    rewards <- list(...)$rewards
    #    if (length(rewards) != nrow(Q)) {
    #        stop("The rewards vector should represent the states of Q")
    #    }
    #}
    names(rewards) <- colnames(Q)

    # Acquire more info about the tree.
    tree_extra <- list(
        states = tip_states
        , n_states = nrow(Q)
        , n_edges = length(tree$edge.length)
        , n_tips = nrow(tip_states)
        , n_nodes = nrow(tip_states) + tree$Nnode
        , rewards = rewards
    )

    # Reorder the tree so the root is the first row of the matrix.
    # We save the original order to make sure we have the result
    # in same order of the tree;
    tree <- reorder(tree, 'pruningwise')

    # Step 1
    # Compute Eigen values and Eigen vectors for the transition rate matrix
    # Q, assuming Q is symmetrical
    Q_eigen <- eigen(Q, TRUE, only.values = FALSE)
    # The inverse of Q_eigen vectors
    Q_eigen[['vectors_inv']] <- solve(Q_eigen$vectors)

    # Step 2
    # Compute P(tp), the transistion probability, for each edge length t of T
    # Q = U X diag(d1,...,dm) X U**-1
    # U are the eigenvectors of Q
    # d1,...,dm are the eigenvalues of Q
    # diag(d1,...,dm) is a diagonal matrix with d1,...,dm on it's main
    # diagonal
    # For now I'm doing this just where it is needed, in the
    # fractional_likelihood calculation. I don't know yet if I'm going
    # to need it somewhere else. It might be useful to compute it here
    # and use it in all different places if that's the case.
    MAP <- list()

    MAP[['Q']] <- Q
    MAP[['prior']] <- prior

    # Transistion probabilities
    MAP[['tp']] <- transition_probabilities(Q_eigen, tree$edge.length, omp)

    # Step 3
    # Employing the eigen decomposition above compute E(h, tp*) for
    # each edge b* in the set of interest Omega using equation 2.4
    # (expected number of markov transitions) and equation 2.12
    # (expected markov rewards).
    MAP[['h']] <- func_H(Q, Q_eigen, tree, tree_extra, omp)

    # Step 4 and 5
    # Traverse the tree once and calculate Fu and Sb for each node u and
    # each edge b;
    # Compute the data likelihood Pr(D) as the dot product of Froot and root
    # distribution pi.
    MAP[['fl']] <- fractional_likelihoods(tree, tree_extra, Q, Q_eigen
                                          , prior, MAP$tp, tol)

    # Posterior restricted moment for branches
    # This is the "per branch" expected value for lmt and emr
    MAP[['prm']] <- posterior_restricted_moment(tree, tree_extra, MAP, omp)

    # This is the global mean, not sure why we need it..
    MAP[['ev']] <- expected_value(tree, Q, MAP)

    # Let's set the elements back to the original tree
    tree[['Q']] <- Q
    tree[['prior']] <- prior
    tree[['logL']] <- logL

    tree[['mapped.edge']] <- MAP[['ev']]$emr
    tree[['mapped.edge.lmt']] <- MAP[['ev']]$lmt

    # Return the tree in the original order
    return (sfreemap.reorder(tree, 'cladewise'))
}

# The final answer!
expected_value <- function(tree, Q, map) {

    likelihood <- map[['fl']][['L']]
    # posterior restricted moment...
    prm <- map[['prm']]

    EV = list()
    #EV[['lmt']] <- apply(prm[['lmt']], 2, sum) / likelihood
    #EV[['emr']] <- apply(prm[['emr']], 2, sum) / likelihood
    EV[['lmt']] <- prm[['lmt']] / likelihood
    EV[['emr']] <- prm[['emr']] / likelihood

    # the rownames of the mapped objects
    mapped_names <- paste(tree$edge[,1], ",", tree$edge[,2], sep="")
    rownames(EV[['lmt']]) <- rownames(EV[['emr']]) <- mapped_names
    colnames(EV[['emr']]) <- colnames(EV[['lmt']]) <- colnames(Q)

    return (EV)
}
