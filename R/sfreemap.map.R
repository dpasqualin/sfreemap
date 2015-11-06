# Input
#   tree    a phylogenetic tree as an object of class "phylo" (from package
#           ape)
sfreemap.map <- function(tree, tip_states, Q=NULL, type="standard", model="SYM", method="empirical", ...) {

    # Am I running on windows? Windows does not have support for the kind of
    # parallelism we are using
    i_am_windows <- Sys.info()['sysname'] == 'Windows'

    # Should this program run in parallel?
    parallel <- !i_am_windows
    if (hasArg(parallel)) {
        parallel <- list(...)$parallel
        if (all(parallel, i_am_windows)) {
            warning('parallel mode is not available on windows.', call. = FALSE)
            parallel <- FALSE
        }
    }

    # how many omp threads should be created?
    omp <- 1
    if (hasArg(omp)) {
        omp <- list(...)$omp
    }

    # Defining the prior distribution for the root node of the tree,
    # also known as "pi"
    prior <- "equal"
    if (hasArg(prior)) {
        prior <- list(...)$prior
    }

    # available types and models
    dna_models <- names(getModels())
    standard_models <- c('SYM', 'ER', 'ARD')
    valid_models <- list(
        "standard" = standard_models
        , "dna" = dna_models
    )

    # check for tree object class
    if (all(!inherits(tree, "phylo"), !inherits(tree, "multiPhylo"))) {
        stop("'tree' should be an object of class 'phylo' or 'multiPhylo'")
    }

    # check for type
    if (! type %in% names(valid_models)) {
        stop('Unknown type', type)
    }

    # check for model
    if (! model %in% valid_models[[type]]) {
        stop('Unknown model ', model)
    }

    if (all(inherits(prior, "list"), !inherits(Q, "list"))) {
        stop("if 'prior' is a list 'Q' should be a list with of same size.")
    }

    if (all(inherits(prior, "list"), inherits(Q, "list"), length(prior) != length(Q))) {
        stop("if 'prior' and 'Q' are lists, their number of elements must match.")
    }

    if (all(inherits(Q, "list"), inherits(tree, "multiPhylo"), length(Q) != length(tree))) {
        stop("if 'Q' is a list and 'tree' is a 'multiPhylo' object, their number of elements must match.")
    }

    # a helper function to call sfreemap.map multiple times, combining
    # trees with Q rate matrices and priors
    call_multiple <- function(idx, tree, tip_states, Q, prior) {
        if (inherits(tree, "multiPhylo")) {
            tree <- tree[[idx]]
        }

        if (all(inherits(tip_states, "matrix"), ncol(tip_states) > 1)) {
            tip_states <- tip_states[,idx]
        }

        if (all(!is.null(Q), inherits(Q, "list"))) {
            Q <- Q[[idx]]
        }

        if (inherits(prior, "list")) {
            prior <- prior[[idx]]
        }

        params <- list(
            "tree" <- tree
            , "tip_states" <- tip_states
            , "Q" = Q
            , "prior" = prior
            , "model" = model
            , "type" = type
            , "method" = method
            , "..." = ...
        )

        return (do.call(sfreemap.map, params))
    }

    # helper to decide whether to call 'call_multiple' in serial or parallel
    serial_or_parallel <- function(times, tree, tip_states, Q, prior) {
        if (parallel) {
            mtrees <- mclapply(1:times, call_multiple, tree, tip_states, Q
                                      , prior, mc.cores=detectCores())
        } else {
            mtrees <- lapply(1:times, call_multiple, tree, tip_states, Q, prior)
        }
        return (fix_return(mtrees))
    }

    # with some combination of parameters we might have a list of multiPhylo
    # objects (a list of a list), so we need to convert it to a single
    # multiPhylo object
    fix_return <- function(mtrees) {
        if (inherits(mtrees[[1]], "multiPhylo")) {
            mtrees <- c(mapply(c, mtrees))
        }
        class(mtrees) <- "multiPhylo"
        return (mtrees)
    }

    # Everything below these tests assume the program is running on with a
    # single tree, single rate matrix and single tip label. So here we check
    # parameters and call sfreemap.map multiple times if needed.
    if (inherits(tree, "multiPhylo")) {
        # if 'multiPhylo', call sfreemap.map for each tree
        return(serial_or_parallel(length(tree), tree, tip_states, Q, prior))
    } else if (inherits(Q, "list")) {
        # if multiple rate matrix, call sfreemap.map for each one.
        # serial_or_parallell will handle the case when we have an equal number
        # of rate matrices and trees, where sfreemap.map should match tree 1
        # with rate matrix 1, 2 with 2, and so on..
        return(serial_or_parallel(length(Q), tree, tip_states, Q, prior))
    } else if (inherits(prior, "list")) {
        # we can run multiple priors on trees and Q rate matrices
        return(serial_or_parallel(length(prior), tree, tip_states, Q, prior))
    } else if (all(inherits(tip_states, "matrix"), ncol(tip_states) > 1)) {
        # if single tree but multiple tip_states (dna type), call sfreemap.map
        # for each set of tip label
        return(serial_or_parallel(ncol(tip_states), tree, tip_states, Q, prior))
    } else if (!is.rooted(tree)) {
        # all trees must be rooted
        stop("'tree' must be rooted")
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
    n_simulation <- 10
    if (hasArg(n_simulation)) {
        n_simulation <- list(...)$n_simulation
    }

    # A single numeric value or a vector containing the (normal)
    # sampling variances for the MCMC
    vQ <- 0.1
    if (hasArg(vQ)) {
        vQ <- list(...)$vQ
    }

    # FIXME: This function should not exist, it's too complicated.
    # The problem is that for function Q_dna the tip_states should not be a
    # matrix, but a character vector instead. For other situtations it should be
    # a matrix. One option is to adapt Q_dna to work with a matrix..
    check_tip_states <- function() {
        if (type == "dna") {
            possible_states <- c("a", "c", "t", "g", "-")
        } else {
            possible_states <- NULL
        }
        tip_states <- build_states_matrix(tree$tip.label, tip_states, possible_states)
        return(tip_states)
    }

    if (all(!is.null(Q), is.matrix(Q))) {
        tip_states <- check_tip_states()
        QP <- Q_matrix(tree, tip_states, Q, model, prior, tol)
    } else if (type == 'standard') {
        tip_states <- check_tip_states()
        # standard data type has currently two ways of estimating the rate
        # matrix
        if (method == "empirical") {
            QP <- Q_empirical(tree, tip_states, prior, model, tol, omp)
        } else if (method == "mcmc") {
            QP <- Q_mcmc(tree, tip_states, prior, model, gamma_prior, tol, burn_in
                         , sample_freq, vQ, n_simulation, omp)
            Q <- lapply(QP, function(x) x$Q)
            prior <- lapply(QP, function(x) x$prior)
            # TODO: how to pass on logL?
            return(serial_or_parallel(length(QP), tree, tip_states, Q, prior))
        }
    # Estimating Q when using nucleotide data
    } else if (all(type == "dna", is.null(Q))) {
        QP <- Q_dna(tip_states, tree, model, tol)
        tip_states <- check_tip_states()
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
