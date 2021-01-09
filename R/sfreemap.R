# Input
#   tree    a phylogenetic tree as an object of class "phylo" (from package
#           ape)
sfreemap <- function(tree, tip_states, Q=NULL, type="standard", model="SYM", method="empirical", ...) {

    # Am I running on windows? Windows does not have support for the kind of
    # parallelism we are using

    # Should this program run in parallel?
    if (hasArg(parallel)) {
        parallel <- list(...)$parallel
        if (isTRUE(parallel) && !support_parallel_mode()) {
            ## parallel requested but parallel mode not supported
            warning('parallel mode is not available on this machine.', call. = FALSE)
            parallel <- FALSE
        }
    } else {
        parallel <- support_parallel_mode()
    }

    # When running in parallel, choose how many cores do use.
    # default to all cores available on the machine
    mc.cores <- detectCores()
    if (hasArg(mc.cores)) {
        tmp <- list(...)$mc.cores
        # ignore things like NULL
        if (is.numeric(tmp)) {
            mc.cores <- tmp
        }
    }
    if (mc.cores == 1) {
        parallel <- FALSE
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

    # a helper function to call sfreemap multiple times, combining
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

        return (do.call(sfreemap, params))
    }

    # helper to decide whether to call 'call_multiple' in serial or parallel
    serial_or_parallel <- function(times, tree, tip_states, Q, prior) {
        if (parallel) {
            mtrees <- mclapply(1:times, call_multiple, tree, tip_states, Q
                                      , prior, mc.cores=mc.cores)
        } else {
            mtrees <- lapply(1:times, call_multiple, tree, tip_states, Q, prior)
        }
        return (fix_return(mtrees))
    }

    # with some combination of parameters we might have a list of multiPhylo
    # objects (a list of a list), so we need to convert it to a single
    # multiPhylo object
    fix_return <- function(mtrees) {
        tmp <- mtrees[[1]]
        if (inherits(tmp, "multiPhylo") || inherits(tmp, "list")) {
            mtrees <- c(mapply(c, mtrees))
        }

        if (length(mtrees) > 1) {
            class(mtrees) <- c("sfreemap", "multiPhylo")
        } else {
            mtrees <- mtrees[[1]]
        }

        return (mtrees)
    }

    # Everything below these tests assume the program is running on with a
    # single tree, single rate matrix and single tip label. So here we check
    # parameters and call sfreemap multiple times if needed.
    if (inherits(tree, "multiPhylo")) {
        # if 'multiPhylo', call sfreemap for each tree
        return(serial_or_parallel(length(tree), tree, tip_states, Q, prior))
    } else if (inherits(Q, "list")) {
        # if multiple rate matrix, call sfreemap for each one.
        # serial_or_parallell will handle the case when we have an equal number
        # of rate matrices and trees, where sfreemap should match tree 1
        # with rate matrix 1, 2 with 2, and so on..
        return(serial_or_parallel(length(Q), tree, tip_states, Q, prior))
    } else if (inherits(prior, "list")) {
        # we can run multiple priors on trees and Q rate matrices
        return(serial_or_parallel(length(prior), tree, tip_states, Q, prior))
    } else if (all(inherits(tip_states, "matrix"), ncol(tip_states) > 1)) {
        # if single tree but multiple tip_states (dna type), call sfreemap
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

    # FIXME: this was not tested yet, not sure it it works
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

    # FIXME: this was not tested yet, not sure it it works
    # A single numeric value or a vector containing the (normal)
    # sampling variances for the MCMC
    vQ <- 0.1
    if (hasArg(vQ)) {
        vQ <- list(...)$vQ
    }

    # We set the class here so we can use functions like reorder.sfreemap
    class(tree) <- c("sfreemap", "phylo")

    if (all(!is.null(Q), is.matrix(Q))) {
        QP <- Q_matrix(tree, tip_states, Q, model, prior, tol, type)
    } else if (type == 'standard') {
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
    }

    # NOTE: It's important to notice that below this point it is garantee that
    # we are dealing with a single tree (a "phylo" object), a single Q matrix,
    # a single prior and a single character.

    if (type == "dna") {
        states <- c("a", "c", "t", "g", "-")
    } else {
        states <- NULL
    }
    tip_states <- build_states_matrix(tree$tip.label, tip_states, states)

    # Set the final value
    Q <- QP$Q
    prior <- QP$prior
    logL <- QP$logL

    # Vector with rewards
    rewards <- rep(1,nrow(Q))
    if (hasArg(rewards)) {
        rewards <- list(...)$rewards
        if (length(rewards) != nrow(Q)) {
            stop("The rewards vector should represent the states")
        }
    }
    names(rewards) <- colnames(Q)

    # Defining the transitions of interest. By default, all transitions.
    QL <- matrix(1, nrow=nrow(Q), ncol=ncol(Q))
    diag(QL) <- 0
    if (hasArg(QL)) {
        QL <- list(...)$QL
        if (!all(dim(Q) == dim(QL))) {
            stop("QL must have same dimensions as Q")
        }
    }
    rownames(QL) <- colnames(QL) <- rownames(Q)

    # Acquire more info about the tree.
    tree_extra <- list(
        states = tip_states
        , n_states = nrow(Q)
        , n_edges = length(tree$edge.length)
        , n_tips = nrow(tip_states)
        , n_nodes = nrow(tip_states) + tree$Nnode
    )

    # Reorder the tree so the root is the first row of the matrix.
    # We save the original order to make sure we have the result
    # in same order of the tree;
    tree <- reorder(tree, order='pruningwise')
    # Let's set the elements back to the original tree
    tree[['Q']] <- Q
    tree[['prior']] <- prior
    tree[['logL']] <- logL

    # Step 1
    # Compute Eigen values and Eigen vectors for the transition rate matrix
    # Q, assuming Q is symmetrical
    Q_eigen <- eigen(Q, TRUE, only.values = FALSE, symmetric = TRUE)
    # The inverse of Q_eigen vectors
    Q_eigen[['vectors_inv']] <- solve(Q_eigen$vectors)

    MAP <- list()
    # Step 2
    # Compute P(tp), the transistion probability, for each edge length t of T
    MAP[['tp']] <- transition_probabilities(Q_eigen, tree$edge.length, omp)

    # Step 4 and 5
    MAP[['fl']] <- fractional_likelihoods(tree, tree_extra, Q, Q_eigen
                                          , prior, MAP$tp, tol)

    # FIXME: for some unknown reason if I remove this line I get a 'out of
    # bounds' error in posterior_restricted_moment(). This line doesn't need to
    # exist
    MAP[['h']] <- list()

    # Build dwelling times
    tree[['mapped.edge']] <- tree[['mapped.edge.lmt']] <- tname <- NULL
    states <- rownames(QL)
    n_states <- tree_extra$n_states
    multiplier <- matrix(0, nrow=n_states, ncol=n_states)

    for (i in 1:n_states) {
        for (j in 1:n_states) {

            if (i == j) {
                value <- rewards[i] # dwelling times
            } else {
                value <- QL[i,j] * Q[i,j] # number of transitions
            }

            if (value == 0) {
                next
            }
            multiplier[i,j] <- value

            # Step 3
            h <- func_H(multiplier, Q_eigen, tree, tree_extra, omp)
            prm <- posterior_restricted_moment(h, tree, tree_extra, MAP, omp)
            ev <- expected_value(tree, Q, MAP, prm)

            if (i == j) {
                # dwelling times
                tree[['mapped.edge']] <- cbind(tree[['mapped.edge']], ev)
            } else {
                # number of transitions
                state_from <- states[i]
                state_to <- states[j]
                tname <- c(tname, paste(state_from, state_to, sep=','))
                tree[['mapped.edge.lmt']] <- cbind(tree[['mapped.edge.lmt']], ev)
            }

            multiplier[i,j] <- 0
        }
    }
    bname <- paste(tree$edge[,1], ",", tree$edge[,2], sep="")
    colnames(tree[['mapped.edge']]) <- names(rewards)
    colnames(tree[['mapped.edge.lmt']]) <- tname
    rownames(tree[['mapped.edge']]) <- rownames(tree[['mapped.edge.lmt']]) <- bname

    # Return the tree in the original order
    return (reorder(tree, order='cladewise'))
}

# The final answer!
expected_value <- function(tree, Q, map, prm) {

    likelihood <- map[['fl']][['L']]

    ret <- prm / likelihood

    # the rownames of the mapped objects
    names(ret) <- paste(tree$edge[,1], ",", tree$edge[,2], sep="")

    return(ret)
}
