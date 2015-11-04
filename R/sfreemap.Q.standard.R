# Return the transition matrix Q calculated empirically
Q_empirical <- function(tree, tip_states, prior, model, tol, omp) {

    states <- tip_states/rowSums(tip_states)

    n_states <- ncol(states)

    # Reorder tree and create a copy called bt. Not sure why
    # phytools need this...
    new_tree <- bt <- reorder.phylo(tree, "cladewise")

    if (!is.binary.tree(bt)) {
        # NOTE: find out what this does...
        bt <- multi2di(bt)
    }

    # NOTE: Add comments to this snipped. It was adapted from phytools
    # and I need to understand what it is doing
    pars <- getPars(bt, states, model, Q=NULL, new_tree, tol, n_states, omp)
    L <- pars$L
    Q <- pars$Q
    logL <- pars$loglik

    # Set prior (pi)
    if (prior[1] == "equal") {
        # set equal
        prior <- setNames(rep(1/n_states, n_states), colnames(L))
    } else if (prior[1] == "estimated") {
        # set from stationary distribution
        prior <- statdist(Q)
    } else {
        # obtain from input
        prior <- prior/sum(prior)
    }

    return (list(Q=Q, prior=prior, logL=logL))

}

# Return the transition matrix Q if it has been passed as a matrix
Q_matrix <- function(tree, tip_states, Q, model, prior, tol) {

    states <- tip_states/rowSums(tip_states)

    n_states <- ncol(states)

    # Reorder tree and create a copy called bt. Not sure why
    # phytools need this...
    new_tree <- bt <- reorder.phylo(tree, "cladewise")

    XX <- getPars(bt, states, model, Q=Q, new_tree, tol, n_states)
    # NOTE: do we need this?
    L <- XX$L
    logL <- XX$loglik

    # Set the priors
    if (prior[1] == "equal") {
        prior <- setNames(rep(1/n_states,n_states), colnames(L)) # set equal
    } else if (prior[1] == "estimated") {
        prior <- statdist(Q) # set from stationary distribution
    } else {
        prior <- prior/sum(prior) # obtain from input
    }

    return (list(Q=Q, prior=prior, logL=logL))
}

# Return the transition matrix Q calculated using a markov chain
Q_mcmc <- function(tree, tip_states, prior, model, gamma_prior, tol, burn_in, sample_freq, vQ, n_simulations, omp) {

    states <- tip_states/rowSums(tip_states)

    n_states <- ncol(states)

    # Reorder tree and create a copy called bt. Not sure why
    # phytools need this copy...
    # TODO: As we set the order of the tree in the main functio, maybe we don't # need to do it here again.
    new_tree <- bt <- reorder.phylo(tree, "cladewise")

    # Define prior
    if (gamma_prior$use.empirical) {
        qq <- apeAce(bt, states, model)$rates
        gamma_prior$alpha <- qq * gamma_prior$beta
    }

    XX <- mcmcQ(bt, states, model, new_tree, tol, n_states, burn_in, sample_freq, n_simulations, vQ, gamma_prior, omp)
    # We compute the likelihood again anyway, maybe we could skip this

    join_Q_with_prior <- function(idx, Q_list, prior_list, logL_list) {
        return (list(Q=Q_list[[idx]]
                     , prior=prior_list[[idx]]
                     , logL=logL_list[[idx]])
                )
    }

    Q <- lapply(XX, function(x) x$Q)
    L <- lapply(XX, function(x) x$L)
    logL <- lapply(XX, function(x) x$loglik)

    if (prior[1] == "equal") {
        prior <- setNames(rep(1/n_states,n_states), colnames(L)) # set equal
        prior <- lapply(1:n_simulations, function(x,y) y, y=prior)
    } else if (prior[1] == "estimated") {
        prior <- lapply(Q, statdist) # set from stationary distribution
    } else {
        prior <- prior/sum(prior) # obtain from input
        prior <- lapply(1:n_simulations, function(x,y) y, y=prior)
    }

    return(lapply(1:n_simulations, join_Q_with_prior
                    , Q_list=Q
                    , prior_list=prior
                    , logL_list=logL))

}

# This function creates a matrix with rownames being states, colnames being
# the tip labels and values being 1, if the tip label has the correspondent
# states and 0 otherwise.
# "possible_states" is useful mainly because DNA can always have A, C, T and G,
# but data might not have one of this letters.
build_states_matrix <- function(tip_labels, tip_states, possible_states=NULL) {
    if (is.null(possible_states)) {
        possible_states <- tip_states
    }
    res <- to.matrix(tip_states, sort(unique(possible_states)))
    return (res[tip_labels,])
}
