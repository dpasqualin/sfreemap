sfreemap.map_posterior_distribution <- function(base_tree, trees, scale=TRUE) {
    n_trees <- length(trees)
    states <- colnames(base_tree$mapped.edge)
    n_states <- length(states)
    # all but the root node
    all_nodes <- unique(base_tree$edge[,2])[-1]
    n_nodes <- max(all_nodes)

    result_dim <- c(n_trees, n_states, n_nodes)
    result <- list(
        emr = array(NA, result_dim, dimnames=list(NULL, states))
        , lmt = array(NA, result_dim, dimnames=list(NULL, states))
    )

    # correspondent nodes of base_tree in tree
    # NA when there is no correspondent
    mymatch <- function(tree) {
        # match of internal nodes
        internal <- matchNodes(base_tree, tree, method='descendants')
        # match of tips
        ta <- base_tree$tip.label
        tb <- tree$tip.label
        tips <- cbind(seq_along(ta), match(ta, tb))
        # concatenate internal and tip nodes
        match <- rbind(internal, tips)

        # this names are useful to index the match vector by starting node
        rownames(match) <- match[,1]
        tree[['match']] <- match

        return(tree)
    }

    trees <- mclapply(trees, mymatch)
    class(trees) <- 'multiPhylo'

    for (node in all_nodes) {
        tree_number <- 1 # index for result
        for (tree in trees) {
            # correspondent node
            cn <- tree$match[as.character(node),2]

            if (!is.na(cn)) {
                # search for branch ending in "cn" on tree and add corresponding
                # values for states
                branch_names <- rownames(tree$mapped.edge)
                pattern <- paste(',', cn, '$', sep='')
                branch <- grepl(pattern, branch_names)
                result$emr[tree_number,,node] <- tree$mapped.edge[branch,]
                result$lmt[tree_number,,node] <- tree$mapped.edge.lmt[branch,]
            }
            tree_number <- tree_number + 1 # index for result
        }
    }

    if (isTRUE(scale)) {
        result$emr <- freq_to_prob(result$emr)
    }

    return(result)
}
