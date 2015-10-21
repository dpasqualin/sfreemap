sfreemap.map_posterior_distribution <- function(base_tree, trees, scale=TRUE) {

    # all but the root node
    all_nodes <- unique(base_tree$edge[,2])

    # create tree dimentional matrix for the result
    states <- colnames(base_tree$mapped.edge)
    tree_names <- 1:length(trees)
    result_dim <- c(length(trees), length(states), max(all_nodes))
    result <- list(
        base_tree = base_tree
        , emr = array(NA, result_dim, dimnames=list(tree_names, states))
        , lmt = array(NA, result_dim, dimnames=list(tree_names, states))
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
        colnames(match) <- c('base_tree', 'me')

        tree[['match']] <- match

        return(tree)
    }

    trees <- mclapply(trees, mymatch, mc.cores=detectCores())
    class(trees) <- 'multiPhylo'

    for (node in all_nodes) {
        tree_number <- 1 # index for result
        for (tree in trees) {
            # correspondent node
            cn <- tree$match[as.character(node),2]

            if (!is.na(cn)) {
                # search for branch ending in "cn" on tree and add corresponding
                # values for states
                emr <- tree$mapped.edge[tree$edge[,2]==cn,]
                result$emr[tree_number,,node] <- emr

                if (!is.null(tree$mapped.edge.lmt)) {
                    lmt <- tree$mapped.edge.lmt[tree$edge[,2]==cn,]
                    result$lmt[tree_number,,node] <- lmt
                }
            }
            tree_number <- tree_number + 1 # index for result
        }
    }

    if (isTRUE(scale)) {
        result$emr <- freq_to_prob(result$emr)
    }

    return(result)
}
