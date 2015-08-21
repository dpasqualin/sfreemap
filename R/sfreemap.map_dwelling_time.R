sfreemap.map_dwelling_time <- function(base_tree, trees, scale=TRUE) {
    n_trees <- length(trees)
    states <- colnames(base_tree$mapped.edge)
    n_states <- length(states)
    # all but the root node
    all_nodes <- unique(base_tree$edge[,2])[-1]

    res <- list()

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

    trees <- lapply(trees, mymatch)
    class(trees) <- 'multiPhylo'

    for (node in all_nodes) {
        result <- matrix(0, nrow=n_trees, ncol=n_states)
        colnames(result) <- states
        count <- 1 # index for result
        for (tree in trees) {
            # correspondent node
            cn <- tree$match[as.character(node),2]

            if (is.na(cn)) {
                # there is no corresponding node on tree
                result[count, ] <- NA
            } else {
                # search for branch ending in "cn" on tree and add corresponding
                # values for states
                branch_names <- rownames(tree$mapped.edge)
                pattern <- paste(',', cn, '$', sep='')
                branch <- grepl(pattern, branch_names)
                if (isTRUE(scale)) {
                    result[count,] <- rule_of_three(tree$mapped.edge[branch,])
                } else {
                    result[count,] <- tree$mapped.edge[branch,]
                }
            }
            count <- count + 1 # index for result
        }
        node_name <- paste('n', as.character(node), sep='')
        res[[node_name]] <- result
    }
    return(res)
}
