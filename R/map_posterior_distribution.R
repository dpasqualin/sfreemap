map_posterior_distribution <- function(base_tree, trees, scale.branches=TRUE, scale.trees=FALSE, parallel=TRUE) {

    # all but the root node
    all_nodes <- unique(base_tree$edge[,2])

    # create tree dimentional matrix for the result
    states <- colnames(base_tree$mapped.edge)
    transitions <- colnames(base_tree$mapped.edge.lmt)
    tree_names <- 1:length(trees)
    node_names <- 1:max(all_nodes)
    result_emr_dim <- c(length(trees), length(states), max(all_nodes))
    result_lmt_dim <- c(length(trees), length(transitions), max(all_nodes))

    result <- list(
        base_tree = base_tree
        , emr = array(NA, result_emr_dim, dimnames=list(tree_names, states))
        , lmt = array(NA, result_lmt_dim, dimnames=list(tree_names, transitions))
        , mr = array(NA, result_lmt_dim, dimnames=list(tree_names, transitions))
    )

    # correspondent nodes of base_tree in tree
    # NA when there is no correspondence
    mymatch <- function(tree) {

        if (is.numeric(scale.trees)) {
            tree <- sfreemap.rescale(tree, scale.trees)
        }

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

    prev_class <- class(trees)
    if (isTRUE(parallel) && !on_windows()) {
        trees <- mclapply(trees, mymatch, mc.cores=detectCores())
    } else {
        trees <- lapply(trees, mymatch)
    }
    class(trees) <- prev_class

    for (node in all_nodes) {
        tree_number <- 1 # index for result
        for (tree in trees) {
            # correspondent node
            cn <- tree$match[as.character(node),2]
            # get branch that ends on cn
            edge <- which(tree$edge[,2]==cn)

            if (!is.na(cn)) {
                # add dwelling times data
                emr <- tree$mapped.edge[edge,]
                result$emr[tree_number,,node] <- emr

                if (!is.null(tree$mapped.edge.lmt)) {
                    # add number of transitions
                    lmt <- tree$mapped.edge.lmt[edge,]
                    result$lmt[tree_number,,node] <- lmt

                    # add mutation rate
                    mr <- lmt/tree$edge.length[edge]
                    result$mr[tree_number,,node] <- mr
                }
            }
            tree_number <- tree_number + 1 # index for result
        }
    }

    if (isTRUE(scale.branches)) {
        result$emr <- freq_to_prob(result$emr)
    }

    return(result)
}
