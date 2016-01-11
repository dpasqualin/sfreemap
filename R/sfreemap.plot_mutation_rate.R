plot_mutation_rate <- function(tree, combination=3) {
    if (!inherits(tree, "multiPhylo")) {
        stop("tree must be a 'multiPhylo' object")
    }
    if (!is.integer(length(tree))/combination) {
        stop("number of trees must be a multiple of 'combinatio'")
    }

    groups <- list()
    for (i in 1:combination) {
        # all trees of the ith group
        group <- tree[seq(i, length(tree), combination)]
        # mutation rate for every branch on every tree of the group
        groups[[i]] <- c(sapply(group, function(t) rowSums(t$mapped.edge.lmt)/t$edge.length))
    }
    groups <- sapply(groups, cbind)
    colnames(groups) <- 1:combination

    to_plot <- data.frame(melt(groups))
    colnames(to_plot) <- c('n', 'color', 'value')
    to_plot$color <- as.character(to_plot$color)

    p <- ggplot(to_plot, aes(x=value, fill=color)) + geom_histogram()
    print(p)

    return(data.frame(to_plot))
}
