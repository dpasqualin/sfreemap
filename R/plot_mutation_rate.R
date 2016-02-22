plot_mutation_rate <- function(tree, combination=3) {
    if (!inherits(tree, "multiPhylo")) {
        stop("tree must be a 'multiPhylo' object")
    }
    if (!is.integer(length(tree))/combination) {
        stop("number of trees must be a multiple of 'combination'")
    }


    groups <- list()
    for (i in 1:combination) {
        # all trees of the ith group
        group <- tree[seq(i, length(tree), combination)]
        # mutation rate for every branch on every tree of the group
        group <- c(sapply(group, function(t) rowSums(t$mapped.edge.lmt)/t$edge.length))
        # remove outliers
        group <- group[group <= boxplot.stats(group)$stats[5]]
        # make sure all groups have the same number of values
        group <- c(group, rep(NA, number_of_edges-length(group)))
        # save
        groups[[i]] <- group
    }
    groups <- sapply(groups, cbind)
    colnames(groups) <- 1:combination

    to_plot <- data.frame(melt(groups))
    colnames(to_plot) <- c('n', 'color', 'value')
    to_plot$color <- as.character(to_plot$color)

    p <- ggplot(to_plot, aes(x=value, fill=color)) +
            geom_density(position='stack', na.rm=TRUE)

    print(p)

    return(data.frame(to_plot))
}
