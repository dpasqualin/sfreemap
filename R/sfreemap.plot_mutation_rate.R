sfreemap.plot_mutation_rate <- function(tree, combination=3) {
    if (!inherits(tree, "multiPhylo")) {
        stop("tree must be a 'multiPhylo' object")
    }
    x <- nrow(tree[[1]]$edge)*combination
    groups <- matrix(0, nrow=nrow(tree[[1]]$edge)*combination, ncol=combination)
    for (i in 1:combination) {
        # all trees of the ith group
        group <- tree[seq(i, length(tree), combination)]
        # mutation rate for every branch on every tree of the group
        groups[,i] <- c(sapply(group, function(t) rowSums(t$mapped.edge.lmt)/t$edge.length))
    }
    colnames(groups) <- 1:combination
    to_plot <- data.frame(melt(groups))
    colnames(to_plot) <- c('n', 'color', 'value')
    to_plot$color <- as.character(to_plot$color)

    p <- ggplot(to_plot, aes(x=value, fill=color)) + geom_bar()
    print(p)

    return(data.frame(groups))
}
