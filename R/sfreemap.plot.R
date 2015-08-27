# TODO: look for p-value on the internet, probably this entire function
# is already implemented somewhere
sfreemap.plot_distribution <- function(node, states=NULL, conf_level=90
	, number_of_ticks=20, ...) {

	# TODO: add sanity check for parameters

	# all states or states passed as argument
	if (is.null(states)) {
		states <- colnames(node)
	} else {
	 	states <- c(states)
	}

	# first divide the dataset
	ticks <- seq(0, 100, 100/number_of_ticks)
	to_plot <- data.frame(x=ticks, alpha=rep(FALSE, length(ticks)))
	for (cont in 1:length(states)) {
		state <- states[cont]
		data <- get_state_data(node, state, conf_level, ticks)

		to_plot[[as.character(state)]] <- data$final_prob
		to_plot$alpha[data$final_idx] <- TRUE
	}

	melted <- melt(to_plot, id=c('x', 'alpha'))

	# graph config
	title <- "Posterior Distribution of Branch Lengths"
	subtitle <- paste("(NA: ", data$final_na_percent, "%)", sep="")
	xlabel <- "Dwelling time (% of branch length)"
	ylabel <- "Probability"

	p <- ggplot(melted, aes(x=x, y=value, fill=variable)) +
			# define the alpha for bars inside and outside HPD
	 		scale_alpha_discrete(range=c(0.3, 0.6), guide=FALSE) +
			# x+2.5 ensures that the beginning of the bar will be at the
			# beginning of the interval, and not at the middle
			geom_bar(stat="identity", position="identity", aes(alpha=alpha, x=x+2.5)) +
			# x and y labels
			xlab(xlabel) + ylab(ylabel) +
			# add title and subtitle
			ggtitle(bquote(atop(.(title), atop(italic(.(subtitle)), "")))) +
			#ggtitle("Distribution of branch length across trees") +
			scale_fill_discrete(name = "States") +
			# set the breaks (ticks) at x axis
			scale_x_discrete(breaks=ticks) +
			# put legend at the bottom
			theme(legend.position="bottom")

	return(p)
}

get_state_data <- function(node, state, conf_level, ticks) {
	data <- node[,as.character(state)]

	final_conf_level <- 0.0
	final_idx <- c()
	final_na_percent <- 0.0

	# group values into the closest tick
	group <- gen_group(data, ticks)
	# just make sure all ticks have a value, even if it's zero

	# translate frequencies into probabilities of seem a particular values
	# in a tick
	prob <- freq_to_prob(group)

	# find the max conf_level we can get with probabilities calculated
	# the maximum conf_level we could find
	partial <- get_max_percent(prob, conf_level)
	# the indexes in prob objected that generated the conf_level above
	final_conf_level <- partial[['conf_level']]
	final_idx <- partial[['prob_idx']]

	final_na_percent <- round(tail(prob, 1), 2)
	# remove NA values from final_prob
	final_prob <- prob[1:(length(prob)-1)]

	return(list(
		final_conf_level=final_conf_level
		, final_idx=final_idx
		, final_prob=final_prob
		, final_na_percent=final_na_percent
	))
}

gen_group <- function(data, ticks) {
	group <- table(findInterval(data, ticks), useNA='always')
	# make sure we have all possible groups
	ng <- names(group)
	res <- sapply(1:length(ticks), function(x) {
		return (ifelse(x %in% ng, group[[as.character(x)]], 0))
	})
	res <- c(res, tail(group,1)) # add na again...
	names(res) <- c(1:length(ticks), 'NA')
	return(res)
}

# prob that came from group_by_break
get_max_percent <- function(prob, limit) {
	res <- list(conf_level=NULL, prob_idx=NULL)
	prob_no_na <- prob[1:(length(prob)-1)] # just remove NA value
	len <- length(prob_no_na)
	max_percent <- 0
	min_interval <- Inf
	for (i in seq_along(prob_no_na)) {
		for (j in 1:(len-i+1)) {
			interval <- j:(j+i-1)
			percent <- sum(prob_no_na[interval])
			n <- length(interval)
			if (percent >= limit
				 	&& percent > max_percent
					&& (min_interval == Inf || n <= min_interval)) {
				min_interval <- length(interval)
				res$conf_level <- max_percent <- percent
				res$prob_idx <- interval
			}
		}
	}
	return(res)
}

# TODO this probably doesn't make sense, think about removing it
# receives a tree modified by sfreemap
sfreemap.pie_plot <- function(tree, percent_only=FALSE) {
    get_rows <- function(x) {
        y <- x$mapped.edge
        rownames(y) <- x$edge[,1]
        y <- y[as.character(length(x$tip)+1:x$Nnode),]
        return(y)
    }

    get_percentage <- function(row) {
        return(row/sum(row))
    }

    do_the_plot <- function(percent) {
        if (class(tree) == 'multiPhylo') {
            t <- tree[[1]]
        } else {
            t <- tree
        }
        states <- colnames(t$Q)

        if (is.null(t$node.label)) {
            t$node.label <- unique(t$edge[,1])
        }

        l_values <- colnames(percent) <- states
        l_colors <- palette()[1:length(l_values)]

        plot.phylo(t, no.margin=TRUE
                    , show.tip.label=TRUE, show.node.label=TRUE
                    , label.offset=0.02)
        nodelabels(pie=percent, piecol=l_colors, cex=0.6)
        legend('topright', legend=l_values, text.col=l_colors)
    }

    if (class(tree) == 'multiPhylo') {
        # get the percentage for each tree
        all_percent <- lapply(tree, sfreemap.pie_plot, percent_only=TRUE)
        # get a mean of all percentages
        percent <- Reduce('+', all_percent)/length(all_percent)
    } else if (class(tree) == 'phylo') {
        if (is.null(tree$mapped.edge)) {
            stop("tree should contain mapped states on edges")
        }
        percent <- t(apply(get_rows(tree), 1, get_percentage))
    } else {
        stop('tree should be object of class \"phylo\"')
    }

    if (isTRUE(percent_only)) {
        return (percent)
    } else {
        do_the_plot(percent)
    }
}
