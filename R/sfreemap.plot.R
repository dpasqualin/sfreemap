# TODO: look for p-value on the internet, probably this entire function
# is already implemented somewhere
sfreemap.plot_distribution <- function(node, state, conf.level=90, ...) {

	# Just a reminder of something we could do, which is to call this function
	# Starting from a specific confidence level that is not 100. This could be
	# useful do determine the highest possible fit
	#start.conf.level <- 100
	#if (hasArg(start.conf.level)) {
	#	start.conf.level <- list(...)$start.conf.level
	#}

	# TODO: add sanity check for parameters

	data <- node[,as.character(state)]

	# FIXME: this were arbitrarily defined and need to be reviwed
	min_ticks <- 10
	max_ticks <- length(data)/2

	begin <- min(data, na.rm=TRUE)
	end <- max(data, na.rm=TRUE)

	final_precision <- 0.0
	final_len <- Inf
	final_idx <- c()
	final_na_percent <- 0.0

	for (number_of_ticks in min_ticks:max_ticks) {
		# first divide the dataset
		tick <- (end - begin) / number_of_ticks
		ticks <- seq(begin, end, tick)
		# group values into the closest tick
		group <- gen_group(data, ticks)
		# just make sure all ticks have a value, even if it's zero

		# translate frequencies into probabilities of seem a particular values
		# in a tick
		prob <- freq_to_prob(group)

		# find the max precision we can get with probabilities calculated
		# the maximum precision we could find
		partial <- get_max_percent(prob, conf.level)
		# the indexes in prob objected that generated the precision above
		p_idx <- partial[['prob_idx']]

		if (is.null(p_idx)) {
			# we've found nothing with conf.level, let's try again
			# with different ticks
			next
		}
		# the number of indexes necessary for the precision will tell the
		# length necessary to get it
		total_len <- tick * length(p_idx)
		# update final values if we found a more suitable result
		if (total_len < final_len) {
			final_len <- total_len
			final_precision <- partial[['precision']]
			final_idx <- p_idx
			final_ticks <- ticks
			final_prob <- prob
		}
	}

	if (final_precision > 0) {
		# main plot
		p_x <- final_ticks
		p_y <- final_prob[1:(length(final_prob)-1)] # don't plot NA
		plot(p_x, p_y, type='l'
				, xlab="Branch length"
				, ylab="Probability"
				, main="Distribution of branch length across trees")

		# lines
		l_x <- final_ticks[c(min(final_idx), max(final_idx))]
		l_y <- rep(max(p_y), 2)

		# TODO: would be nice to have a shade inside the graph instead of a line
		# on top to show the interval found
		#pol_x <- rep(0, length(final_idx))
		#pol_y <- p_y[final_idx]
		#polygon(pol_x, pol_y, col='red')

		# add lines for confidence level
		par(pch=21, col="red")
		lines(l_x, l_y, type="o")

		# add grid
		g_x <- length(final_ticks)/10
		g_y <- length(final_prob)/10
		grid(g_x, g_y)

		# add NA percentage
		na_percent <- round(tail(final_prob,1), 2)
		mtext(paste('NA:', na_percent, '%'), 4)

		return(list(
			px=p_x
			, py=p_y
			, lx=l_x
			, ly=l_y
			, probabilities=final_prob
			, ticks=final_ticks
			, precision=final_precision))
	} else {
		return(NULL)
	}
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
	len <- length(prob)
	res <- list(precision=NULL, prob_idx=NULL)
	prob_no_na <- prob[1:(length(prob)-1)] # just remove NA value
	for (i in seq_along(prob_no_na)) {
		for (j in 1:(len-i+1)) {
			interval <- j:(j+i-1)
			percent <- sum(prob[interval])
			if (percent >= limit) {
				res$precision <- percent
				res$prob_idx <- interval
				return(res)
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
