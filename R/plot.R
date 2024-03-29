plot_distribution_chart <- function(map, nodes=NULL, trees=NULL, states=NULL, conf_level=95
	                                   , number_of_ticks=20, type='emr') {

    if (!type %in% c('lmt', 'emr', 'mr')) {
        stop(paste('Unrecognized type', type))
    }

    # Filtering data for type and trees
    map <- map[[type]]
    if (!is.null(trees)) {
        map <- map[trees,,]
    }

    # check whether states exists
    if (!is.null(states) && !all(states %in% names(map[1,,1]))) {
        stop(paste('Unrecognized state:', states))
    }

    # Filtering by nodes
    # NOTE: this is a bit of a hack here that gives an object that doesn't make
    # much sense by it's own, but works nicely in this function
    if (is.null(nodes) && inherits(map, "array")) {
        nodes <- apply(map, 2, c)
    } else if (inherits(map, "array")) {
        nodes <- apply(map[,,nodes], 2, c)
    } else {
        nodes <- t(map)
    }

	# all states, or states passed as argument
	if (is.null(states)) {
		states <- colnames(nodes)
	} else {
	 	states <- c(states)
	}

	# first divide the dataset
	ticks <- get_ticks(map, type, number_of_ticks)

	to_plot <- data.frame(x=ticks, alpha=rep(FALSE, length(ticks)))
	for (state in states) {
		data <- get_state_data(nodes, state, conf_level, ticks)

		to_plot[[as.character(state)]] <- data$final_prob
		to_plot$alpha[data$final_idx] <- TRUE
	}

	# branch posterior probability
	# when NA = 0, bpp = 100%
	bpp <- 100-data$final_na_percent

	# graph config
	title <- "Posterior Distribution of Branch Lengths"
	subtitle <- paste("Branch posterior probability: ", bpp, "%", sep="")
	ylabel <- "Probability"
    # label angle
    angle <- 0
	if (type == 'emr') {
		xlabel <- "Dwelling time (% of branch length)"
		# shift bars do it fits within intervals
        to_plot$x <- to_plot$x + 2.5
	} else if (type == 'lmt') {
		angle <- 90
		xlabel <- "Expected number of state transitions"
	} else if (type == 'mr') {
        angle <- 90
        xlabel <- "Mutation rate"
    }

	melted <- melt(to_plot, id=c('x', 'alpha'))

	# format x labels, limiting the number of digits
	fmt <- function(){
	  function(x) format(x, digits=5)
	}

	p <- ggplot(melted, aes_string(x='x', y='value', fill='variable')) +
			# define the alpha for bars inside and outside HPD
	 		scale_alpha_discrete(range=c(0.3, 0.6), guide=FALSE) +
			# x+2.5 ensures that the beginning of the bar will be at the
			# beginning of the interval, and not at the middle
			geom_bar(stat="identity", position="identity", aes_string(alpha='alpha')) +
			# x and y labels
			xlab(xlabel) + ylab(ylabel) +
			# add title and subtitle
			ggtitle(bquote(atop(.(title), atop(italic(.(subtitle)), "")))) +
			# legend name
			scale_fill_discrete(name = "States") +
			# set the breaks (ticks) at x axis
			scale_x_continuous(breaks=ticks, labels=fmt()) +
			# put legend at the bottom
			theme(legend.position="bottom", text=element_text(size=15)
					, axis.text.x=element_text(angle=angle))

	return(p)
}

plot_distribution_tree <- function(map, state='all', type='emr'
                                , conf_level=95, number_of_ticks=20
                                , tip_states=NULL, cex=0.9, ftype="i"
                                , lwd=3, tip.label=NULL) {

    if (!type %in% c('lmt', 'emr', 'mr')) {
        stop(paste('Unrecognized type', type))
    }

    tree <- map$base_tree
    map <- map[[type]]
	ticks <- get_ticks(map, type, number_of_ticks)
	tree$maps <- list()

    if (!state %in% c('all', names(map[1,,1]))) {
        stop ('Unrecognized state')
    }

    # all but the root node
    all_nodes <- unique(tree$edge[,2])

    if (type %in% c('mr', 'lmt')) {
        color_names <- format(round(ticks,1), nsmall=1, trim=TRUE, scientific=FALSE)
    } else if (type == 'emr') {
        color_names <- as.character(ticks)
    }

	for (node in all_nodes) {

		data <- get_state_data(map[,,node], state, conf_level, ticks, na.rm=FALSE)

		if (data$final_conf_level >= conf_level) {
			value <- freq_to_prob(data$final_prob[data$final_idx])
			# workaround to set <NA> as a character
			tmp <- as.character(as.numeric(names(value)))
			tmp[is.na(tmp)] <- 'NA'
			names(value) <- tmp

            # we want to show the more relevant probabilities at the end of the
            # branch. The more relevant probabilities are the ones either closer
            # to zero or to one hundred percent
            tmp <- as.numeric(tmp)
            names(value) <- color_names[data$final_idx]
            if (tmp[1] < 100 - tail(tmp, n=1)) {
               value <- value[order(tmp, decreasing=TRUE)]
            }
		} else {
			# 100% unknown
			value <- 100
			names(value) <- 'NA'
		}

		# scale maps to branch length
        b_number <- which(tree$edge[,2]==node)
		b_len <- tree$edge.length[b_number]
		tree$maps[[b_number]] <- (value*b_len)/100.0
    }

	# color gradient
	colors <- get_color_pallete(color_names)

	# make room for the legend
	ylim <- c(-2, length(tree$tip.label))

    # add tip state
    if (!is.null(tip_states)) {
        tree$tip.label <- join_tip_states(tree, tip_states)
    }

    if (!is.null(tip.label)) {
        if (length(tip.label) == length(tree$tip.label)) {
            tree$tip.label <- tip.label
        } else {
            stop("tip.label doesn't match the number of tips")
        }
    }

	plotSimmap(tree, colors=colors, fsize=cex, ftype=ftype, ylim=ylim, lwd=lwd)

	add_subtitle(colors=colors, cex=cex)

	return(tree)
}

join_tip_states <- function(tree, tip_states) {
    # Define the tip states as a matrix
    if (!is.matrix(tip_states)) {
        tip_states <- build_states_matrix(tree$tip.label, tip_states)
    }

    # FIXME: this is a hell of a work around, there should be a better way
    # to do it
    a <- t(apply(tip_states, 1, function(x) ifelse(x==1, names(x), '')))
    b <- apply(a, 1, paste, collapse='')
    b <- b[tree$tip.label]
    res <- apply(cbind(names(b), b), 1, paste, collapse=' - ')
    names(res) <- NULL
    return (res)
}

get_color_pallete <- function(color_names) {
    red <- c(1,1,2,2,3,4,5,7,9,12,16,21,28,37,48,64,84,111,147,193,255)
    green <- c(0,48,92,130,163,191,214,232,245,252,255,252,245,232,214,191,163,130,92,48,0)
    blue <- c(255,193,147,111,84,64,48,37,28,21,16,12,9,7,5,4,3,2,2,1,1)

    colors <- rgb(red, green, blue, maxColorValue=255)

    # set color names
    names(colors) <- color_names
    colors['NA'] <- '#B3B3B3FF' # grey 30%

    return (colors)
}

get_state_data <- function(node, state, conf_level, ticks, na.rm=TRUE) {
    if (state == 'all') {
        data <- apply(node, 1, sum)
    } else {
    	data <- node[,as.character(state)]
    }

	final_conf_level <- 0.0
	final_idx <- c()
	final_na_percent <- 0.0

	# group values into the closest ticknames(b) <- NULL
	group <- gen_group(data, ticks)

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
	final_prob <- prob
	if (isTRUE(na.rm)) {
		# remove NA values from final_prob
		final_prob <- prob[1:(length(prob)-1)]
	}

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
	names(res) <- c(ticks, 'NA')
	return(res)
}

# prob that came from group_by_break
get_max_percent <- function(prob, limit) {
	res <- list(conf_level=0, prob_idx=c())
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

get_ticks <- function(data, type, number_of_ticks) {
	# first divide the dataset
	if (type == 'emr') {
		ticks <- seq(0, 100, 100/number_of_ticks)
	} else if (type %in% c('lmt', 'mr')) {
		begin <- min(data, na.rm=TRUE)
		end <- max(data, na.rm=TRUE)
		ticks <- seq(begin, end, (end-begin)/number_of_ticks)
	} else {
		stop(paste('Unrecognized type:', type))
	}
	return(ticks)
}
