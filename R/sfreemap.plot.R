# TODO: look for p-value on the internet, probably this entire function
# is already implemented somewhere
sfreemap.plot_distribution <- function(node, states=NULL, conf_level=90
	, number_of_ticks=20, type='emr', ...) {

	# TODO: add sanity check for parameters

	# all states or states passed as argument
	if (is.null(states)) {
		states <- colnames(node)
	} else {
	 	states <- c(states)
	}

	# first divide the dataset
	ticks <- get_ticks(node, type, number_of_ticks)

	to_plot <- data.frame(x=ticks, alpha=rep(FALSE, length(ticks)))
	for (cont in 1:length(states)) {
		state <- states[cont]
		data <- get_state_data(node, state, conf_level, ticks)

		to_plot[[as.character(state)]] <- data$final_prob
		to_plot$alpha[data$final_idx] <- TRUE
	}
	print(to_plot)

	# branch posterior probability
	# when NA = 0, bpp = 100%
	bpp <- 100-data$final_na_percent

	# graph config
	title <- "Posterior Distribution of Branch Lengths"
	subtitle <- paste("Branch posterior probability: ", bpp, "%", sep="")
	ylabel <- "Probability"
	if (type == 'emr') {
		xlabel <- "Dwelling time (% of branch length)"
		# shift bars do it fits within intervals
		to_plot$x <- to_plot$x + 2.5
		# label angle
		angle <- 0
	} else if (type == 'lmt') {
		# label angle
		angle <- 90
		xlabel <- "Expected number of state transitions"
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

sfreemap.plot_tree <- function(base_tree, trees, state, type='emr'
															, conf_level=90, number_of_ticks=20) {

	map <- sfreemap.map_posterior_distribution(base_tree, trees, scale=TRUE)
	ticks <- get_ticks(node, type, number_of_ticks)

	tree <- base_tree
	tree$maps <- list()
	for (node in 1:nrow(tree$edge)) {
		b_len <- tree$edge.length[node]

		data <- get_state_data(map$emr[,,node], state, conf_level, ticks, na.rm=FALSE)

		if (data$final_conf_level >= conf_level) {
			value <- freq_to_prob(data$final_prob[data$final_idx])
			# workaround to set <NA> as a character
			tmp <- as.character((as.numeric(names(value))-1) * 5)
			tmp[is.na(tmp)] <- 'NA'
			names(value) <- tmp
		} else {
			# 100% unknown
			value <- 100
			names(value) <- 'NA'
		}

		# scale maps to branch length
		tree$maps[[node]] <- (value*b_len)/100.0
	}

	# color grandient
	colfunc <- colorRampPalette(c("#A2A2FF", "blue"))
	colors <- colfunc(length(ticks))

	color_names <- as.numeric(unique(unlist(sapply(tree$maps,names))))
	color_names <- as.character(sort(color_names, na.last=TRUE))
	color_names[is.na(color_names)] <- 'NA'

	names(colors) <- color_names
	colors['NA'] <- 'red'

print(sort(unique(unlist(sapply(tree$maps,names)))))

	# make room for the legend
	ylim <- c(-4, length(tree$tip.label))

	plotSimmap(tree, colors=colors, fsize=0.7, ftype="i", ylim=ylim)
	sfreemap.add.legend(colors=colors)

	return(tree)
}

get_state_data <- function(node, state, conf_level, ticks, na.rm=TRUE) {
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
	names(res) <- c(1:length(ticks), 'NA')
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

get_ticks <- function(node, type, number_of_ticks) {
	# first divide the dataset
	if (type == 'emr') {
		ticks <- seq(0, 100, 100/number_of_ticks)
	} else if (type == 'lmt') {
		begin <- min(node)
		end <- max(node)
		ticks <- seq(begin, end, (end-begin)/number_of_ticks)
	} else {
		stop(paste('unrecognized type:', type))
	}
	return(ticks)
}
