# adds legend to an open stochastic map style plot
# modified from phytools (written by Liam J. Revell 2013)
sfreemap.add.legend <- function(leg=NULL, colors, prompt=FALSE
                                , vertical=FALSE, ...) {

    shape <- "rect"
    if (hasArg(shape)) {
      shape <- list(...)$shape
    }

    if (hasArg(fsize)) {
        fsize <- list(...)$fsize
    } else {
        fsize <- ifelse(isTRUE(vertical), 0.7, 1)
    }

    h <- fsize * strheight(LETTERS[1])
    usr <- par()$usr
    w <- h*(usr[2]-usr[1]) / (usr[4]-usr[3]) * 1.5

    if (prompt) {
        cat("Click where you want to draw the legend\n")
        pos <- unlist(locator(1))
        x <- pos[1]
        y <- pos[2]
    } else {
        x <- ifelse(hasArg(x), list(...)$x, -4)
        y <- ifelse(hasArg(y), list(...)$y, -4)
    }

    if (is.null(leg)) {
      leg <- names(colors)
    }

    if (vertical) {
        y <- y+0:(length(leg)-1)*1.5*h
        x <- rep(x,length(y))
        text(x+w/2, y+0.5, leg, pos=4, cex=fsize/par()$cex)
    } else {
        x <- x + (0:(length(leg)-2) * w)
        # add NA a bit separate from the rest
        x <- c(x, w*length(leg))
        y <- rep(y,length(x))

        labels <- leg

        text(x, y+0.5, labels, pos=3, cex=0.7*fsize)
    }

    if (shape=="rect") {
        rects <- cbind(rep(w, length(x)), rep(h*1, length(x)))
        y <- y + h/2
        symbols(x, y, rectangles=rects, bg=colors, add=TRUE, inches=FALSE)
    } else if (shape=="square") {
        symbols(x, y, squares=rep(w,length(x)), bg=colors, add=TRUE, inches=FALSE)
    } else if (shape=="circle") {
        symbols(x, y, circles=rep(w,length(x)), bg=colors, add=TRUE, inches=FALSE)
    } else {
        stop(paste("shape=\"",shape,"\" is not a recognized option.",sep=""))
    }
}

sfreemap.read_tips <- function(file, character=1, sep="\t") {

    data <- read.csv(file, sep=sep, header=FALSE, colClasses = "character"
                         , fill=FALSE, strip.white = TRUE)

    # data representing the state or states (when ambiguous) in which every taxa # is in
    taxa_state <- data[,character+1]             # first column is the label
    taxa_labels <- names(taxa_state) <- data[,1] # use taxa as state name
    # we have to split ambiguous states to find out how many there are
    possible_states <- paste(taxa_state, collapse='')
    possible_states <- unique(strsplit(possible_states, split='')[[1]])

    # check whether the taxon is in state
    taxon_in_state <- Vectorize(function(taxon, state) {
        res <- if (grepl(state, taxa_state[taxon])) 1 else 0
        return (res)
    })

    res <- outer(taxa_labels, possible_states, taxon_in_state)
    rownames(res) <- taxa_labels
    colnames(res) <- possible_states

    return (res)
}

# return t1 only with tips that are in t2 too
# optionally reroot at node 'reroot'
sfreemap.pruning <- function(t1, t2, reroot=NULL) {
    tips_to_remove <- t1$tip.label[!t1$tip.label %in% t2$tip.label]
    t <- drop.tip(t1, tips_to_remove)
    if (!is.null(reroot)) {
        if (reroot %in% t$tip.label) {
            t <- root(t, reroot)
        } else {
            msg <- paste('trying to root tree but tip', reroot, 'doesn\'t exist')
            stop(msg)
        }
    }
    return(t)
}

describe.sfreemap <- function (tree, ...) {
    if (inherits(tree, "phylo")) {
        lmt <- colSums(tree$mapped.edge.lmt)
        emr <- colSums(tree$mapped.edge)
    } else if (inherits(tree, "multiPhylo")) {
        lmt <- t(sapply(tree, function(x) colSums(x$mapped.edge.lmt)))
        lmt <- colMeans(lmt)

        emr <- t(sapply(tree, function(x) colSums(x$mapped.edge)))
        emr <- colMeans(emr)
    } else {
        stop ("tree must be an object of type \"phylo\" or \"multiPhylo\"")
    }

    emr_total <- sum(emr)
    names(emr_total) <- 'total'
    emr <- c(emr, emr_total)
    emr <- rbind(emr, prop=emr/emr_total)

    return (list(transitions=lmt, dwelling_times=emr))
}

# Alias for describe.sfreemap
summary.sfreemap <- function (object, ...) {
    describe.sfreemap(object, ...)
}

# function reorders sfreemap tree
# based on reorderSimmap, written by Liam Revell 2011, 2013
reorder.sfreemap <- function(x, ...) {

    tree <- x

    order <- "cladewise"
    if (hasArg(order)) {
        order <- list(...)$order
    }

    index.only <- FALSE
    if (hasArg(index.only)) {
        index.only <- list(...)$index.only
    }

    index <- reorder.phylo(tree, order, index.only=TRUE, ...)
    if (!index.only) {
        if (inherits(index, "phylo")) {
            ## bug workaround, from phytools package
            index <- whichorder(index$edge[,2], tree$edge[,2])
        }
        tree$edge <- tree$edge[index,]
        tree$edge.length <- tree$edge.length[index]
        if (!is.null(tree$mapped.edge)) {
            tree$mapped.edge <- tree$mapped.edge[index,]
            tree$mapped.edge.lmt <- tree$mapped.edge.lmt[index,]
        }
        attr(tree, "order") <- order
        return (tree)
    } else {
        return (index)
    }
}

# function whichorder
# written by Liam Revell 2011, 2013
whichorder <- function(x,y) sapply(x,function(x,y) which(x==y),y=y)

freq_to_prob <- function(x) {

    if (is.null(dim(x))) {
        # vector
        result <- (x*100)/sum(x)
    } else if (length(dim(x)) == 2) {
        # matrix
        result <- (x*100)/rowSums(x)
    } else if (length(dim(x)) == 3) {
        # 3-dimentional array
        result <- x
        for (i in 1:dim(x)[3]) {
            tmp <- x[,,i]
            result[,,i] <- (tmp*100) / rowSums(tmp)
        }
    } else {
        stop("dim(x) cannot be greater than 3")
    }
    return (result)
}

sfreemap.read.fasta <- function(file, ensure_dna=FALSE) {
    nucleo_data <- read.fasta(file)
    nucleo_data <- sapply(nucleo_data, cbind)
    if (ensure_dna) {
        # converts RNA to DNA
        nucleo_data[nucleo_data == 'u'] <- 't'
    }
    rownames(nucleo_data) <- paste('t', 1:nrow(nucleo_data), sep='')
    return (t(nucleo_data))
}

## based on function to rescale simmap style trees
## written by Liam J. Revell 2012, 2013, 2014, 2015
sfreemap.rescale <- function(tree, height=NULL, parallel=FALSE) {
    if (inherits(tree, "multiPhylo")){
        tree <- unclass(tree)
        if (parallel==TRUE && !on_windows()) {
            tree <- mclapply(tree, sfreemap.rescale, height, mc.cores=detectCores())
        } else {
            tree <- lapply(tree, sfreemap.rescale, height)
        }
        class(tree) <- "multiPhylo"
        return (tree)
    } else if (inherits(tree, "phylo")) {
        max_height <- max(nodeHeights(tree))
        if (is.null(height)) {
            height <- max_height
        }
        if (height != max_height) {
            s <- height/max_height
            tree$edge.length <- tree$edge.length * s
            if (inherits(tree, "sfreemap")) {
                tree$mapped.edge <- tree$mapped.edge * s
                tree$mapped.edge.lmt <- tree$mapped.edge.lmt * s
            }
        }
        return(tree)
    } else {
        message("tree should be an object of class \"phylo\" or \"multiPhylo\"")
    }
}

# Return true if we are running on windows and false otherwise
on_windows <- function() {
    return(Sys.info()['sysname'] == 'Windows')
}

# This function creates a matrix with rownames being states, colnames being
# the tip labels and values being 1, if the tip label has the correspondent
# states and 0 otherwise.
# "possible_states" is useful mainly because DNA can always have A, C, T and G,
# but data might not have one of this letters.
build_states_matrix <- function(tip_labels, tip_states, possible_states=NULL) {
    if (is.null(possible_states)) {
        possible_states <- tip_states
    }
    res <- to.matrix(tip_states, sort(unique(possible_states)))
    return (res[tip_labels,])
}
