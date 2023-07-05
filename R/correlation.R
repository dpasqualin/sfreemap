# This is a class to plot correlation between mapping executions.
# Works for sfreemap, simmap and any other method that returns a phylo type tree
# with mapped.edge matrix
# Based on the plot found in http://www.inside-r.org/r-doc/graphics/pairs
correlation <- function(map, state, name) {

    # sanity check
    if (!inherits(map, 'phylo')) {
        stop ('map object should be of class "phylo"')
    }

    if (!state %in% colnames(map$mapped.edge)) {
        stop ('Unrecognized state', state)
    }

    # get data
    data <- matrix(map$mapped.edge[,state], ncol=1)
    colnames(data) <- name
    rownames(data) <- NULL

    # add tip information
    nodes <- rownames(map$mapped.edge)
    nodes <- sapply(strsplit(nodes, ','), function(x) x[2])
    tip_label <- map$tip.label[as.numeric(nodes)]
    is_tip <- rep(FALSE, length(tip_label))
    is_tip[!is.na(tip_label)] <- TRUE

    obj <- list(
        data = data
        , is_tip = is_tip
    )

    class(obj) <- append('correlation', class(obj))

    return(obj)
}

# Allow to sum up two correlation objects
"+.correlation" <- function(c1, c2) {
    obj <- c1

    obj$data <- cbind(obj$data, c2$data)
    # it is a tip if node is a tip in any of the mapping analysed
    obj$is_tip <- obj$is_tip | c2$is_tip

    return(obj)
}

# print correlation object as a list
print.correlation <- function(x, ...) {
    NextMethod()
}

# Plots the correlation matrix
plot.correlation <- function(x, y=NULL, ...) {

    # panel.smooth function is built in.
    # panel.cor puts correlation in upper panels, size proportional to correlation
    # from: http://www.r-bloggers.com/scatterplot-matrices-in-r/
    panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)
    {
        usr <- par("usr"); on.exit(par(usr))
        par(usr = c(0, 1, 0, 1))
        r <- abs(cor(x, y))
        txt <- format(c(r, 0.123456789), digits=digits)[1]
        txt <- paste(prefix, txt, sep="")
        if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
        text(0.5, 0.5, txt, cex = cex.cor * r)
    }

    # add histogram on the diagonals
    # from: https://stat.ethz.ch/R-manual/R-devel/library/graphics/html/pairs.html
    panel.hist <- function(x, ...) {
        usr <- par("usr")
        on.exit(par(usr))
        par(usr = c(usr[1:2], 0, 1.5))
        h <- hist(x, plot = FALSE)
        breaks <- h$breaks
        nB <- length(breaks)
        y <- h$counts
        y <- y/max(y)
        rect(breaks[-nB], 0, breaks[-1], y, col = "grey", ...)
    }

    # correlation object
    cor <- x

    # Transform to data.frame
    data <- data.frame(cor$data, stringsAsFactors=FALSE)

    # build formula
    formula <- paste('~', paste(colnames(data), collapse='+'))
    formula <- as.formula(formula)

    # define color, green when it is a tip when grey when isn't
    color <- ifelse(cor$is_tip, 'green', 'grey')

    # package car has a scatter plot Matrix function that can print things like
    # an ellipse showing 95% confidence level, but I'm not sure it this fits
    # here. https://cran.r-project.org/web/packages/car/car.pdf
    pairs(formula
        , data=data
        , pch = 21
        , bg = color
        , upper.panel=panel.cor
        , diag.panel=panel.hist
    )
}
