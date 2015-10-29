# This is a class to plot correlation between mapping executions.
# Works for sfreemap, simmap and any other method that returns a phylo type tree
# with mapped.edge matrix
# Based on the plot found in http://www.inside-r.org/r-doc/graphics/pairs
sfreemap.correlation <- function(state) {
    data <- list(
        state = state
        , data = NULL
        , color = NULL
    )

    class(data) <- append(class(data), 'sfreemap.correlation')
    return(data)
}

# This method adds more data to the object obj.
# map is a sfreemap.map object, name is the name that will be shown in the plot
add.data <- function(obj, map, name) {
    UseMethod('add.data', obj)
}
add.data.sfreemap.correlation <- function(obj, map, name) {

    data <- matrix(map$mapped.edge[,obj$state], ncol=1)
    colnames(data) <- name
    rownames(obj$data) <- NULL

    if (is.null(obj$data)) {
        obj$data <- data
    } else {
        obj$data <- cbind(obj$data, data)
    }

    if (is.null(obj$color)) {
        # add color information
        nodes <- rownames(map$mapped.edge)
        nodes <- sapply(strsplit(nodes, ','), function(x) x[2])
        tip_label <- color <- map$tip.label[as.numeric(nodes)]
        color[is.na(tip_label)] <- 'grey'
        color[!is.na(tip_label)] <- 'green'
        obj$color <- color
    }

    return(obj)
}

# Actualy plots the correlation matrix
# TODO: add title, axis labels and stuff..
plot <- function(obj) {
    UseMethod('plot', obj)
}
plot.sfreemap.correlation <- function(obj) {

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

    # Transform to data.frame
    data <- data.frame(obj$data, stringsAsFactors=FALSE)

    # build formula
    formula <- paste('~', paste(colnames(data), collapse='+'))
    formula <- as.formula(formula)

    # package car has a scatterplotMatrix function that can print things like
    # an ellipse showing 95% confidence level, but I'm not sure it this fits
    # here. https://cran.r-project.org/web/packages/car/car.pdf
    pairs(formula
        , data=data
        , pch = 21
        , bg = obj$color
        , upper.panel=panel.cor
        , diag.panel=panel.hist
    )
}
