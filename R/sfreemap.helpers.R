# THESE ARE MOSTLY EXTERNAL FUNCTIONS HONESTLY STOLEN AND
# SLIGHTLY CHANGE FROM PHYTOOLS. WE HAD TO DO THIS BECAUSE
# THIS FUNCTIONS ARE NOT VISIBLE FROM OUTSIDE PHYTOOLS

# PHYTOOLS
# get pars
# written by Liam J. Revell 2013
getPars<-function(bt,xx,model,Q,tree,tol,m,omp,liks=TRUE){
    XX<-apeAce(bt,xx,model,omp,fixedQ=Q,output.liks=liks)
    N<-length(bt$tip.label)
    II<-XX$index.matrix
    lvls<-XX$states
    if(liks){
        L<-XX$lik.anc
        rownames(L)<-N+1:nrow(L)
        if(!is.binary.tree(tree)){
            ancNames<-matchNodes(tree,bt)
            L<-L[as.character(ancNames[,2]),]
            rownames(L)<-ancNames[,1]
        }
        L<-rbind(xx,L)
        rownames(L)[1:N]<-1:N
    } else L<-NULL
    if(any(XX$rates<tol)){
        #message(paste("\nWarning: some elements of Q not numerically distinct from 0; setting to",tol,"\n"))
        XX$rates[XX$rates<tol]<-tol
    }
    Q<-matrix(XX$rates[II],m,m,dimnames=list(lvls,lvls))
    diag(Q)<--rowSums(Q,na.rm=TRUE)
    return(list(Q=Q,L=L,loglik=XX$loglik))
}

# PHYTOOLS
# function uses numerical optimization to solve for the stationary distribution
# written by Liam J. Revell 2013
statdist <- function(Q){
    foo <- function(theta,Q){
        Pi<-c(theta[1:(nrow(Q)-1)],1-sum(theta[1:(nrow(Q)-1)]))
        sum((Pi%*%Q)^2)
    }
    k<-nrow(Q)
    if (nrow(Q)>2) {
        fit <- optim(rep(1/k,k-1), foo, Q=Q, control=list(reltol=1e-16))
        return(setNames(c(fit$par[1:(k-1)], 1-sum(fit$par[1:(k-1)])), rownames(Q)))
    } else {
        fit<-optimize(foo,interval=c(0,1),Q=Q)
        return(setNames(c(fit$minimum,1-fit$minimum),rownames(Q)))
    }
}

# PHYTOOLS
# function for conditional likelihoods at nodes, from ace(...,type="discrete")
# modified (only very slightly) from E. Paradis et al. 2013
apeAce <- function(tree,x,model,omp,fixedQ=NULL,...){
    if(hasArg(output.liks)) output.liks<-list(...)$output.liks
    else output.liks<-TRUE
    ip<-0.1
    nb.tip<-length(tree$tip.label)
    nb.node<-tree$Nnode
    if(is.matrix(x)){
        x<-x[tree$tip.label,]
        nl<-ncol(x)
        lvls<-colnames(x)
    } else {
        x<-x[tree$tip.label]
        if(!is.factor(x)) x<-factor(x)
        nl<-nlevels(x)
        lvls<-levels(x)
        x<-as.integer(x)
    }
    if(is.null(fixedQ)){
        if(is.character(model)){
            rate<-matrix(NA,nl,nl)
            if(model=="ER") np<-rate[]<-1
            if(model=="ARD"){
                np<-nl*(nl-1)
                rate[col(rate)!=row(rate)]<-1:np
            }
            if (model=="SYM") {
                np<-nl*(nl-1)/2
                sel<-col(rate)<row(rate)
                rate[sel]<-1:np
                rate<-t(rate)
                rate[sel]<-1:np
            }
        } else {
            if(ncol(model)!=nrow(model)) stop("the matrix given as 'model' is not square")
            if(ncol(model)!=nl) stop("the matrix 'model' must have as many rows as the number of categories in 'x'")
            rate<-model
            np<-max(rate)
        }
        Q<-matrix(0,nl,nl)
    } else {
        rate<-matrix(NA,nl,nl)
        np<-nl*(nl-1)
        rate[col(rate)!=row(rate)]<-1:np
        Q<-fixedQ
    }
    index.matrix<-rate
    tmp<-cbind(1:nl,1:nl)
    index.matrix[tmp]<-NA
    rate[tmp]<-0
    rate[rate==0]<-np+1
    liks<-matrix(0,nb.tip+nb.node,nl)
    TIPS<-1:nb.tip
    if(is.matrix(x)) liks[TIPS,]<-x
    else liks[cbind(TIPS,x)]<-1
    phy<-reorder(tree,"pruningwise")

    dev<-function(p,output.liks=FALSE,fixedQ=NULL){
        if(any(is.nan(p))||any(is.infinite(p))) return(1e50)
        comp<-numeric(nb.tip+nb.node)
        if(is.null(fixedQ)){
            Q[]<-c(p,0)[rate]
            diag(Q)<--rowSums(Q)
        } else Q<-fixedQ

        Q_eigen <- eigen(Q, TRUE, only.values = FALSE)
        Q_eigen[['vectors_inv']] <- solve(Q_eigen$vectors)
        tb <- transition_probabilities(Q_eigen, phy$edge.length, omp)

        for(i in seq(from=1,by=2,length.out=nb.node)){
            j<-i+1L
            anc<-phy$edge[i,1]
            des1<-phy$edge[i,2]
            des2<-phy$edge[j,2]
            v.l<-tb[,,i]%*%liks[des1,]
            v.r<-tb[,,j]%*%liks[des2,]
            v<-v.l*v.r
            comp[anc]<-sum(v)
            liks[anc,]<-v/comp[anc]
        }
        if(output.liks) return(liks[-TIPS,])
        dev<--2*sum(log(comp[-TIPS]))
        if(is.na(dev)) Inf else dev
    }
    if(is.null(fixedQ)){
        out<-nlminb(rep(ip,length.out=np),function(p) dev(p),lower=rep(0,np),upper=rep(1e50,np))
        obj<-list()
        obj$loglik<--out$objective/2
        obj$rates<-out$par
        obj$index.matrix<-index.matrix
        if(output.liks){
            obj$lik.anc<-dev(obj$rates,TRUE)
            colnames(obj$lik.anc)<-lvls
        }
        obj$states<-lvls
    } else {
        out<-dev(rep(ip,length.out=np),fixedQ=Q)
        obj<-list()
        obj$loglik<--out/2
        obj$rates<-fixedQ[sapply(1:np,function(x,y) which(x==y),index.matrix)]
        obj$index.matrix<-index.matrix
        if(output.liks){
            obj$lik.anc<-dev(obj$rates,TRUE,fixedQ=Q)
            colnames(obj$lik.anc)<-lvls
        }
        obj$states<-lvls
    }
    return(obj)
}

# PHYTOOLS
# mcmc for Q used in Q="mcmc"
# written by Liam J. Revell 2013
mcmcQ<-function(bt,xx,model,tree,tol,m,burnin,samplefreq,nsim,vQ,prior){
    update<-function(x){
        ## x<-exp(log(x)+rnorm(n=np,mean=0,sd=sqrt(vQ)))
        x<-abs(x+rnorm(n=np,mean=0,sd=sqrt(vQ)))
        return(x)
    }
    # get model matrix
    if(is.character(model)){
        rate<-matrix(NA,m,m)
        if(model=="ER"){
            np<-rate[]<-1
            diag(rate)<-NA
        }
        if(model=="ARD"){
            np<-m*(m-1)
            rate[col(rate)!=row(rate)]<-1:np
        }
        if (model=="SYM") {
            np<-m*(m-1)/2
            sel<-col(rate)<row(rate)
            rate[sel]<-1:np
            rate<-t(rate)
            rate[sel]<-1:np
        }
    } else {
        if(ncol(model)!=nrow(model)) stop("the matrix given as 'model' is not square")
        if(ncol(model)!=m) stop("the matrix 'model' must have as many rows as the number of categories in 'x'")
        rate<-model
        np<-max(rate)
    }
    # burn-in
    p<-rgamma(np,prior$alpha,prior$beta)
    Q<-matrix(p[rate],m,m)
    diag(Q)<--rowSums(Q,na.rm=TRUE)
    yy<-getPars(bt,xx,model,Q,tree,tol,m)
    #cat("Running MCMC burn-in. Please wait....\n")
    for(i in 1:burnin){
        pp<-update(p)
        Qp<-matrix(pp[rate],m,m)
        diag(Qp)<--rowSums(Qp,na.rm=TRUE)
        zz<-getPars(bt,xx,model,Qp,tree,tol,m,FALSE)
        p.odds<-exp(zz$loglik+sum(dgamma(pp,prior$alpha,prior$beta,log=TRUE))-
            yy$loglik-sum(dgamma(p,prior$alpha,prior$beta,log=TRUE)))
        if(p.odds>=runif(n=1)){
            yy<-zz
            p<-pp
        }
    }
    # now run MCMC generation, sampling at samplefreq
    #cat(paste("Running",samplefreq*nsim,"generations of MCMC, sampling every",samplefreq,"generations. Please wait....\n"))
    XX<-vector("list",nsim)
    for(i in 1:(samplefreq*nsim)){
        pp<-update(p)
        Qp<-matrix(pp[rate],m,m)
        diag(Qp)<--rowSums(Qp,na.rm=TRUE)
        zz<-getPars(bt,xx,model,Qp,tree,tol,m,FALSE)
        p.odds<-exp(zz$loglik+sum(dgamma(pp,prior$alpha,prior$beta,log=TRUE))-
            yy$loglik-sum(dgamma(p,prior$alpha,prior$beta,log=TRUE)))
        if(p.odds>=runif(n=1)){
            yy<-zz
            p<-pp
        }
        if(i%%samplefreq==0){
            Qi<-matrix(p[rate],m,m)
            diag(Qi)<--rowSums(Qi,na.rm=TRUE)
            XX[[i/samplefreq]]<-getPars(bt,xx,model,Qi,tree,tol,m,TRUE)
        }
    }
    return(XX)
}

# adds legend to an open stochastic map style plot
# modified from phytools (written by Liam J. Revell 2013)
sfreemap.add.legend <- function(leg=NULL,colors,prompt=FALSE,vertical=FALSE,...) {
    if (hasArg(shape)) {
      shape <- list(...)$shape
    } else {
      shape <- "square"
    }

    if (prompt) {
      cat("Click where you want to draw the legend\n")
      pos <- unlist(locator(1))
      x <- pos[1]
      y <- pos[2]
    } else {
      x <- ifelse(hasArg(x), list(...)$x, 0)
      y <- ifelse(hasArg(y), list(...)$y, -4)
    }

    if (hasArg(fsize)) {
      fsize <- list(...)$fsize
    } else {
      fsize <- 1
    }

    if (is.null(leg)) {
      leg <- names(colors)
    }

    h <- fsize*strheight(LETTERS[1])

    usr <- par()$usr
    w <- h*(usr[2]-usr[1]) / (usr[4]-usr[3])

    if (vertical) {
      y <- y-0:(length(leg)-1)*1.5*h
      x <- rep(x+w/2,length(y))
      text(x+w, y, leg, pos=4, cex=fsize/par()$cex)
    } else {
      x <- x + (0:(length(leg)-2) * w)
      # add NA a bit separate from the rest
      x <- c(x, w*length(leg))
      y <- rep(y,length(x)) - h

      labels <- leg
      suppressWarnings(labels[as.numeric(labels) %% 10 != 0] <- '')

      text(x, y+0.5, labels, pos=3, cex=0.7*fsize)
    }

    if (shape=="square") {
      symbols(x,y,squares=rep(w,length(x)),bg=colors, add=TRUE,inches=FALSE)
    } else if (shape=="circle") {
      symbols(x,y,circles=rep(w,length(x)),bg=colors, add=TRUE, inches=FALSE)
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

sfreemap.describe <- function (tree) {
    if (class(tree) == 'phylo') {
        lmt <- colSums(tree$mapped.edge.lmt)
        emr <- colSums(tree$mapped.edge)
    } else if (class(tree) == 'multiPhylo') {
        lmt <- t(sapply(tree, function(x) colSums(x$mapped.edge.lmt)))
        lmt <- colMeans(lmt)

        emr <- t(sapply(tree, function(x) colSums(x$mapped.edge)))
        emr <- colMeans(emr)
    } else {
        stop ("tree must be an object of type \"phylo\" or \"multiPhylo\"")
    }
    return (list(transitions=lmt, dwelling_times=emr))
}

# function reorders simmap tree
# written Liam Revell 2011, 2013
sfreemap.reorder <- function(tree, order='cladewise') {
    x <- reorder(tree, order)
    o <- whichorder(x$edge[,2], tree$edge[,2])
    x$mapped.edge <- tree$mapped.edge[o,]
    # NOTE: make.simmap (phytools) use this, but I'm not sure why it's
    # useful.
    # x$maps <- tree$maps[o]
    return(x)
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
