# Return the transition matrix Q calculated empirically
Q_empirical <- function(tree, tip_states, prior, model, tol, omp) {

    tip_states <- build_states_matrix(tree$tip.label, tip_states)

    states <- tip_states/rowSums(tip_states)

    n_states <- ncol(states)

    # Reorder tree and create a copy called bt. Not sure why
    # phytools need this...
    new_tree <- bt <- reorder.phylo(tree, "cladewise")

    if (!ape::is.binary(bt)) {
        # NOTE: find out what this does...
        bt <- multi2di(bt)
    }

    # NOTE: Add comments to this snipped. It was adapted from phytools
    # and I need to understand what it is doing
    pars <- getPars(bt, states, model, Q=NULL, new_tree, tol, n_states, omp)
    L <- pars$L
    Q <- pars$Q
    logL <- pars$loglik

    # Set prior (pi)
    if (prior[1] == "equal") {
        # set equal
        prior <- setNames(rep(1/n_states, n_states), colnames(L))
    } else if (prior[1] == "estimated") {
        # set from stationary distribution
        prior <- statdist(Q)
    } else {
        # obtain from input
        prior <- prior/sum(prior)
    }

    return (list(Q=Q, prior=prior, logL=logL))

}

# Return the transition matrix Q if it has been passed as a matrix
Q_matrix <- function(tree, tip_states, Q, model, prior, tol, type) {

    if (type == "dna") {
        states <- c("a", "c", "t", "g", "-")
    } else {
        states <- NULL
    }

    tip_states <- build_states_matrix(tree$tip.label, tip_states, states)

    states <- tip_states/rowSums(tip_states)

    n_states <- ncol(states)

    # Reorder tree and create a copy called bt. Not sure why
    # phytools need this...
    new_tree <- bt <- reorder.phylo(tree, "cladewise")

    XX <- getPars(bt, states, model, Q=Q, new_tree, tol, n_states)
    # NOTE: do we need this?
    L <- XX$L
    logL <- XX$loglik

    # Set the priors
    if (prior[1] == "equal") {
        prior <- setNames(rep(1/n_states,n_states), colnames(L)) # set equal
    } else if (prior[1] == "estimated") {
        prior <- statdist(Q) # set from stationary distribution
    } else {
        prior <- prior/sum(prior) # obtain from input
    }

    return (list(Q=Q, prior=prior, logL=logL))
}

# Return the transition matrix Q calculated using a markov chain
Q_mcmc <- function(tree, tip_states, prior, model, gamma_prior, tol, burn_in, sample_freq, vQ, n_simulations, omp) {

    tip_states <- build_states_matrix(tree$tip.label, tip_states)

    states <- tip_states/rowSums(tip_states)

    n_states <- ncol(states)

    # Reorder tree and create a copy called bt. Not sure why phytools need this copy.
    # TODO: As we set the order of the tree in the main function, maybe we don't need
    # to do it here again.
    new_tree <- bt <- reorder.phylo(tree, "cladewise")

    # Define prior
    if (gamma_prior$use.empirical) {
        qq <- apeAce(bt, states, model)$rates
        gamma_prior$alpha <- qq * gamma_prior$beta
    }

    XX <- mcmcQ(bt, states, model, new_tree, tol, n_states, burn_in, sample_freq, n_simulations, vQ, gamma_prior, omp)
    # We compute the likelihood again anyway, maybe we could skip this

    join_Q_with_prior <- function(idx, Q_list, prior_list, logL_list) {
        return (list(Q=Q_list[[idx]]
                     , prior=prior_list[[idx]]
                     , logL=logL_list[[idx]])
                )
    }

    Q <- lapply(XX, function(x) x$Q)
    L <- lapply(XX, function(x) x$L)
    logL <- lapply(XX, function(x) x$loglik)

    if (prior[1] == "equal") {
        prior <- setNames(rep(1/n_states,n_states), colnames(L)) # set equal
        prior <- lapply(1:n_simulations, function(x,y) y, y=prior)
    } else if (prior[1] == "estimated") {
        prior <- lapply(Q, statdist) # set from stationary distribution
    } else {
        prior <- prior/sum(prior) # obtain from input
        prior <- lapply(1:n_simulations, function(x,y) y, y=prior)
    }

    return(lapply(1:n_simulations, join_Q_with_prior
                    , Q_list=Q
                    , prior_list=prior
                    , logL_list=logL))

}


# These are mostly helper functions honestly stolen and slightly changed from phytools.
# We had to copy them because they are not exported.

# PHYTOOLS
# get pars
# written by Liam J. Revell 2013
getPars<-function(bt,xx,model,Q,tree,tol,m,omp=1,liks=TRUE){
    XX<-apeAce(bt,xx,model,omp,fixedQ=Q,output.liks=liks)
    N<-length(bt$tip.label)
    II<-XX$index.matrix
    lvls<-XX$states
    if(liks){
        L<-XX$lik.anc
        rownames(L)<-N+1:nrow(L)
        if(!ape::is.binary(tree)){
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
.getSEs <- function(out)
{
    h <- out$hessian
    if (any(diag(h) == 0)) {
        warning("The likelihood gradient seems flat in at least one dimension (gradient null):\ncannot compute the standard-errors of the transition rates.\n")
        se <- rep(NaN, nrow(h))
    } else {
        se <- sqrt(diag(solve(h)))
    }
    se
}

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
    phy<-reorder.phylo(tree,"pruningwise")

    dev<-function(p,output.liks=FALSE,fixedQ=NULL){
        if(any(is.nan(p))||any(is.infinite(p))) return(1e+50)
        comp<-numeric(nb.tip+nb.node)
        if(is.null(fixedQ)){
            Q[]<-c(p,0)[rate]
            diag(Q)<--rowSums(Q)
        } else Q<-fixedQ

        Q_eigen <- eigen(Q, only.values = FALSE)
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
        out <- nlminb(rep(ip, length.out = np), function(p) dev(p),
                      lower = rep(0, np), upper = rep(1e50, np))
        obj<-list()
        obj$loglik <- -out$objective/2
        obj$rates <- out$par
        oldwarn <- options("warn")
        options(warn = -1)
        out.nlm <- try(nlm(function(p) dev(p), p = obj$rates, iterlim = 1,
                           stepmax = 0, hessian = TRUE), silent = TRUE)
        options(oldwarn)
        obj$se <- if (inherits(out.nlm, "try-error")) {
            warning("model fit suspicious: gradients apparently non-finite")
            rep(NaN, np)
        } else .getSEs(out.nlm)
        obj$index.matrix <- index.matrix
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
mcmcQ<-function(bt,xx,model,tree,tol,m,burnin,samplefreq,nsim,vQ,prior,omp){
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
    yy<-getPars(bt,xx,model,Q,tree,tol,m,omp)
    #cat("Running MCMC burn-in. Please wait....\n")
    for(i in 1:burnin){
        pp<-update(p)
        Qp<-matrix(pp[rate],m,m)
        diag(Qp)<--rowSums(Qp,na.rm=TRUE)
        zz<-getPars(bt,xx,model,Qp,tree,tol,m,omp,FALSE)
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
        zz<-getPars(bt,xx,model,Qp,tree,tol,m,omp,FALSE)
        p.odds<-exp(zz$loglik+sum(dgamma(pp,prior$alpha,prior$beta,log=TRUE))-
            yy$loglik-sum(dgamma(p,prior$alpha,prior$beta,log=TRUE)))
        if(p.odds>=runif(n=1)){
            yy<-zz
            p<-pp
        }
        if(i%%samplefreq==0){
            Qi<-matrix(p[rate],m,m)
            diag(Qi)<--rowSums(Qi,na.rm=TRUE)
            XX[[i/samplefreq]]<-getPars(bt,xx,model,Qi,tree,tol,m,omp,TRUE)
        }
    }
    return(XX)
}
