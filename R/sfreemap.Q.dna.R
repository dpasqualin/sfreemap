Q_dna <- function(tip_states, tree, model, tol) {

    contrast = matrix(data = c(1,0,0,0,0,
                               0,1,0,0,0,
                               0,0,1,0,0,
                               0,0,0,1,0,
                               1,0,1,0,0,
                               0,1,0,1,0,
                               0,0,0,0,1,
                               1,1,1,1,0,
                               1,1,1,1,1)
                        , ncol = 5
                        , byrow = TRUE)
    dimnames(contrast) = list(
        c("a","c","g","t","r","y","-","n","?"),
        c("a", "c", "g", "t", "-")
    )

    # FIXME:
    # Workaround because pml does not uses inherits, so objects with multiple
    # inheritance causes a warning message. I've sent a pull request for them,
    # when it's accepted and published on cran we can remove this
    class(tree) <- "phylo"

    m <- getModels()[[model]]
    data <- phyDat(tip_states, type="USER", contrast=contrast)
    fit <- pml(tree, data)
    ctrl <- pml.control(trace=0)
    fit <- suppressWarnings(
        optim.pml(fit, optQ=m$optQ, optBf=m$optBf, subs=m$subs, control=ctrl)
    )

    # Reconstruct Q matrix
    levels <- attr(fit$data, "levels")
    nc <- attr(fit$data, "nc")
    QM <- matrix(0, nc, nc, dimnames = list(levels, levels))
    QM[lower.tri(QM)] <- fit$Q
    QM <- QM+t(QM)

    # Value smaller than tol are set to tol
    QM[QM<tol] <- tol

    # set diagonal
    diag(QM) <- 0
    diag(QM) <- rowSums(QM) * -1

    return (list(Q=QM, prior=fit$bf, logL=fit$logLik))
}

# from package phangorn
# helps to get parameters for different DNA models
getModels <- function() {
    models <- list(
         JC = list(optQ=FALSE, optBf=FALSE,   subs=c(0, 0, 0, 0, 0, 0)),
         F81 = list(optQ=FALSE, optBf=TRUE,   subs=c(0, 0, 0, 0, 0, 0)),
         K80 = list(optQ=TRUE, optBf=FALSE,   subs=c(0, 1, 0, 0, 1, 0)),
         HKY = list(optQ=TRUE, optBf=TRUE,    subs=c(0, 1, 0, 0, 1, 0)),
         TrNe = list(optQ=TRUE, optBf=FALSE,  subs=c(0, 1, 0, 0, 2, 0)),
         TrN = list(optQ=TRUE, optBf=TRUE,    subs=c(0, 1, 0, 0, 2, 0)),
         TPM1 = list(optQ=TRUE, optBf=FALSE,  subs=c(0, 1, 2, 2, 1, 0)),
         K81 = list(optQ=TRUE, optBf=FALSE,   subs=c(0, 1, 2, 2, 1, 0)),
         TPM1u = list(optQ=TRUE, optBf=TRUE,  subs=c(0, 1, 2, 2, 1, 0)),
         TPM2 = list(optQ=TRUE, optBf=FALSE,  subs=c(1, 2, 1, 0, 2, 0)),
         TPM2u = list(optQ=TRUE, optBf=TRUE,  subs=c(1, 2, 1, 0, 2, 0)),
         TPM3 = list(optQ=TRUE, optBf=FALSE,  subs=c(1, 2, 0, 1, 2, 0)),
         TPM3u = list(optQ=TRUE, optBf=TRUE,  subs=c(1, 2, 0, 1, 2, 0)),
         TIM1e = list(optQ=TRUE, optBf=FALSE, subs=c(0, 1, 2, 2, 3, 0)),
         TIM1 = list(optQ=TRUE, optBf=TRUE,   subs=c(0, 1, 2, 2, 3, 0)),
         TIM2e = list(optQ=TRUE, optBf=FALSE, subs=c(1, 2, 1, 0, 3, 0)),
         TIM2 = list(optQ=TRUE, optBf=TRUE,   subs=c(1, 2, 1, 0, 3, 0)),
         TIM3e = list(optQ=TRUE, optBf=FALSE, subs=c(1, 2, 0, 1, 3, 0)),
         TIM3 = list(optQ=TRUE, optBf=TRUE,   subs=c(1, 2, 0, 1, 3, 0)),
         TVMe = list(optQ=TRUE, optBf=FALSE,  subs=c(1, 2, 3, 4, 2, 0)),
         TVM = list(optQ=TRUE, optBf=TRUE,    subs=c(1, 2, 3, 4, 2, 0)),
         SYM = list(optQ=TRUE, optBf=FALSE,   subs=c(1, 2, 3, 4, 5, 0)),
         GTR = list(optQ=TRUE, optBf=TRUE,    subs=c(1, 2, 3, 4, 5, 0))
    )
    return (models)
}
