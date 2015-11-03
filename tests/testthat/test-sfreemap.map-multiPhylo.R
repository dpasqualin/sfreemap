context('sfreemap.map')

test_that ('works for multiphylo, standard, empirical', {

    sm <- sfreemap.map(sfreemap.trees[1:10], sfreemap.tips, type="standard", method="empirical")

    expect_true(inherits(sm, "multiPhylo"))
    expect_equal(length(sm), 10)

    expect_false(any(sapply(sm, function(x) any(is.nan(x$mapped.edge)))))
    expect_false(any(sapply(sm, function(x) any(is.nan(x$mapped.edge.lmt)))))

})

test_that ('works for multiphylo, standard, empirical, Q given', {

    sm1 <- sfreemap.map(sfreemap.trees[1:10], sfreemap.tips, type="standard", method="empirical")
    Q <- lapply(sm1, function(x) x$Q)
    sm2 <- sfreemap.map(sfreemap.trees[1:10], sfreemap.tips, Q=Q, type="standard", method="empirical")

    expect_true(inherits(sm1, "multiPhylo"))
    expect_true(inherits(sm2, "multiPhylo"))
    expect_equal(sm1, sm2)
})

test_that ('works for multiphylo, standard, mcmc', {

    sm <- sfreemap.map(sfreemap.trees[1:5], sfreemap.tips, type="standard", method="mcmc", n_simulations=10)

    expect_true(inherits(sm, "multiPhylo"))
    expect_equal(length(sm), 50)

    expect_false(any(sapply(sm, function(x) any(is.nan(x$mapped.edge)))))
    expect_false(any(sapply(sm, function(x) any(is.nan(x$mapped.edge.lmt)))))

})

test_that ('works for multiphylo, standard, mcmc, Q given', {

    sm1 <- sfreemap.map(sfreemap.trees[1:5], sfreemap.tips, type="standard", method="mcmc", n_simulations=10)
    Q <- lapply(sm1, function(x) x$Q)
    sm2 <- sfreemap.map(sfreemap.trees[1:5], sfreemap.tips, Q=Q, type="standard", method="mcmc", n_simulations=10)

    expect_true(inherits(sm1, "multiPhylo"))
    expect_true(inherits(sm2, "multiPhylo"))
    expect_equal(sm1, sm2)
})
