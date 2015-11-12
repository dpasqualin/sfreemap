context('sfreemap.map')

test_that ('works for phylo, standard, empirical', {

    sm <- sfreemap.map(sfreemap.corals.trees[[1]], sfreemap.corals.tips, type="standard", method="empirical")

    expect_true(inherits(sm, "phylo"))

    expect_false(any(is.nan(sm$mapped.edge)))
    expect_false(any(is.nan(sm$mapped.edge.lmt)))

})

test_that ('works for phylo, standard, empirical, Q given', {

    sm1 <- sfreemap.map(sfreemap.corals.trees[[1]], sfreemap.corals.tips, type="standard", method="empirical")
    sm2 <- sfreemap.map(sfreemap.corals.trees[[1]], sfreemap.corals.tips, Q=sm1$Q, type="standard", method="empirical")

    expect_true(inherits(sm1, "phylo"))
    expect_true(inherits(sm2, "phylo"))
    expect_equal(sm1, sm2)
})

test_that ('works for phylo, standard, mcmc', {

    sm <- sfreemap.map(sfreemap.corals.trees[[1]], sfreemap.corals.tips, type="standard", method="mcmc", n_simulation=10)

    expect_true(inherits(sm, "multiPhylo"))
    expect_equal(length(sm), 10)

    expect_false(any(sapply(sm, function(x) any(is.nan(x$mapped.edge)))))
    expect_false(any(sapply(sm, function(x) any(is.nan(x$mapped.edge.lmt)))))
})

test_that ('works for phylo, dna, single taxa', {

    sm <- sfreemap.map(sfreemap.primates.trees[[1]], sfreemap.primates.dna.tips[,1], type="dna")

    expect_true(inherits(sm, "phylo"))

    expect_false(any(is.nan(sm$mapped.edge)))
    expect_false(any(is.nan(sm$mapped.edge.lmt)))

})

test_that ('works for phylo, dna, multiple taxa', {

    sm <- sfreemap.map(sfreemap.primates.trees[[1]], sfreemap.primates.dna.tips[,1:10], type="dna")

    expect_true(inherits(sm, "multiPhylo"))
    expect_equal(length(sm), 10)

    expect_false(any(sapply(sm, function(x) any(is.nan(x$mapped.edge)))))
    expect_false(any(sapply(sm, function(x) any(is.nan(x$mapped.edge.lmt)))))

})

test_that ('works for phylo, dna, multiple taxa, Q given', {

    sm1 <- sfreemap.map(sfreemap.primates.trees[[1]], sfreemap.primates.dna.tips[,1:10], type="dna")
    Q <- lapply(sm1, function(x) x$Q)
    sm2 <- sfreemap.map(sfreemap.primates.trees[[1]], sfreemap.primates.dna.tips[,1:10], Q=Q, type="dna")

    expect_true(inherits(sm1, "multiPhylo"))
    expect_true(inherits(sm2, "multiPhylo"))

    expect_equal(sm1, sm2)
})
