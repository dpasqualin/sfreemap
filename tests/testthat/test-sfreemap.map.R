context('sfreemap.map')

trees_file <- system.file("extdata", "trees.data", package="sfreemapc")
tips_file <- system.file("extdata", "tips.data", package="sfreemapc")

trees <- read.tree(trees_file)
tips <- sfreemap.read_tips(tips_file, 1)

test_that ('return correct value for a single tree', {

    tree <- trees[[1]]
    sm <- sfreemap.map(tree, tips, Q='empirical')

    expect_false(is.null(sm$mapped.edge))
    expect_false(is.null(sm$mapped.edge.lmt))

    # we use signif because all.equal might get things wrong otherwise
    lmt_result <- signif(colSums(sm$mapped.edge.lmt), 6)
    emr_result <- signif(colSums(sm$mapped.edge), 6)

    # the real values
    lmt_expected <- c(1.133170, 0.135389)
    emr_expected <- c(3.895470, 0.229746)
    names(emr_expected) <- names(lmt_expected) <- c('a', 'b')

    expect_equal(lmt_result, lmt_expected)
    expect_equal(emr_result, emr_expected)
})

test_that ('we can describe the sfreemap.map result for a single tree', {

    tree <- trees[[1]]
    sm <- sfreemap.map(tree, tips, Q='empirical')

    desc <- sfreemap.describe(sm)

    transitions_result <- signif(desc$transitions, 6)
    dwelling_result <- signif(desc$dwelling_times, 6)

    transitions_expected <- c(1.133170, 0.135389)
    dwelling_expected <- c(3.895470, 0.229746)
    names(transitions_expected) <- names(dwelling_expected) <- c('a', 'b')

    expect_equal(transitions_result, transitions_expected)
    expect_equal(dwelling_result, dwelling_expected)

})
