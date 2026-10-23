test_that("mini.obo is parsed with the expected structure", {
    d <- read_obo(mini_obo_path(), verbose = FALSE)
    expect_s3_class(d, "go_dag")
    expect_equal(nrow(d$nodes), 23)
    expect_equal(sum(d$edges$type == "part_of"), 2)
    expect_true(all(c("GO:0009765", "GO:0009535") %in% d$edges$child[d$edges$type == "part_of"]))
    expect_false("GO:0000003" %in% d$nodes$id)
    expect_false("GO:9900001" %in% d$nodes$id)
    expect_true("GO:0009415" %in% d$nodes$id)
    expect_setequal(d$roots, c("GO:0008150", "GO:0003674", "GO:0005575"))
})

test_that("ontology filtering works", {
    bp <- read_obo(mini_obo_path(), ont = "BP", verbose = FALSE)
    expect_true(all(bp$nodes$namespace == "BP"))
    expect_equal(length(bp$roots), 1)
    expect_false(any(bp$edges$child %in% c("GO:0005634", "GO:0003674")))
})

test_that("dag helpers return complete ancestor and descendant sets", {
    d <- mini_dag("BP")
    anc <- dag_ancestors(d, "GO:0009414")
    expect_true(all(c("GO:0009414", "GO:0009415", "GO:0008150") %in% anc))
    desc <- dag_descendants(d, "GO:0015979")
    expect_true(all(c("GO:0019684", "GO:0009765") %in% desc))
    expect_false("GO:0009409" %in% desc)
})

test_that("a local OBO can be supplied through the package option", {
    old <- options(PruneTopGODAGviz.obo = mini_obo_path())
    on.exit(options(old), add = TRUE)
    d <- go_dag("BP", verbose = FALSE)
    expect_s3_class(d, "go_dag")
    expect_gt(nrow(d$nodes), 10)
    expect_true(all(d$nodes$namespace == "BP"))
})

test_that("GO.db is available as a fallback", {
    skip_if_not_installed("GO.db")
    d <- go_dag("BP", source = "GO.db", verbose = FALSE)
    expect_s3_class(d, "go_dag")
    expect_gt(nrow(d$nodes), 1000)
    expect_true(all(d$edges$type == "is_a"))
})
