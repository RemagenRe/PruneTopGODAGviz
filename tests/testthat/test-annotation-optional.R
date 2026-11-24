test_that("minimal two-column annotation works when ont is supplied", {
    f <- tempfile(fileext = ".tsv")
    writeLines(c("transcript_id\tgo_id",
                 "T1\tGO:0008150", "T2\tGO:0008150"), f)
    a <- read_go_annotation(f, ont = "BP", keytype = "transcript",
                            verbose = FALSE)
    expect_s3_class(a, "go_annotation")
    expect_equal(nrow(a$pairs), 2L)
    expect_equal(a$terms$name, "GO:0008150")
    expect_equal(unique(a$terms$ontology), "BP")
})

test_that("minimal annotation without ontology is rejected clearly", {
    f <- tempfile(fileext = ".tsv")
    writeLines(c("transcript_id\tgo_id", "T1\tGO:0008150"), f)
    expect_error(read_go_annotation(f, keytype = "transcript", verbose = FALSE),
                 "No ontology information")
})

test_that("term name is optional when ontology is explicit", {
    f <- tempfile(fileext = ".tsv")
    writeLines(c("transcript_id\tgo_id\tterm_name",
                 "T1\tGO:0008150\tbiological_process"), f)
    a <- read_go_annotation(f, ont = "BP", keytype = "transcript",
                            verbose = FALSE)
    expect_equal(a$terms$name, "biological_process")
})

test_that("ont filters a multi-ontology file without noise", {
    f <- tempfile(fileext = ".tsv")
    writeLines(c("transcript_id\tgo_id\tnamespace",
                 "T1\tGO:0008150\tbiological_process",
                 "T2\tGO:0003674\tmolecular_function"), f)
    expect_silent(a <- read_go_annotation(f, ont = "BP",
                                          keytype = "transcript",
                                          verbose = FALSE))
    expect_equal(a$pairs$transcript, "T1")
    expect_equal(a$terms$go_id, "GO:0008150")
})

test_that("explicit invalid column names still fail", {
    f <- tempfile(fileext = ".tsv")
    writeLines(c("transcript_id\tgo_id", "T1\tGO:0008150"), f)
    expect_error(read_go_annotation(f, col_name = "not_a_column",
                                    ont = "BP", keytype = "transcript",
                                    verbose = FALSE),
                 "not found")
})
