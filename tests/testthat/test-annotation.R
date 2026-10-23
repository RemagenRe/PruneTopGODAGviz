test_that("the packaged example annotation is read correctly", {
    f <- system.file("extdata", "example_go_annotation.csv",
                     package = "PruneTopGODAGviz")
    expect_true(file.exists(f))
    a <- read_go_annotation(f, ont = "BP", verbose = FALSE)
    expect_s3_class(a, "go_annotation")
    expect_gt(nrow(a$pairs), 100)
    expect_true(all(grepl("^GO:[0-9]{7}$", a$terms$go_id)))
    expect_true(all(a$terms$ontology == "BP"))
    expect_true(all(a$pairs$go_id %in% a$terms$go_id))
    expect_identical(sort(unique(a$pairs$transcript)), a$universe$BP)
})

test_that("a leading empty column and quoted commas are handled", {
    tmp <- tempfile(fileext = ".csv")
    lines <- c(
        ',ensembl_gene_id,ensembl_transcript_id,go_id,name_1006,namespace_1003',
        '0,AT1G01010,AT1G01010.1,GO:0005634,"nucleus, membrane-bounded",cellular_component',
        '1,AT1G01010,AT1G01010.1,GO:0006355,regulation of DNA-templated transcription,biological_process')
    writeLines(lines, tmp)
    a <- read_go_annotation(tmp, verbose = FALSE)
    expect_equal(nrow(a$pairs), 2)
    expect_equal(a$terms$name[a$terms$go_id == "GO:0005634"],
                 "nucleus, membrane-bounded")
    expect_equal(a$keytype, "transcript")
})

test_that("gene identifiers can be used as keytype", {
    tmp <- tempfile(fileext = ".csv")
    lines <- c(
        "ensembl_gene_id,ensembl_transcript_id,go_id,name_1006,namespace_1003",
        "AT1G01010,AT1G01010.2,GO:0005634,nucleus,cellular_component",
        "AT1G01020,AT1G01020.1,GO:0005634,nucleus,cellular_component")
    writeLines(lines, tmp)
    a <- read_go_annotation(tmp, keytype = "gene", verbose = FALSE)
    expect_equal(a$keytype, "gene")
    expect_setequal(a$pairs$transcript, c("AT1G01010", "AT1G01020"))
})

test_that("invalid input is reported", {
    expect_error(read_go_annotation(tempfile(), verbose = FALSE),
                 "not found")
    tmp <- tempfile(fileext = ".csv")
    writeLines(c("a,b,c", "1,2,3"), tmp)
    expect_error(read_go_annotation(tmp, verbose = FALSE),
                 "Could not detect")
})

test_that("ontology filtering keeps only the requested ontology", {
    a <- mini_annotation()
    expect_true(all(a$terms$ontology %in% c("BP", "MF", "CC")))
    expect_equal(length(a$universe$BP), 80)
    expect_equal(length(a$universe$MF), 30)
    expect_equal(length(a$universe$CC), 70)
})
