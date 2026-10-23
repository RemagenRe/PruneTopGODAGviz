mini_obo_path <- function() {
    system.file("extdata", "mini.obo", package = "PruneTopGODAGviz")
}

mini_dag <- function(ont = "BP") {
    read_obo(mini_obo_path(), ont = ont, verbose = FALSE)
}

mini_annotation <- function() {
    dag <- read_obo(mini_obo_path(), verbose = FALSE)
    tx <- sprintf("T%03d", 1:120)
    photo <- tx[1:30]
    stress <- tx[31:60]
    other <- tx[61:120]

    add <- function(df, ids, go) {
        rbind(df, data.frame(transcript = rep(ids, each = length(go)),
                             go_id = rep(go, times = length(ids)),
                             stringsAsFactors = FALSE))
    }
    pairs <- data.frame(transcript = character(), go_id = character(),
                        stringsAsFactors = FALSE)
    pairs <- add(pairs, photo, c("GO:0015979", "GO:0019684", "GO:0009765",
                                 "GO:0006091", "GO:0009507", "GO:0009579",
                                 "GO:0009535"))
    pairs <- add(pairs, stress, c("GO:0006950", "GO:0009409", "GO:0009651",
                                  "GO:0009415", "GO:0009414"))
    pairs <- add(pairs, other[1:20], c("GO:0009987", "GO:0008152", "GO:0005634"))
    pairs <- add(pairs, other[21:40], c("GO:0003674", "GO:0003824", "GO:0016787"))
    pairs <- add(pairs, other[41:60], c("GO:0005575", "GO:0005737"))
    pairs <- add(pairs, photo[1:10], c("GO:0003677", "GO:0005634"))
    pairs <- unique(pairs)

    terms <- dag$nodes[, c("id", "name", "namespace")]
    names(terms) <- c("go_id", "name", "ontology")
    terms <- terms[terms$go_id %in% pairs$go_id, , drop = FALSE]
    terms <- unique(terms)
    rownames(terms) <- NULL

    universe <- lapply(c("BP", "MF", "CC"), function(o) {
        g <- terms$go_id[terms$ontology == o]
        sort(unique(pairs$transcript[pairs$go_id %in% g]))
    })
    names(universe) <- c("BP", "MF", "CC")

    out <- list(pairs = pairs, terms = terms, universe = universe,
                keytype = "transcript",
                source = list(file = "synthetic", size = NA, mtime = NA,
                              n_rows = nrow(pairs), n_bad = 0L, sep = ","))
    class(out) <- "go_annotation"
    out
}

mini_genes <- function() {
    sprintf("T%03d", c(1:30, 31:45))
}

mini_pruned <- function(top_n = 4) {
    a <- mini_annotation()
    e <- enrich_go(mini_genes(), a, ont = "BP", verbose = FALSE)
    prune_go_dag(e, dag = mini_dag("BP"), top_n = top_n, verbose = FALSE)
}

expect_layout_clean <- function(lay) {
    n <- lay$nodes
    for (i in seq_len(nrow(n) - 1L)) {
        for (j in (i + 1L):nrow(n)) {
            dx <- abs(n$x[i] - n$x[j]) - (n$w[i] + n$w[j]) / 2
            dy <- abs(n$y[i] - n$y[j]) - (n$h[i] + n$h[j]) / 2
            testthat::expect_false(dx < 0 && dy < 0,
                                   info = paste("overlap", n$id[i], n$id[j]))
        }
    }
    testthat::expect_true(all(is.finite(n$x)) && all(is.finite(n$y)))
}
