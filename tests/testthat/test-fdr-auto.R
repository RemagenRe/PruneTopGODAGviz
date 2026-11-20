test_that("auto FDR breaks are ordered, bounded and keep contrast", {
    p <- c(1e-6, 1e-4, 0.01, 0.05, 0.004)
    br <- PruneTopGODAGviz:::.auto_fdr_breaks(p)
    expect_length(br, 3)
    expect_true(br[[1]] > br[[2]] && br[[2]] > br[[3]])
    expect_lte(br[["blue"]], max(p))
    expect_gte(br[["red"]], min(p))
    expect_gte(br[["red"]], 1e-300)
    expect_lte(br[["blue"]], 1)

    pos <- PruneTopGODAGviz:::.fdr_position(p, br[[1]], br[[2]], br[[3]])
    expect_gte(length(unique(round(pos, 2))), 4)
    expect_gt(diff(range(pos)), 0.6)
})

test_that("auto FDR breaks separate terms that all sit far below 1e-8", {
    # Typical of GO tables: every displayed term is extremely significant.
    p <- c(2.6e-88, 1.0e-75, 1.0e-69, 3.0e-60, 8.5e-57)
    br <- PruneTopGODAGviz:::.auto_fdr_breaks(p)
    expect_true(br[[1]] > br[[2]] && br[[2]] > br[[3]])
    pos <- PruneTopGODAGviz:::.fdr_position(p, br[[1]], br[[2]], br[[3]])
    expect_equal(length(unique(round(pos, 2))), length(p))
    expect_gt(diff(range(pos)), 0.9)
})

test_that("auto FDR breaks widen a narrow range and cope with one value", {
    narrow <- PruneTopGODAGviz:::.auto_fdr_breaks(c(1e-4, 2e-4, 1e-3))
    expect_true(narrow[[1]] > narrow[[2]] && narrow[[2]] > narrow[[3]])
    expect_gt(narrow[[1]] / narrow[[3]], 9)

    one <- PruneTopGODAGviz:::.auto_fdr_breaks(rep(1e-7, 5))
    expect_true(one[[1]] > one[[2]] && one[[2]] > one[[3]])
})

test_that("auto FDR breaks stay valid for extreme or missing values", {
    wide <- PruneTopGODAGviz:::.auto_fdr_breaks(c(1e-300, 1e-12, 0.02))
    expect_true(wide[[1]] > wide[[2]] && wide[[2]] > wide[[3]])
    expect_lte(wide[["blue"]], 1)
    expect_gte(wide[["red"]], 1e-300)

    filtered <- PruneTopGODAGviz:::.auto_fdr_breaks(
        c(NA, Inf, 0, -1, 0.02, 1e-5))
    expect_equal(unname(filtered[["blue"]]), 0.02)

    # Values above 1 are invalid but must not break the ordering.
    over <- PruneTopGODAGviz:::.auto_fdr_breaks(c(1.5, 2, 0.5))
    expect_true(over[[1]] > over[[2]] && over[[2]] > over[[3]])
    expect_lte(over[["blue"]], 1)

    expect_error(PruneTopGODAGviz:::.auto_fdr_breaks(c(NA, NaN)),
                 "No FDR values")
})

test_that("auto FDR breaks are recomputed after the terms change", {
    p <- c(1e-5, 1e-4, 1e-3, 2e-2, 0.2)
    before <- PruneTopGODAGviz:::.auto_fdr_breaks(p)
    after <- PruneTopGODAGviz:::.auto_fdr_breaks(p[-1])
    expect_false(isTRUE(all.equal(unname(before), unname(after))))
    expect_gt(after[["red"]], before[["red"]])

    after_top <- PruneTopGODAGviz:::.auto_fdr_breaks(p[-length(p)])
    expect_lt(after_top[["blue"]], before[["blue"]])
})

test_that("auto FDR breaks accept the p.adjust of a pruned graph", {
    p <- mini_pruned()$nodes$p.adjust
    br <- PruneTopGODAGviz:::.auto_fdr_breaks(p)
    expect_true(br[[1]] > br[[2]] && br[[2]] > br[[3]])
})

test_that("colour fills differ across terms after an auto fit", {
    p <- c(2.6e-88, 1.0e-75, 1.0e-69, 3.0e-60, 8.5e-57, 1.2e-40)
    br <- PruneTopGODAGviz:::.auto_fdr_breaks(p)
    nodes <- data.frame(
        id = paste0("GO:", seq_along(p)),
        selected = TRUE,
        p.adjust = p,
        stringsAsFactors = FALSE)
    colors <- list(low = "#2878b5", mid = "#ffffff", high = "#c82423",
                   ancestor = "#9E9E9E")
    fills <- PruneTopGODAGviz:::.vis_colours(nodes, colors, br)
    expect_equal(length(fills), nrow(nodes))
    expect_gte(length(unique(fills)), 4)
})

test_that("break labels switch to scientific notation only when needed", {
    expect_equal(PruneTopGODAGviz:::.fmt_break(c(0.05, 0.01, 0.0001)),
                 c("0.05", "0.01", "0.0001"))
    expect_equal(PruneTopGODAGviz:::.fmt_break(3.2e-5), "3.2e-05")
    expect_equal(PruneTopGODAGviz:::.fmt_break(1e-30), "1.0e-30")
})
test_that("typed FDR breaks accept scientific notation and fall back", {
    expect_true(PruneTopGODAGviz:::.fdr_breaks_ok(
        c("9.7e-68", "7.6e-74", "1.2e-111")))
    expect_true(PruneTopGODAGviz:::.fdr_breaks_ok(c("0.05", "0.01", "0.0001")))
    expect_equal(
        PruneTopGODAGviz:::.fdr_breaks_from_text(c("0.05", "0.01", "0.0001")),
        c(0.05, 0.01, 0.0001))
    expect_equal(
        PruneTopGODAGviz:::.fdr_breaks_from_text(c("9e-68", "1e-74", "1e-111")),
        c(9e-68, 1e-74, 1e-111))

    expect_false(PruneTopGODAGviz:::.fdr_breaks_ok(c("0.01", "0.05", "0.0001")))
    expect_false(PruneTopGODAGviz:::.fdr_breaks_ok(c("abc", "0.01", "0.0001")))
    expect_false(PruneTopGODAGviz:::.fdr_breaks_ok(c("0.05", "0.01")))
    expect_equal(
        PruneTopGODAGviz:::.fdr_breaks_from_text(c("abc", "0.01", "0.0001")),
        c(0.05, 0.01, 0.0001))
})