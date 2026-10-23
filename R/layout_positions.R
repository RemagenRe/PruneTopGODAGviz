#' Lay out a GO DAG from interactive coordinates
#'
#' Builds a publication layout from node coordinates captured in the Shiny
#' editor (for example with `visNetwork::visGetPositions()`).  Box sizes are
#' recomputed from the publication text metrics and overlapping boxes are
#' separated, so dragged layouts can be exported without collisions.
#'
#' @param dag A `pruned_go_dag` object.
#' @param positions Node coordinates: a data frame with `id`, `x`, `y`, a named
#'   list of `c(x = , y = )`, or a two-column matrix with node ids as row names.
#' @param fontsize Text size in points.
#' @param lineheight Line height multiplier.
#' @param family Font family used for measuring text.
#' @param node_pad Inner padding of a box, in inches.
#' @param scale Divisor converting interactive units to inches (96 px per inch).
#' @param relax Run a collision relaxation pass.
#' @param relax_iter Maximum number of relaxation iterations.
#' @param gap Minimum gap between boxes, in inches.
#' @param rankdir Layout direction used for nodes without coordinates.
#' @param verbose Print a summary.
#'
#' @return A `go_dag_layout` object.
#' @export
layout_from_positions <- function(dag, positions, fontsize = 7,
                                  lineheight = 0.9, family = "Arial",
                                  node_pad = 0.06, scale = 1 / 96,
                                  relax = TRUE, relax_iter = 80, gap = 0.04,
                                  rankdir = "TB", verbose = TRUE) {
    if (!inherits(dag, "pruned_go_dag")) {
        stop("'dag' must be a pruned_go_dag object.", call. = FALSE)
    }
    family <- .font_family(family)
    nodes <- .node_boxes(dag$nodes, fontsize, lineheight, family, node_pad)
    pos <- .positions_df(positions)
    pos <- pos[!duplicated(pos$id), , drop = FALSE]

    x <- stats::setNames(rep(NA_real_, nrow(nodes)), nodes$id)
    y <- stats::setNames(rep(NA_real_, nrow(nodes)), nodes$id)
    hit <- match(pos$id, nodes$id)
    ok <- !is.na(hit) & is.finite(pos$x) & is.finite(pos$y)
    x[hit[ok]] <- pos$x[ok] * scale
    y[hit[ok]] <- pos$y[ok] * scale

    if (anyNA(x) || anyNA(y)) {
        base <- layout_go_dag(dag, rankdir = rankdir, fontsize = fontsize,
                              lineheight = lineheight, family = family,
                              node_pad = node_pad, verbose = FALSE)
        bx <- stats::setNames(base$nodes$x, base$nodes$id)
        by <- stats::setNames(base$nodes$y, base$nodes$id)
        miss <- is.na(x) | is.na(y)
        x[miss] <- bx[names(x)[miss]]
        y[miss] <- by[names(y)[miss]]
    }

    x <- x - min(x - nodes$w / 2)
    y <- y - min(y - nodes$h / 2)

    if (isTRUE(relax)) {
        iter <- min(relax_iter,
                    max(20L, ceiling(2e5 / max(1, nrow(nodes) * (nrow(nodes) - 1) / 2))))
        r <- .relax_boxes(unname(x), unname(y), nodes$w, nodes$h,
                          gap = gap, iter = iter)
        x[] <- r$x
        y[] <- r$y
    }

    nodes$x <- unname(x[nodes$id])
    nodes$y <- unname(y[nodes$id])

    ed <- data.frame(
        parent = dag$edges$parent, child = dag$edges$child,
        type = dag$edges$type %||% rep("is_a", nrow(dag$edges)),
        n_skipped = dag$edges$n_skipped %||% rep(0L, nrow(dag$edges)),
        stringsAsFactors = FALSE)
    pi_ <- match(ed$parent, nodes$id)
    ci_ <- match(ed$child, nodes$id)
    routes <- vector("list", nrow(ed))
    routes_o <- vector("list", nrow(ed))
    for (i in seq_len(nrow(ed))) {
        if (is.na(pi_[i]) || is.na(ci_[i])) {
            routes[[i]] <- NULL
            routes_o[[i]] <- NULL
            next
        }
        if (rankdir == "TB") {
            start <- c(nodes$x[pi_[i]], nodes$y[pi_[i]] + nodes$h[pi_[i]] / 2)
            end <- c(nodes$x[ci_[i]], nodes$y[ci_[i]] - nodes$h[ci_[i]] / 2)
        } else {
            start <- c(nodes$x[pi_[i]] + nodes$w[pi_[i]] / 2, nodes$y[pi_[i]])
            end <- c(nodes$x[ci_[i]] - nodes$w[ci_[i]] / 2, nodes$y[ci_[i]])
        }
        routes[[i]] <- rbind(start, end)
        if (rankdir == "TB" && abs(end[1] - start[1]) > 1e-6) {
            mid <- (start[2] + end[2]) / 2
            routes_o[[i]] <- rbind(start, c(start[1], mid), c(end[1], mid), end)
        } else if (rankdir == "LR" && abs(end[2] - start[2]) > 1e-6) {
            mid <- (start[1] + end[1]) / 2
            routes_o[[i]] <- rbind(start, c(mid, start[2]), c(mid, end[2]), end)
        } else {
            routes_o[[i]] <- routes[[i]]
        }
    }
    ed$route <- routes
    ed$route_ortho <- routes_o

    bbox <- list(xmin = min(nodes$x - nodes$w / 2),
                 xmax = max(nodes$x + nodes$w / 2),
                 ymin = min(nodes$y - nodes$h / 2),
                 ymax = max(nodes$y + nodes$h / 2))
    bbox$width <- bbox$xmax - bbox$xmin
    bbox$height <- bbox$ymax - bbox$ymin

    out <- list(nodes = nodes, edges = ed, segments = NULL, bbox = bbox,
                rankdir = rankdir, edge_style = "hybrid", source = "positions",
                params = list(level_gap = NA_real_, sibling_gap = gap,
                              node_pad = node_pad, fontsize = fontsize,
                              lineheight = lineheight, family = family,
                              crossings = NA_integer_))
    class(out) <- "go_dag_layout"
    if (verbose) {
        message(sprintf("Layout from positions: %d nodes, %.2f x %.2f in",
                        nrow(nodes), bbox$width, bbox$height))
    }
    out
}

.positions_df <- function(positions) {
    if (is.null(positions) || !length(positions)) {
        return(data.frame(id = character(), x = numeric(), y = numeric(),
                          stringsAsFactors = FALSE))
    }
    if (is.data.frame(positions)) {
        if (!all(c("id", "x", "y") %in% names(positions))) {
            stop("'positions' must have columns id, x and y.", call. = FALSE)
        }
        return(data.frame(id = as.character(positions$id),
                          x = as.numeric(positions$x),
                          y = as.numeric(positions$y), stringsAsFactors = FALSE))
    }
    if (is.matrix(positions)) {
        return(data.frame(id = rownames(positions),
                          x = as.numeric(positions[, 1]),
                          y = as.numeric(positions[, 2]), stringsAsFactors = FALSE))
    }
    if (is.list(positions)) {
        ids <- names(positions)
        xy <- t(vapply(positions, function(p) {
            nm <- names(p)
            if (!is.null(nm) && all(c("x", "y") %in% nm)) {
                px <- p[["x"]]; py <- p[["y"]]
            } else {
                px <- p[[1]]; py <- p[[2]]
            }
            c(as.numeric(px), as.numeric(py))
        }, numeric(2)))
        return(data.frame(id = ids, x = xy[, 1], y = xy[, 2],
                          stringsAsFactors = FALSE))
    }
    stop("Unsupported 'positions' format.", call. = FALSE)
}

.relax_boxes <- function(x, y, w, h, gap = 0.04, iter = 80) {
    n <- length(x)
    if (n < 2L) return(list(x = x, y = y))
    for (it in seq_len(iter)) {
        moved <- FALSE
        for (i in seq_len(n - 1L)) {
            for (j in (i + 1L):n) {
                dx <- x[j] - x[i]
                dy <- y[j] - y[i]
                ox <- (w[i] + w[j]) / 2 + gap - abs(dx)
                oy <- (h[i] + h[j]) / 2 + gap - abs(dy)
                if (ox > 0 && oy > 0) {
                    if (ox < oy) {
                        s <- if (dx >= 0) 1 else -1
                        x[i] <- x[i] - s * ox / 2
                        x[j] <- x[j] + s * ox / 2
                    } else {
                        s <- if (dy >= 0) 1 else -1
                        y[i] <- y[i] - s * oy / 2
                        y[j] <- y[j] + s * oy / 2
                    }
                    moved <- TRUE
                }
            }
        }
        if (!moved) break
    }
    list(x = x, y = y)
}
