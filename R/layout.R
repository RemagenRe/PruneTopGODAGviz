#' Lay out a pruned GO DAG
#'
#' Deterministic, collision-free layered layout.  Box sizes are computed from
#' the actual text metrics, so labels never overlap and never leave the canvas.
#' Long edges are routed through virtual nodes, which keeps them away from the
#' boxes and reduces crossings.
#'
#' @param dag A `pruned_go_dag` object.
#' @param rankdir `"TB"` (top to bottom) or `"LR"` (left to right).
#' @param level_gap Gap between levels, in inches.
#' @param sibling_gap Minimum gap between neighbouring boxes, in inches.
#' @param node_pad Inner padding of a box, in inches.
#' @param fontsize Text size in points.
#' @param lineheight Multiplier for the line height.
#' @param family Font family used for measuring the text.
#' @param edge_style `"hybrid"` (smooth short edges, orthogonal long edges),
#'   `"orthogonal"`, `"bezier"` or `"straight"`.
#' @param verbose Print a summary.
#'
#' @return A `go_dag_layout` object with `nodes`, `edges`, `segments`, `bbox`,
#'   `rankdir` and the layout parameters.
#' @export
layout_go_dag <- function(dag, rankdir = c("TB", "LR"), level_gap = 0.28,
                          sibling_gap = 0.12, node_pad = 0.06, fontsize = 7,
                          lineheight = 0.9, family = "Arial",
                          edge_style = c("hybrid", "orthogonal", "bezier", "straight"),
                          verbose = TRUE) {
    if (!inherits(dag, "pruned_go_dag")) {
        stop("'dag' must be a pruned_go_dag object.", call. = FALSE)
    }
    rankdir <- match.arg(rankdir)
    edge_style <- match.arg(edge_style)
    family <- .font_family(family)
    nodes <- dag$nodes
    edges <- dag$edges
    if (!nrow(nodes)) stop("No terms to lay out.", call. = FALSE)
    nodes <- .node_boxes(nodes, fontsize, lineheight, family, node_pad)

    lev <- .dag_levels(list(nodes = nodes, edges = edges))
    nodes$level <- as.integer(lev[nodes$id])
    levels <- sort(unique(nodes$level))

    aug <- .add_dummies(nodes, edges, lev, levels)
    an <- aug$nodes
    ae <- aug$edges

    order_lv <- lapply(levels, function(l) {
        ids <- an$id[an$level == l]
        ids[order(an$name[match(ids, an$id)], ids)]
    })
    r0 <- integer(nrow(an)); names(r0) <- an$id
    for (l in seq_along(order_lv)) r0[order_lv[[l]]] <- seq_along(order_lv[[l]])
    for (l in seq_along(levels)[-1]) {
        ids <- order_lv[[l]]
        b <- vapply(ids, function(v) {
            nb <- ae$parent[ae$child == v]
            nb <- nb[nb %in% names(r0)]
            if (!length(nb)) NA_real_ else mean(r0[nb])
        }, numeric(1))
        order_lv[[l]] <- ids[order(ifelse(is.na(b), r0[ids], b), r0[ids])]
        r0[order_lv[[l]]] <- seq_along(order_lv[[l]])
    }

    .rank <- function(ol) {
        r <- integer(nrow(an))
        names(r) <- an$id
        for (l in seq_along(ol)) r[ol[[l]]] <- seq_along(ol[[l]])
        r
    }
    .crossings <- function(ol) {
        r <- .rank(ol)
        total <- 0
        for (l in seq_along(levels)[-length(levels)]) {
            e <- ae[ae$parent %in% ol[[l]] & ae$child %in% ol[[l + 1]], ,
                    drop = FALSE]
            if (nrow(e) < 2) next
            a <- r[e$parent]
            b <- r[e$child]
            total <- total + sum(outer(a, a, "<") & outer(b, b, ">"))
        }
        total
    }
    .bary <- function(ids, r, use_parent = TRUE) {
        vapply(ids, function(v) {
            nb <- if (use_parent) {
                ae$parent[ae$child == v]
            } else {
                ae$child[ae$parent == v]
            }
            nb <- nb[nb %in% names(r)]
            if (!length(nb)) NA_real_ else mean(r[nb])
        }, numeric(1))
    }

    best <- order_lv
    best_n <- .crossings(order_lv)
    for (seed in list(order_lv, lapply(order_lv, rev))) {
        ol <- seed
        for (iter in seq_len(10)) {
            r <- .rank(ol)
            if (iter %% 2 == 1) {
                for (l in seq_along(ol)[-1]) {
                    ids <- ol[[l]]
                    b <- .bary(ids, r, TRUE)
                    ol[[l]] <- ids[order(ifelse(is.na(b), r[ids], b), r[ids])]
                }
            } else {
                for (l in rev(seq_along(ol))[-1]) {
                    ids <- ol[[l]]
                    b <- .bary(ids, r, FALSE)
                    ol[[l]] <- ids[order(ifelse(is.na(b), r[ids], b), r[ids])]
                }
            }
        }
        n <- .crossings(ol)
        if (n < best_n) {
            best_n <- n
            best <- ol
        }
    }
    order_lv <- best

    x <- stats::setNames(numeric(nrow(an)), an$id)
    for (l in seq_along(levels)) {
        ids <- order_lv[[l]]
        ww <- an$w[match(ids, an$id)]
        dd <- an$dummy[match(ids, an$id)]
        if (length(ids) == 1L) {
            x[ids] <- 0
        } else {
            gap <- ifelse(dd[-length(dd)] | dd[-1], 0.04, sibling_gap)
            x[ids] <- cumsum(c(0, (ww[-length(ww)] + ww[-1]) / 2 + gap))
        }
    }

    for (pass in seq_len(4)) {
        for (l in seq_along(levels)) {
            ids <- order_lv[[l]]
            target <- vapply(ids, function(v) {
                nb <- if (pass %% 2 == 1) {
                    ae$parent[ae$child == v]
                } else {
                    ae$child[ae$parent == v]
                }
                nb <- nb[nb %in% names(x)]
                if (!length(nb)) NA_real_ else mean(x[nb])
            }, numeric(1))
            target[is.na(target)] <- x[ids][is.na(target)]
            newx <- x[ids] + (target - x[ids])
            if (length(ids) > 1L) {
                ww <- an$w[match(ids, an$id)]
                dd <- an$dummy[match(ids, an$id)]
                gap <- ifelse(dd[-length(dd)] | dd[-1], 0.04, sibling_gap)
                for (i in 2:length(ids)) {
                    newx[i] <- max(newx[i],
                                   newx[i - 1] + (ww[i - 1] + ww[i]) / 2 + gap[i - 1])
                }
            }
            x[ids] <- newx + mean(target - newx)
        }
    }

    maxh <- vapply(seq_along(levels), function(l) {
        hh <- an$h[match(order_lv[[l]], an$id)]
        if (all(hh == 0)) 0 else max(hh)
    }, numeric(1))
    ytop <- cumsum(c(0, maxh[-length(maxh)] + level_gap))
    y <- stats::setNames(numeric(nrow(an)), an$id)
    for (l in seq_along(levels)) {
        ids <- order_lv[[l]]
        y[ids] <- ytop[l] + maxh[l] / 2
    }

    an$x <- unname(x[an$id])
    an$y <- unname(y[an$id])
    an$ord <- unname(.rank(order_lv)[an$id])

    real <- an[!an$dummy, , drop = FALSE]
    real <- real[match(nodes$id, real$id), , drop = FALSE]
    real$n_lines <- NULL
    nodes$x <- real$x
    nodes$y <- real$y
    nodes$ord <- real$ord

    if (rankdir == "LR") {
        colw <- vapply(seq_along(levels), function(l) {
            hh <- an$w[match(order_lv[[l]], an$id)]
            max(hh)
        }, numeric(1))
        xcol <- cumsum(c(0, (colw[-length(colw)] + colw[-1]) / 2 + level_gap)) +
            colw / 2
        ymax <- max(an$x - an$w / 2) + 0.1
        lvl0 <- an$level - min(levels) + 1L
        an$y <- ymax - an$x
        an$x <- xcol[lvl0]
        nodes$x <- an$x[match(nodes$id, an$id)]
        nodes$y <- an$y[match(nodes$id, an$id)]
    }

    if (rankdir == "LR") {
        flow_top <- stats::setNames(xcol - colw / 2, levels)
        flow_bottom <- stats::setNames(xcol + colw / 2, levels)
    } else {
        flow_top <- stats::setNames(ytop, levels)
        flow_bottom <- stats::setNames(ytop + maxh, levels)
    }

    cx <- (min(nodes$x - nodes$w / 2) + max(nodes$x + nodes$w / 2)) / 2
    nodes$x <- nodes$x - cx
    an$x <- an$x - cx
    if (rankdir == "LR") {
        flow_top <- flow_top - cx
        flow_bottom <- flow_bottom - cx
    }

    ed <- .edge_routes(edges, nodes, an, flow_top, flow_bottom, levels,
                       gap = 0.03, rankdir = rankdir)

    if (rankdir == "LR") {
        bbox <- list(xmin = min(nodes$x - nodes$w / 2),
                     xmax = max(nodes$x + nodes$w / 2),
                     ymin = min(nodes$y - nodes$h / 2),
                     ymax = max(nodes$y + nodes$h / 2))
    } else {
        bbox <- list(xmin = min(nodes$x - nodes$w / 2),
                     xmax = max(nodes$x + nodes$w / 2),
                     ymin = min(nodes$y - nodes$h / 2),
                     ymax = max(nodes$y + nodes$h / 2))
    }
    bbox$width <- bbox$xmax - bbox$xmin
    bbox$height <- bbox$ymax - bbox$ymin

    out <- list(nodes = nodes, edges = ed, segments = NULL,
                bbox = bbox, rankdir = rankdir, edge_style = edge_style,
                params = list(level_gap = level_gap, sibling_gap = sibling_gap,
                              node_pad = node_pad, fontsize = fontsize,
                              lineheight = lineheight, family = family,
                              crossings = best_n))
    class(out) <- "go_dag_layout"
    if (verbose) {
        message(sprintf("Layout: %d nodes, %d levels, crossings = %d, %d virtual nodes",
                        nrow(nodes), length(levels), best_n, sum(an$dummy)))
    }
    out
}

.add_dummies <- function(nodes, edges, lev, levels) {
    nodes$dummy <- FALSE
    if (!nrow(edges)) {
        return(list(nodes = nodes, edges = transform(
            edges, edge = seq_len(nrow(edges)), dummy = FALSE)))
    }
    dummy_id <- character()
    dummy_level <- integer()
    ae <- vector("list", nrow(edges))
    for (i in seq_len(nrow(edges))) {
        p <- edges$parent[i]
        ch <- edges$child[i]
        lp <- lev[[p]]
        lc <- lev[[ch]]
        if (is.na(lp) || is.na(lc) || lc - lp <= 1L) {
            ae[[i]] <- data.frame(parent = p, child = ch, type = edges$type[i],
                                  edge = i, dummy = FALSE,
                                  stringsAsFactors = FALSE)
            next
        }
        prev <- p
        rows <- list()
        for (l in seq.int(lp + 1L, lc - 1L)) {
            did <- sprintf(".d%d_%d", i, l)
            dummy_id <- c(dummy_id, did)
            dummy_level <- c(dummy_level, l)
            rows[[length(rows) + 1L]] <- data.frame(
                parent = prev, child = did, type = edges$type[i],
                edge = i, dummy = TRUE, stringsAsFactors = FALSE)
            prev <- did
        }
        rows[[length(rows) + 1L]] <- data.frame(
            parent = prev, child = ch, type = edges$type[i],
            edge = i, dummy = FALSE, stringsAsFactors = FALSE)
        ae[[i]] <- do.call(rbind, rows)
    }
    ae <- do.call(rbind, ae)
    if (!length(dummy_id)) {
        return(list(nodes = nodes, edges = ae))
    }
    dn <- data.frame(
        id = dummy_id, name = "", namespace = NA_character_,
        definition = NA_character_, obsolete = FALSE,
        label_l1 = "", label_l2 = "", label_l3 = "",
        selected = FALSE, n_collapsed = 0L, level = dummy_level,
        w = 0.02, h = 0, dummy = TRUE, stringsAsFactors = FALSE)
    common <- intersect(names(dn), names(nodes))
    an <- rbind(nodes[, common, drop = FALSE], dn[, common, drop = FALSE])
    an <- an[match(c(nodes$id, dn$id), an$id), , drop = FALSE]
    list(nodes = an, edges = ae)
}

.edge_routes <- function(edges, nodes, an, flow_top, flow_bottom, levels,
                         gap = 0.03, rankdir = "TB") {
    ed <- data.frame(
        parent = edges$parent, child = edges$child,
        type = edges$type %||% rep("is_a", nrow(edges)),
        n_skipped = edges$n_skipped %||% rep(0L, nrow(edges)),
        stringsAsFactors = FALSE)
    pi_ <- match(ed$parent, nodes$id)
    ci_ <- match(ed$child, nodes$id)
    routes <- vector("list", nrow(ed))
    routes_o <- vector("list", nrow(ed))
    for (i in seq_len(nrow(ed))) {
        p <- edges$parent[i]
        ch <- edges$child[i]
        lp <- an$level[match(p, an$id)]
        lc <- an$level[match(ch, an$id)]
        if (rankdir == "TB") {
            pts <- rbind(c(nodes$x[pi_[i]], nodes$y[pi_[i]] + nodes$h[pi_[i]] / 2))
        } else {
            pts <- rbind(c(nodes$x[pi_[i]] + nodes$w[pi_[i]] / 2, nodes$y[pi_[i]]))
        }
        mids <- numeric()
        for (l in seq.int(lp + 1L, lc - 1L)) {
            did <- sprintf(".d%d_%d", i, l)
            k <- match(did, an$id)
            if (is.na(k)) next
            pts <- rbind(pts, c(an$x[k], an$y[k]))
            mids <- c(mids, (flow_bottom[match(l - 1L, levels)] +
                                 flow_top[match(l, levels)]) / 2)
        }
        if (rankdir == "TB") {
            pts <- rbind(pts, c(nodes$x[ci_[i]], nodes$y[ci_[i]] - nodes$h[ci_[i]] / 2))
        } else {
            pts <- rbind(pts, c(nodes$x[ci_[i]] - nodes$w[ci_[i]] / 2, nodes$y[ci_[i]]))
        }
        routes[[i]] <- pts

        po <- pts[1, , drop = FALSE]
        for (k in 2:nrow(pts)) {
            lvl <- lp + k - 1L
            mid <- (flow_bottom[match(lvl - 1L, levels)] +
                        flow_top[match(lvl, levels)]) / 2
            if (rankdir == "TB") {
                po <- rbind(po, c(po[nrow(po), 1], mid),
                            c(pts[k, 1], mid), c(pts[k, 1], pts[k, 2]))
            } else {
                po <- rbind(po, c(mid, po[nrow(po), 2]),
                            c(mid, pts[k, 2]), c(pts[k, 1], pts[k, 2]))
            }
        }
        routes_o[[i]] <- po
    }
    ed$route <- routes
    ed$route_ortho <- routes_o
    ed
}

#' @export
print.go_dag_layout <- function(x, ...) {
    cat("<go_dag_layout>", x$rankdir, x$edge_style, "\n")
    cat(sprintf("  canvas: %.2f x %.2f in\n", x$bbox$width, x$bbox$height))
    cat("  nodes :", nrow(x$nodes), "\n")
    invisible(x)
}
#' Rescale a GO DAG layout uniformly
#'
#' Multiplies every geometric quantity of a layout (node positions and box sizes,
#' edge end points and edge routes, the bounding box and the gap parameters) by
#' `factor`, so a figure can be rendered on a canvas of a different physical size
#' without changing its proportions. Combine it with `fontsize * factor` in
#' [plot_go_dag()] for a true uniform scaling of the whole figure.
#'
#' @param layout A `go_dag_layout` object.
#' @param factor Positive numeric scale factor (1 leaves the layout untouched).
#' @return A rescaled `go_dag_layout` object.
#' @export
scale_go_layout <- function(layout, factor) {
    if (!inherits(layout, "go_dag_layout")) {
        stop("'layout' must be a go_dag_layout object.", call. = FALSE)
    }
    if (!is.numeric(factor) || length(factor) != 1L ||
        !is.finite(factor) || factor <= 0) {
        stop("'factor' must be a single positive finite number.", call. = FALSE)
    }
    if (abs(factor - 1) < 1e-9) {
        return(layout)
    }
    numeric_cols <- c("x", "y", "xend", "yend", "w", "h", "top", "bot", "yy")
    route_cols <- c("route", "route_ortho")
    scale_frame <- function(df) {
        if (is.null(df) || !length(df) || !nrow(df)) {
            return(df)
        }
        for (nm in intersect(numeric_cols, names(df))) {
            if (is.numeric(df[[nm]])) {
                df[[nm]] <- df[[nm]] * factor
            }
        }
        for (nm in intersect(route_cols, names(df))) {
            df[[nm]] <- lapply(df[[nm]], function(m) {
                if (is.matrix(m) && is.numeric(m)) m * factor else m
            })
        }
        df
    }
    layout$nodes <- scale_frame(layout$nodes)
    layout$edges <- scale_frame(layout$edges)
    if (!is.null(layout$segments)) {
        layout$segments <- scale_frame(layout$segments)
    }
    if (!is.null(layout$bbox)) {
        bbox <- layout$bbox
        for (nm in intersect(c("xmin", "xmax", "ymin", "ymax", "width", "height"),
                             names(bbox))) {
            if (is.numeric(bbox[[nm]])) {
                bbox[[nm]] <- bbox[[nm]] * factor
            }
        }
        layout$bbox <- bbox
    }
    if (!is.null(layout$params)) {
        for (nm in intersect(c("level_gap", "sibling_gap", "node_pad"),
                             names(layout$params))) {
            if (is.numeric(layout$params[[nm]])) {
                layout$params[[nm]] <- layout$params[[nm]] * factor
            }
        }
    }
    layout
}

#' Fit a GO DAG layout to a target page size
#'
#' @param dag A `pruned_go_dag` object.
#' @param width_in,height_in Target page size in inches.
#' @param fontsize Text size in points.
#' @param lineheight,family Passed to [layout_go_dag()].
#' @param level_gap,sibling_gap,node_pad Base spacing in inches.
#' @param gap_floor Lower bound for the gap multipliers.
#' @param edge_style Edge routing style passed to [layout_go_dag()].
#' @param rankdirs Layout directions to try.
#' @return A list with `layout`, `fontsize`, `margin`, `text_pt`, `shrunk`,
#'   `gaps`, `rankdir`, `fits` and `size_in`.
#' @keywords internal
fit_go_layout_to_page <- function(dag, width_in, height_in, fontsize,
                                  lineheight = 0.9, family = "Arial",
                                  level_gap = 0.28, sibling_gap = 0.12,
                                  node_pad = 0.06, gap_floor = 0.12,
                                  edge_style = "hybrid",
                                  rankdirs = c("TB", "LR")) {
    size_of <- function(lay) plot_go_dag_size(dag, layout = lay)
    build <- function(fs, bs, mL, mS, rd) {
        layout_go_dag(dag, rankdir = rd, fontsize = fs, lineheight = lineheight,
                      family = family, node_pad = node_pad * bs,
                      level_gap = level_gap * bs * mL,
                      sibling_gap = sibling_gap * bs * mS,
                      edge_style = edge_style, verbose = FALSE)
    }
    probes <- function(fs, bs, rd) {
        s11 <- size_of(build(fs, bs, 1, 1, rd))
        sl <- size_of(build(fs, bs, 2, 1, rd))
        ss <- size_of(build(fs, bs, 1, 2, rd))
        list(w0 = s11$width, h0 = s11$height,
             dWL = sl$width - s11$width, dHL = sl$height - s11$height,
             dWS = ss$width - s11$width, dHS = ss$height - s11$height)
    }
    solve_gaps <- function(pb) {
        det <- pb$dWL * pb$dHS - pb$dWS * pb$dHL
        if (!is.finite(det) || abs(det) < 1e-9) {
            return(c(1, 1))
        }
        r1 <- width_in - pb$w0
        r2 <- height_in - pb$h0
        c(1 + (r1 * pb$dHS - r2 * pb$dWS) / det,
          1 + (r2 * pb$dWL - r1 * pb$dHL) / det)
    }
    attempt_dir <- function(rd) {
        fs <- fontsize
        bs <- 1
        pb <- probes(fs, bs, rd)
        k <- min(if (pb$w0 > 0) width_in / pb$w0 else 1,
                 if (pb$h0 > 0) height_in / pb$h0 else 1)
        if (is.finite(k) && k < 1) {
            bs <- k * 0.97
            fs <- fontsize * bs
            pb <- probes(fs, bs, rd)
        }
        g <- solve_gaps(pb)
        mL <- max(gap_floor, g[1])
        mS <- max(gap_floor, g[2])
        lay <- build(fs, bs, mL, mS, rd)
        sz <- size_of(lay)
        fits <- sz$width <= width_in * 1.02 && sz$height <= height_in * 1.02
        if (!fits) {
            sc <- min(width_in / sz$width, height_in / sz$height)
            if (is.finite(sc) && sc > 0) {
                lay <- scale_go_layout(lay, sc)
                fs <- fs * sc
                sz <- size_of(lay)
                fits <- sz$width <= width_in * 1.02 && sz$height <= height_in * 1.02
            }
        }
        list(layout = lay, fontsize = fs, text_pt = fs, rankdir = rd,
             gaps = c(level = mL, sibling = mS), fits = isTRUE(fits),
             size_in = c(width = sz$width, height = sz$height))
    }
    best <- NULL
    for (rd in rankdirs) {
        cur <- attempt_dir(rd)
        if (is.null(best) || cur$text_pt > best$text_pt * 1.02 ||
            (cur$text_pt > best$text_pt * 0.98 && isTRUE(cur$fits) && !isTRUE(best$fits))) {
            best <- cur
        }
    }
    best$shrunk <- best$text_pt < fontsize - 0.05
    best$margin <- 0.08
    best
}