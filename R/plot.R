#' Draw a pruned GO DAG
#'
#' Renders the graph with ggplot2.  Boxes carry three lines (identifier and
#' statistic, term name, gene counts), significant terms are coloured on a
#' continuous scale, and edges are styled by relation type.  The canvas size is
#' derived from the layout unless explicitly supplied.
#'
#' @param x A `pruned_go_dag` object.
#' @param layout Optional `go_dag_layout`; computed when `NULL`.
#' @param fill Three colours used for the significance gradient.
#' @param stat Statistic used for the colour gradient.
#' @param label_stat Statistic printed in the first line.
#' @param show_name,show_counts,show_p Toggle the three label lines.
#' @param wrap If numeric, re-wrap term names at this width.
#' @param fontsize Text size in points.
#' @param font_family Font family.
#' @param edge_types Style edges by `is_a` / `part_of`.
#' @param legend Show legends.
#' @param title Optional plot title.
#' @param rankdir `"TB"` or `"LR"`.
#' @param lineheight Line height multiplier.
#' @param margin Canvas margin in inches.
#' @param ancestor_fill Fill colour for non-significant ancestor terms.
#' @param bg Canvas background colour.
#' @param edge_col Named colours for edge types (is_a, part_of,
#'   unknown).
#' @param fdr_breaks Three FDR values for the blue, white and red ends of the
#'   colour scale, in decreasing order. Defaults to `c(0.05, 0.01, 0.0001)`.
#' @param ... Passed to [layout_go_dag()].
#'
#' @return A ggplot object.
#' @export
plot_go_dag <- function(x, layout = NULL, fill = c("#2878B5", "#FFFFFF", "#c82423"),
                        stat = NULL, label_stat = NULL, show_name = TRUE,
                        show_counts = TRUE, show_p = TRUE, wrap = NULL,
                        fontsize = 7, font_family = "Arial",
                        edge_types = TRUE, legend = TRUE, title = NULL,
                        rankdir = c("TB", "LR"), lineheight = 0.9,
                        margin = 0.08, ancestor_fill = "#FFFFFF", bg = "white",
                        edge_col = c(is_a = "#000000", part_of = "#2878B5",
                                     unknown = "#000000"),
                        fdr_breaks = c(0.05, 0.01, 0.0001), ...) {
    if (!inherits(x, "pruned_go_dag")) {
        stop("'x' must be a pruned_go_dag object.", call. = FALSE)
    }
    rankdir <- match.arg(rankdir)
    font_family <- .font_family(font_family)
    fdr_breaks <- .check_fdr_breaks(fdr_breaks)
    stat <- stat %||% x$stat
    label_stat <- label_stat %||% stat
    for (s in c(stat, label_stat)) {
        if (!s %in% c("pvalue", "p.adjust")) {
            stop("Unknown statistic: ", s, call. = FALSE)
        }
    }
    layout <- layout %||% layout_go_dag(x, rankdir = rankdir, fontsize = fontsize,
                                        lineheight = lineheight,
                                        family = font_family, verbose = FALSE)

    nodes <- layout$nodes
    node_pad <- layout$params$node_pad
    if (is.numeric(wrap)) {
        nodes$label_l2 <- .wrap_label(nodes$name, width = wrap, max_lines = Inf)
        nodes$n_lines <- 1L + .n_lines(nodes$label_l2) + 1L
    }

    nodes$lab1 <- nodes$id
    if (isTRUE(show_p)) {
        v <- nodes[[label_stat]]
        has <- nodes$selected & !is.na(v)
        nodes$lab1[has] <- paste0(nodes$id[has], " (", .fmt_p(v[has]), ")")
    }
    nodes$lab2 <- if (isTRUE(show_name)) nodes$label_l2 else rep("", nrow(nodes))
    nodes$lab3 <- if (isTRUE(show_counts)) nodes$label_l3 else rep("", nrow(nodes))

    fv <- .fdr_position(nodes$p.adjust, fdr_breaks[1], fdr_breaks[2],
                          fdr_breaks[3])
    fv[!nodes$selected] <- NA_real_
    nodes$fill_value <- fv

    pad <- node_pad
    nodes$yy <- -nodes$y
    nodes$top <- nodes$yy + nodes$h / 2
    nodes$bot <- nodes$yy - nodes$h / 2
    nodes$label_all <- paste(nodes$lab1, nodes$lab2, nodes$lab3, sep = "\n")

    lim <- c(0, 1)

    p <- ggplot2::ggplot() +
        ggplot2::geom_rect(
            data = nodes,
            ggplot2::aes(xmin = .data$x - .data$w / 2,
                         xmax = .data$x + .data$w / 2,
                         ymin = .data$bot, ymax = .data$top,
                         fill = .data$fill_value),
            colour = "black", linewidth = 0.25, show.legend = legend) +
        ggplot2::geom_rect(
            data = nodes,
            ggplot2::aes(xmin = .data$x - .data$w / 2 + 0.02,
                         xmax = .data$x + .data$w / 2 - 0.02,
                         ymin = .data$bot + 0.02,
                         ymax = .data$top - 0.02),
            fill = "white", alpha = 0.4, colour = NA) +
        ggplot2::geom_text(
            data = nodes,
            ggplot2::aes(x = .data$x, y = .data$yy, label = .data$label_all),
            size = fontsize / .pt_to_mm, family = font_family, colour = "black",
            lineheight = lineheight, vjust = 0.5) +
        ggplot2::scale_fill_gradientn(
            colours = fill, values = c(0, 0.5, 1), limits = lim,
            breaks = c(0, 0.5, 1),
            labels = .fmt_break(fdr_breaks),
            na.value = ancestor_fill,
            name = "FDR", guide = "none")

    n_node_layers <- length(p$layers)
    paths <- .edge_paths(layout)
    n_edge_layers <- 0L
    if (nrow(paths)) {
        if (edge_types) {
            p <- p + ggplot2::geom_path(
                data = paths,
                ggplot2::aes(x = .data$x, y = .data$y, group = .data$edge,
                             colour = .data$type, linetype = .data$type),
                linewidth = 0.3, show.legend = FALSE)
            p <- p + ggplot2::scale_colour_manual(
                values = edge_col, guide = "none", drop = TRUE) +
                ggplot2::scale_linetype_manual(
                    values = c(is_a = "solid", part_of = "dashed",
                               unknown = "solid"),
                    guide = "none", drop = TRUE)
            heads <- .arrow_heads(paths)
            for (tp in unique(heads$type)) {
                p <- p + ggplot2::geom_polygon(
                    data = heads[heads$type == tp, , drop = FALSE],
                    ggplot2::aes(x = .data$x, y = .data$y, group = .data$edge),
                    fill = unname(edge_col[tp]), colour = NA,
                    show.legend = FALSE)
            }
        } else {
            p <- p + ggplot2::geom_path(
                data = paths,
                ggplot2::aes(x = .data$x, y = .data$y, group = .data$edge),
                colour = "black", linewidth = 0.3)
            heads <- .arrow_heads(paths)
            p <- p + ggplot2::geom_polygon(
                data = heads,
                ggplot2::aes(x = .data$x, y = .data$y, group = .data$edge),
                fill = "black", colour = NA, show.legend = FALSE)
        }
        n_edge_layers <- length(p$layers) - n_node_layers
    }
    xr <- c(layout$bbox$xmin - margin, layout$bbox$xmax + margin)
    yr <- c(-layout$bbox$ymax - margin, -layout$bbox$ymin + margin)
    legend_extra <- 0
    if (isTRUE(legend)) {
        u <- fontsize / 72
        lw <- max(1.90, u * 12)
        lh <- max(1.15, u * 8)
        pos <- .free_legend_rect(layout, lw, lh)
        if (!isTRUE(pos$inside)) {
            yr[2] <- yr[2] + lh + 0.16
            legend_extra <- lh + 0.16
        }
        for (ly in .inset_legend(pos$x, pos$y, lw, lh, fill, edge_col,
                                 fontsize, font_family, bg, fdr_breaks)) {
            p <- p + ly
        }
    }
    p <- p +
        ggplot2::coord_fixed(xlim = xr, ylim = yr, expand = FALSE) +
        ggplot2::labs(title = title) +
        ggplot2::theme_void(base_family = font_family) +
        ggplot2::theme(
            legend.position = "none",
            legend.title = ggplot2::element_text(size = 9, face = "plain"),
            legend.text = ggplot2::element_text(size = 8.5, face = "plain"),
            legend.key.size = grid::unit(0.5, "cm"),
            legend.key.width = grid::unit(0.8, "cm"),
            legend.margin = ggplot2::margin(1, 1, 1, 1),
            plot.title = ggplot2::element_text(size = 8, face = "plain",
                                               hjust = 0),
            plot.margin = ggplot2::margin(3, 3, 3, 3),
            panel.background = ggplot2::element_rect(fill = bg,
                                                     colour = NA),
            plot.background = ggplot2::element_rect(fill = bg,
                                                    colour = NA))
    if (n_edge_layers > 0L) {
        idx_edge <- seq.int(n_node_layers + 1L, n_node_layers + n_edge_layers)
        idx_rest <- if (n_node_layers + n_edge_layers < length(p$layers)) {
            seq.int(n_node_layers + n_edge_layers + 1L, length(p$layers))
        } else integer()
        p$layers <- c(p$layers[idx_edge], p$layers[seq_len(n_node_layers)],
                      p$layers[idx_rest])
    }
    attr(p, "go_layout") <- layout
    attr(p, "go_size") <- list(
        width = layout$bbox$width + 2 * margin,
        height = layout$bbox$height + 2 * margin +
            (if (is.null(title)) 0 else 0.3) + legend_extra)
    p
}

.free_legend_rect <- function(layout, width, height, pad = 0.06) {
    b <- layout$bbox
    xmin <- b$xmin
    xmax <- b$xmax
    ymin <- -b$ymax
    ymax <- -b$ymin
    n <- layout$nodes
    nx0 <- n$x - n$w / 2
    nx1 <- n$x + n$w / 2
    ny0 <- -n$y - n$h / 2
    ny1 <- -n$y + n$h / 2
    free <- function(x0, y0) {
        !any(nx1 > x0 - pad & nx0 < x0 + width + pad &
                 ny1 > y0 - pad & ny0 < y0 + height + pad)
    }
    cand <- list(
        c(xmin + pad, ymax - pad - height),
        c(xmax - width - pad, ymax - pad - height),
        c(xmin + pad, ymin + pad),
        c(xmax - width - pad, ymin + pad))
    xs <- seq(xmin + pad, xmax - width - pad, length.out = 15)
    ys <- seq(ymin + pad, ymax - height - pad, length.out = 15)
    for (y in ys) {
        for (x in xs) cand[[length(cand) + 1L]] <- c(x, y)
    }
    for (p in cand) {
        if (free(p[1], p[2])) {
            return(list(x = p[1], y = p[2], inside = TRUE))
        }
    }
    list(x = xmin + pad, y = ymax + 0.10, inside = FALSE)
}

.inset_legend <- function(x0, y0, width, height, fill, edge_col, fontsize,
                          family, bg = "white",
                          fdr_breaks = c(0.05, 0.01, 0.0001)) {
    col_a <- unname(edge_col[["is_a"]])
    col_p <- unname(edge_col[["part_of"]])
    if (is.na(col_a)) col_a <- "#000000"
    if (is.na(col_p)) col_p <- "#2878B5"
    layers <- list(
        ggplot2::annotate("rect", xmin = x0, xmax = x0 + width,
            ymin = y0, ymax = y0 + height, fill = bg, colour = NA))
    bar_x0 <- x0 + 0.10
    bar_x1 <- x0 + 0.18
    u <- fontsize / 72
    bar_y0 <- y0 + 0.14
    bar_y1 <- y0 + height - u * 1.6
    cols <- grDevices::colorRampPalette(fill)(81)
    ys <- seq(bar_y0, bar_y1, length.out = 82)
    for (i in seq_len(80)) {
        layers[[length(layers) + 1L]] <- ggplot2::annotate(
            "rect", xmin = bar_x0, xmax = bar_x1,
            ymin = ys[i], ymax = ys[i + 1L], fill = cols[i], colour = NA)
    }
    ts <- fontsize / .pt_to_mm
    layers[[length(layers) + 1L]] <- ggplot2::annotate(
        "text", x = bar_x0, y = y0 + height - 0.02, label = "FDR",
        hjust = 0, vjust = 1, size = ts, family = family)
    fdr_labels <- .fmt_break(fdr_breaks)
    for (k in list(c(bar_y0, fdr_labels[1]),
                   c((bar_y0 + bar_y1) / 2, fdr_labels[2]),
                   c(bar_y1, fdr_labels[3]))) {
        layers[[length(layers) + 1L]] <- ggplot2::annotate(
            "text", x = bar_x1 + 0.05, y = as.numeric(k[1]), label = k[2],
            hjust = 0, vjust = 0.5, size = ts, family = family)
    }
    rx <- x0 + 0.85
    rx1 <- x0 + 1.15
    tx <- x0 + 1.20
    y_a <- y0 + height - u * 3.0
    y_p <- y0 + height - u * 4.5
    layers[[length(layers) + 1L]] <- ggplot2::annotate(
        "text", x = rx, y = y0 + height - 0.02, label = "Relation",
        hjust = 0, vjust = 1, size = ts, family = family)
    layers[[length(layers) + 1L]] <- ggplot2::annotate(
        "segment", x = rx, xend = rx1, y = y_a, yend = y_a,
        colour = col_a, linewidth = 0.4)
    layers[[length(layers) + 1L]] <- ggplot2::annotate(
        "text", x = tx, y = y_a, label = "is_a", hjust = 0, vjust = 0.5,
        size = ts, family = family)
    layers[[length(layers) + 1L]] <- ggplot2::annotate(
        "segment", x = rx, xend = rx1, y = y_p, yend = y_p,
        colour = col_p, linewidth = 0.4, linetype = "dashed")
    layers[[length(layers) + 1L]] <- ggplot2::annotate(
        "text", x = tx, y = y_p, label = "part_of", hjust = 0, vjust = 0.5,
        size = ts, family = family)
    layers
}

.edge_paths <- function(layout, n = 40L) {
    ed <- layout$edges
    if (!nrow(ed)) {
        return(data.frame(edge = integer(), type = character(),
                          x = numeric(), y = numeric()))
    }
    style <- layout$edge_style
    for (i in seq_len(nrow(ed))) {
        pts <- ed$route[[i]]
        if (is.null(pts) || nrow(pts) < 2L) next
        if (style == "straight") {
            bx <- pts[, 1]
            by <- pts[, 2]
        } else if (style == "orthogonal" ||
                   (style == "hybrid" && nrow(pts) > 2L)) {
            po <- ed$route_ortho[[i]]
            if (is.null(po) || nrow(po) < 2L) {
                bx <- pts[, 1]; by <- pts[, 2]
            } else {
                bx <- po[, 1]; by <- po[, 2]
            }
        } else {
            sm <- .smooth_route(pts, n_per = max(8L, n %/% max(1L, nrow(pts) - 1L)),
                                rankdir = layout$rankdir)
            bx <- sm[, 1]; by <- sm[, 2]
        }
        if (length(bx) >= 2L) {
            last <- length(bx)
            dx <- bx[last - 1] - bx[last]
            dy <- by[last - 1] - by[last]
            dd <- sqrt(dx^2 + dy^2)
            if (dd > 1e-9) {
                bx[last] <- bx[last] + 0.03 * dx / dd
                by[last] <- by[last] + 0.03 * dy / dd
            }
        }
        ed$route[[i]] <- cbind(bx, by)
    }
    out <- vector("list", nrow(ed))
    for (i in seq_len(nrow(ed))) {
        pts <- ed$route[[i]]
        out[[i]] <- data.frame(edge = i, type = ed$type[i],
                               x = pts[, 1], y = -pts[, 2],
                               stringsAsFactors = FALSE)
    }
    do.call(rbind, out)
}

.arrow_heads <- function(paths, length = 0.09, width = 0.055) {
    out <- list()
    for (e in unique(paths$edge)) {
        p <- paths[paths$edge == e, , drop = FALSE]
        if (nrow(p) < 2L) next
        n <- nrow(p)
        tip <- c(p$x[n], p$y[n])
        prev <- c(p$x[n - 1L], p$y[n - 1L])
        d <- tip - prev
        len <- sqrt(sum(d^2))
        if (!is.finite(len) || len < 1e-9) next
        u <- d / len
        perp <- c(-u[2], u[1])
        base <- tip - u * length
        out[[length(out) + 1L]] <- data.frame(
            edge = e, type = p$type[1],
            x = c(tip[1], base[1] + perp[1] * width / 2,
                  base[1] - perp[1] * width / 2),
            y = c(tip[2], base[2] + perp[2] * width / 2,
                  base[2] - perp[2] * width / 2),
            stringsAsFactors = FALSE)
    }
    if (!length(out)) {
        return(data.frame(edge = integer(), type = character(),
                          x = numeric(), y = numeric()))
    }
    do.call(rbind, out)
}

.smooth_route <- function(pts, n_per = 15L, rankdir = "TB") {
    xs <- numeric(); ys <- numeric()
    for (i in seq_len(nrow(pts) - 1L)) {
        p <- pts[i, ]; q <- pts[i + 1L, ]
        if (rankdir == "TB") {
            c1 <- c(p[1], p[2] + (q[2] - p[2]) * 0.5)
            c2 <- c(q[1], q[2] - (q[2] - p[2]) * 0.5)
        } else {
            c1 <- c(p[1] + (q[1] - p[1]) * 0.5, p[2])
            c2 <- c(q[1] - (q[1] - p[1]) * 0.5, q[2])
        }
        t <- seq(0, 1, length.out = n_per)
        if (i > 1L) t <- t[-1]
        bx <- (1 - t)^3 * p[1] + 3 * (1 - t)^2 * t * c1[1] +
            3 * (1 - t) * t^2 * c2[1] + t^3 * q[1]
        by <- (1 - t)^3 * p[2] + 3 * (1 - t)^2 * t * c1[2] +
            3 * (1 - t) * t^2 * c2[2] + t^3 * q[2]
        xs <- c(xs, bx); ys <- c(ys, by)
    }
    cbind(xs, ys)
}
#' Canvas size needed by a GO DAG figure
#'
#' @param x A `pruned_go_dag` object.
#' @param layout Optional layout.
#' @param width,height Override the computed size (inches).
#' @param margin Margin in inches.
#' @param title Logical; reserve space for a title.
#' @param ... Passed to [layout_go_dag()].
#'
#' @return A list with `width` and `height` in inches.
#' @export
plot_go_dag_size <- function(x, layout = NULL, width = NULL, height = NULL,
                             margin = 0.08, title = FALSE, ...) {
    layout <- layout %||% layout_go_dag(x, verbose = FALSE, ...)
    w <- layout$bbox$width + 2 * margin
    h <- layout$bbox$height + 2 * margin +
        if (isTRUE(title)) 0.3 else 0
    list(width = width %||% w, height = height %||% h)
}

#' Save a GO DAG as PDF and/or PNG
#'
#' @param x A `pruned_go_dag` object.
#' @param prefix Output path without extension.
#' @param plot Optional ggplot object; created when `NULL`.
#' @param width,height Canvas size in inches; computed when `NULL`.
#' @param dpi Raster resolution (default 800).
#' @param formats Any of `"pdf"`, `"png"`, `"svg"` and `"tiff"`.
#' @param overwrite Overwrite existing files.
#' @param bg Canvas background colour.
#' @param ... Passed to [plot_go_dag()].
#'
#' @return The written file paths, invisibly.
#' @export
save_go_dag <- function(x, prefix, plot = NULL, width = NULL, height = NULL,
                        dpi = 800, formats = c("pdf", "png"),
                        overwrite = TRUE, bg = "white", ...) {
    p <- plot %||% plot_go_dag(x, ...)
    if (is.null(width) || is.null(height)) {
        sz <- attr(p, "go_size")
        if (is.null(sz)) {
            sz <- plot_go_dag_size(x, width = width, height = height,
                                   title = !is.null(p$labels$title))
        }
        width <- width %||% sz$width
        height <- height %||% sz$height
    }
    px_w <- width * dpi
    px_h <- height * dpi
    if (px_w > 20000 || px_h > 20000 || px_w * px_h > 3.0e8) {
        warning(sprintf("Large raster (%d x %d px). Consider rankdir = 'LR', stronger pruning, or a lower dpi.",
                        as.integer(px_w), as.integer(px_h)), call. = FALSE)
    }
    out <- character()
    if ("pdf" %in% formats) {
        f <- paste0(prefix, ".pdf")
        if (file.exists(f) && !overwrite) stop("File exists: ", f, call. = FALSE)
        if (isTRUE(capabilities("cairo"))) {
            grDevices::cairo_pdf(f, width = width, height = height, bg = bg)
        } else {
            grDevices::pdf(f, width = width, height = height, bg = bg)
        }
        print(p)
        grDevices::dev.off()
        out <- c(out, f)
    }
    if ("png" %in% formats) {
        .check_pkg("ragg", "PNG export")
        f <- paste0(prefix, ".png")
        if (file.exists(f) && !overwrite) stop("File exists: ", f, call. = FALSE)
        ragg::agg_png(f, width = width, height = height, units = "in",
                      res = dpi, background = bg)
        print(p)
        grDevices::dev.off()
        out <- c(out, f)
    }
    if ("svg" %in% formats) {
        .check_pkg("svglite", "SVG export")
        f <- paste0(prefix, ".svg")
        if (file.exists(f) && !overwrite) stop("File exists: ", f, call. = FALSE)
        svglite::svglite(f, width = width, height = height, bg = bg)
        print(p)
        grDevices::dev.off()
        out <- c(out, f)
    }
    if ("tiff" %in% formats) {
        .check_pkg("ragg", "TIFF export")
        f <- paste0(prefix, ".tiff")
        if (file.exists(f) && !overwrite) stop("File exists: ", f, call. = FALSE)
        ragg::agg_tiff(f, width = width, height = height, units = "in",
                       res = dpi, background = bg, compression = "lzw")
        print(p)
        grDevices::dev.off()
        out <- c(out, f)
    }
    invisible(out)
}
