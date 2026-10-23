#' Interactive HTML version of a GO DAG
#'
#' Produces a self-contained HTML widget in which terms can be hovered to show
#' the term name, statistics, counts and the input genes annotated to the term.
#'
#' @param x A `pruned_go_dag` object.
#' @param file Optional output HTML file.
#' @param width,height Widget size in pixels.
#' @param tooltip Any of `"term"`, `"p"`, `"counts"`, `"genes"`.
#' @param max_genes Maximum number of genes listed per tooltip.
#' @param layout Optional layout.
#' @param ... Passed to [plot_go_dag()] for styling.
#'
#' @return A `girafe` HTML widget.
#' @export
plot_go_dag_interactive <- function(x, file = NULL, width = 1000, height = 900,
                                    tooltip = c("term", "p", "counts", "genes"),
                                    max_genes = 20L, layout = NULL, ...) {
    .check_pkg("ggiraph", "interactive HTML output")
    .check_pkg("htmlwidgets", "interactive HTML output")
    if (!inherits(x, "pruned_go_dag")) {
        stop("'x' must be a pruned_go_dag object.", call. = FALSE)
    }
    layout <- layout %||% layout_go_dag(x, verbose = FALSE)
    dots <- list(...)
    fontsize <- dots$fontsize %||% 7
    font_family <- .font_family(dots$font_family %||% "Arial")
    fdr_breaks <- .check_fdr_breaks(dots$fdr_breaks %||%
                                        c(0.05, 0.01, 0.0001))
    lineheight <- dots$lineheight %||% 1.15
    node_pad <- layout$params$node_pad
    line_h <- fontsize * lineheight / 72

    nodes <- layout$nodes
    sets <- split(x$annotation$pairs$transcript, x$annotation$pairs$go_id)
    genes <- lapply(nodes$id, function(id) {
        intersect(sets[[id]] %||% character(), x$gene)
    })
    esc <- function(z) .html_escape(as.character(z %||% ""))

    tip <- vapply(seq_len(nrow(nodes)), function(i) {
        parts <- character()
        if ("term" %in% tooltip) {
            parts <- c(parts, esc(nodes$id[i]), esc(nodes$name[i]))
        }
        if ("p" %in% tooltip && !is.na(nodes$p.adjust[i])) {
            parts <- c(parts, paste0("FDR: ", esc(.fmt_p(nodes$p.adjust[i])),
                                     " | P: ", esc(.fmt_p(nodes$pvalue[i]))))
        }
        if ("counts" %in% tooltip) {
            parts <- c(parts, esc(nodes$label_l3[i]))
        }
        if ("genes" %in% tooltip) {
            g <- genes[[i]]
            if (length(g)) {
                shown <- utils::head(g, max_genes)
                more <- if (length(g) > length(shown)) {
                    paste0(" (+", length(g) - length(shown), " more)")
                } else ""
                parts <- c(parts, paste0("Genes: ", esc(paste(shown, collapse = ", ")),
                                         more))
            } else {
                parts <- c(parts, "Genes: none in input")
            }
        }
        paste(parts, collapse = "<br/>")
    }, character(1))

    nodes$tooltip <- tip
    nodes$yy <- -nodes$y
    nodes$top <- nodes$yy + nodes$h / 2
    nodes$bot <- nodes$yy - nodes$h / 2
    nodes$lab1 <- nodes$id
    has <- nodes$selected & !is.na(nodes$p.adjust)
    nodes$lab1[has] <- paste0(nodes$id[has], " (", .fmt_p(nodes$p.adjust[has]), ")")
    nodes$y1 <- nodes$top - node_pad - line_h / 2
    nodes$y3 <- nodes$bot + node_pad + line_h / 2
    nodes$y2 <- (nodes$y1 + nodes$y3) / 2
    fv <- .fdr_position(nodes$p.adjust, fdr_breaks[1], fdr_breaks[2],
                          fdr_breaks[3])
    fv[!nodes$selected] <- NA_real_
    nodes$fill_value <- fv
    lim <- c(0, 1)

    paths <- .edge_paths(layout)
    p <- ggplot2::ggplot() +
        ggiraph::geom_rect_interactive(
            data = nodes,
            ggplot2::aes(xmin = .data$x - .data$w / 2,
                         xmax = .data$x + .data$w / 2,
                         ymin = .data$bot, ymax = .data$top,
                         fill = .data$fill_value,
                         tooltip = .data$tooltip, data_id = .data$id),
            colour = "black", linewidth = 0.25) +
        ggplot2::geom_text(
            data = nodes,
            ggplot2::aes(x = .data$x, y = .data$y1, label = .data$lab1),
            size = fontsize / .pt_to_mm, family = font_family, colour = "black") +
        ggplot2::geom_text(
            data = nodes,
            ggplot2::aes(x = .data$x, y = .data$y2, label = .data$label_l2),
            size = fontsize / .pt_to_mm, family = font_family, colour = "black",
            lineheight = lineheight) +
        ggplot2::geom_text(
            data = nodes,
            ggplot2::aes(x = .data$x, y = .data$y3, label = .data$label_l3),
            size = fontsize / .pt_to_mm, family = font_family, colour = "black") +
        ggplot2::geom_path(
            data = paths,
            ggplot2::aes(x = .data$x, y = .data$y, group = .data$edge,
                         colour = .data$type, linetype = .data$type),
            linewidth = 0.3,
            arrow = ggplot2::arrow(type = "closed",
                                   length = grid::unit(1.6, "mm"))) +
        ggplot2::scale_fill_gradientn(
            colours = c("#2878B5", "#FFFFFF", "#c82423"), values = c(0, 0.5, 1),
            limits = lim, breaks = c(0, 0.5, 1),
            labels = .fmt_break(fdr_breaks),
            na.value = "#FFFFFF", name = "FDR") +
        ggplot2::scale_colour_manual(
            values = c(is_a = "black", part_of = "#2878b5", unknown = "black"),
            name = "Relation") +
        ggplot2::scale_linetype_manual(
            values = c(is_a = "solid", part_of = "dashed", unknown = "solid"),
            name = "Relation") +
        ggplot2::coord_fixed(
            xlim = c(layout$bbox$xmin - 0.25, layout$bbox$xmax + 0.25),
            ylim = c(-layout$bbox$ymax - 0.25, -layout$bbox$ymin + 0.25),
            expand = FALSE) +
        ggplot2::theme_void(base_family = font_family) +
        ggplot2::theme(legend.position = "right",
                       plot.background = ggplot2::element_rect(fill = "white",
                                                               colour = NA))

    w <- ggiraph::girafe(
        ggobj = p, width_svg = layout$bbox$width + 1.4,
        height_svg = layout$bbox$height + 0.6,
        options = list(
            ggiraph::opts_tooltip(
                css = "background-color:white;color:black;padding:4px;border:1px solid black;"),
            ggiraph::opts_hover(css = "fill:#F6CAE5;"),
            ggiraph::opts_zoom(max = 5),
            ggiraph::opts_sizing(rescale = TRUE)))
    if (!is.null(file)) {
        selfcontained <- nzchar(Sys.which("pandoc"))
        if (!selfcontained) {
            warning("pandoc not found; writing a non-self-contained HTML file.",
                    call. = FALSE)
        }
        htmlwidgets::saveWidget(w, file, selfcontained = selfcontained)
    }
    w
}
