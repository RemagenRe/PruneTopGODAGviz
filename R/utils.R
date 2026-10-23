`%||%` <- function(x, y) if (is.null(x)) y else x

.pt_to_mm <- 72.27 / 25.4

.check_ont <- function(ont) {
    ont <- toupper(trimws(ont))
    if (length(ont) != 1L || !ont %in% c("BP", "MF", "CC")) {
        stop("'ont' must be one of 'BP', 'MF' or 'CC'.", call. = FALSE)
    }
    ont
}

.norm_ontology <- function(x) {
    key <- tolower(trimws(as.character(x)))
    key <- gsub("[^a-z]", "", key)
    out <- rep(NA_character_, length(key))
    out[key %in% c("bp", "biologicalprocess", "p")] <- "BP"
    out[key %in% c("mf", "molecularfunction", "f")] <- "MF"
    out[key %in% c("cc", "cellularcomponent", "c")] <- "CC"
    out
}

.fmt_p <- function(p) {
    out <- ifelse(is.na(p), NA_character_,
                  format.pval(p, digits = 3, eps = 1e-300, scientific = TRUE))
    trimws(out)
}

.check_fdr_breaks <- function(x) {
    x <- as.numeric(x)
    if (length(x) != 3L || anyNA(x) || any(!is.finite(x)) ||
        any(x <= 0) || !(x[1] > x[2] && x[2] > x[3])) {
        stop("FDR colour breaks must be three positive numbers in decreasing order.",
             call. = FALSE)
    }
    x
}

.fmt_break <- function(x) {
    format(as.numeric(x), scientific = FALSE, trim = TRUE, drop0trailing = TRUE)
}

.fdr_position <- function(p, at_blue = 0.05, at_white = 0.01,
                           at_red = 0.0001) {
    p <- as.numeric(p)
    t <- rep(NA_real_, length(p))
    top <- log10(at_blue)
    mid <- log10(at_white)
    bot <- log10(at_red)
    high <- !is.na(p) & p >= at_white
    t[high] <- 0.5 * (top - log10(p[high])) / (top - mid)
    low <- !is.na(p) & p < at_white
    t[low] <- 0.5 + 0.5 * (mid - log10(pmax(p[low], 1e-300))) / (mid - bot)
    pmin(pmax(t, 0), 1)
}

.wrap_label <- function(x, width = 30, max_lines = 2L) {
    vapply(x, function(s) {
        if (is.na(s) || !nzchar(s)) return("")
        w <- strwrap(s, width = width)
        if (length(w) > max_lines) {
            w <- w[seq_len(max_lines)]
            last <- w[length(w)]
            w[length(w)] <- paste0(substr(last, 1, max(1L, nchar(last) - 1L)),
                                   "\u2026")
        }
        paste(w, collapse = "\n")
    }, character(1), USE.NAMES = FALSE)
}

.text_width <- function(x, size = 7, family = "Arial", weight = "normal") {
    x <- as.character(x)
    x[is.na(x)] <- ""
    parts <- unlist(strsplit(x, "\n", fixed = TRUE), use.names = FALSE)
    if (!length(parts)) return(0)
    max(systemfonts::string_width(parts, family = family, size = size,
                                  res = 72, weight = weight)) / 72
}

.n_lines <- function(x) {
    vapply(strsplit(as.character(x), "\n", fixed = TRUE), length, integer(1))
}

.first_non_na <- function(x) {
    x <- x[!is.na(x) & nzchar(x)]
    if (!length(x)) NA_character_ else x[[1]]
}

.check_pkg <- function(pkg, what = pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        stop("Package '", pkg, "' is required for ", what, ".", call. = FALSE)
    }
}

.empty_enrich_result <- function() {
    data.frame(ID = character(), Description = character(),
               GeneRatio = character(), BgRatio = character(),
               pvalue = numeric(), p.adjust = numeric(), qvalue = numeric(),
               geneID = character(), Count = integer(),
               FoldEnrichment = numeric(), RichFactor = numeric(),
               BgCount = integer(), stringsAsFactors = FALSE)
}

.ratio_parts <- function(x) {
    sp <- strsplit(as.character(x), "/", fixed = TRUE)
    num <- suppressWarnings(as.numeric(vapply(sp, `[`, character(1), 1)))
    den <- suppressWarnings(as.numeric(vapply(sp, function(v) {
        if (length(v) >= 2) v[[2]] else NA_character_
    }, character(1))))
    list(num = num, den = den)
}

.node_boxes <- function(nodes, fontsize = 7, lineheight = 0.9,
                           family = "Arial", node_pad = 0.06) {
    if (!"label_l1" %in% names(nodes)) {
        nodes$label_l1 <- nodes$id
        nodes$label_l2 <- .wrap_label(nodes$name, width = 25, max_lines = Inf)
        nodes$label_l3 <- ""
    }
    nodes$n_lines <- 1L + .n_lines(nodes$label_l2) + 1L
    nodes$w <- vapply(seq_len(nrow(nodes)), function(i) {
        .text_width(c(nodes$label_l1[i], strsplit(nodes$label_l2[i], "\n")[[1]],
                      nodes$label_l3[i]),
                    size = fontsize, family = family) + 2 * node_pad
    }, numeric(1))
    nodes$h <- nodes$n_lines * fontsize * lineheight / 72 + 2 * node_pad
    nodes
}

.html_escape <- function(x) {
    x <- gsub("&", "&amp;", x, fixed = TRUE)
    x <- gsub("<", "&lt;", x, fixed = TRUE)
    x <- gsub(">", "&gt;", x, fixed = TRUE)
    gsub("\"", "&quot;", x, fixed = TRUE)
}

.timed <- function(expr, verbose, label) {
    if (!verbose) return(expr)
    t0 <- Sys.time()
    on.exit(message(sprintf("  %-28s %.1f s", label,
                            as.numeric(difftime(Sys.time(), t0, units = "secs")))),
            add = TRUE)
    force(expr)
}
.font_family <- function(family) {
    if (is.null(family) || !nzchar(family)) return("sans")
    if (identical(Sys.info()[["sysname"]], "Darwin") &&
        identical(family, "Arial")) return("sans")
    if (requireNamespace("systemfonts", quietly = TRUE)) {
        fonts <- tryCatch(systemfonts::system_fonts(), error = function(e) NULL)
        if (!is.null(fonts) && "family" %in% names(fonts) &&
            any(tolower(fonts$family) == tolower(family))) {
            return(family)
        }
    }
    "sans"
}