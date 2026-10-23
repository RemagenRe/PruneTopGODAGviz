#' Propagate GO annotations to ancestors
#'
#' Applies the true-path rule: every transcript annotated to a term is also
#' annotated to all ancestors of that term.  This is the convention used by
#' agriGO and similar tools, and is required when the input
#' table contains direct annotations only.
#'
#' @param annotation A `go_annotation` object.
#' @param dag A `go_dag` object, usually for the same ontology.
#' @param ont Optional ontology to propagate; `NULL` propagates all three.
#' @param verbose Print a summary.
#'
#' @return A `go_annotation` object with expanded `pairs` and `terms`.
#' @export
propagate_go_annotation <- function(annotation, dag, ont = NULL, verbose = TRUE) {
    if (!inherits(annotation, "go_annotation")) {
        stop("'annotation' must be a go_annotation object.", call. = FALSE)
    }
    if (!inherits(dag, "go_dag")) {
        stop("'dag' must be a go_dag object.", call. = FALSE)
    }
    onts <- if (is.null(ont)) c("BP", "MF", "CC") else .check_ont(ont)
    out_pairs <- list()
    out_terms <- list()
    for (o in onts) {
        gos <- unique(annotation$pairs$go_id[
            annotation$pairs$go_id %in%
                annotation$terms$go_id[annotation$terms$ontology == o]])
        gos <- intersect(gos, dag$nodes$id)
        if (!length(gos)) next
        anc <- lapply(gos, function(g) dag_ancestors(dag, g))
        names(anc) <- gos
        pairs <- annotation$pairs[annotation$pairs$go_id %in% gos, ,
                                  drop = FALSE]
        expanded <- vector("list", length(gos))
        for (i in seq_along(gos)) {
            g <- gos[i]
            ts <- pairs$transcript[pairs$go_id == g]
            if (!length(ts)) next
            a <- anc[[g]]
            expanded[[i]] <- data.frame(
                transcript = rep(ts, each = length(a)),
                go_id = rep(a, times = length(ts)),
                stringsAsFactors = FALSE)
        }
        expanded <- unique(do.call(rbind, expanded))
        nodes <- dag$nodes[dag$nodes$id %in% unique(expanded$go_id),
                           c("id", "name", "namespace")]
        out_pairs[[o]] <- expanded
        out_terms[[o]] <- data.frame(
            go_id = nodes$id, name = nodes$name,
            ontology = nodes$namespace, def = NA_character_,
            stringsAsFactors = FALSE)
    }
    if (!length(out_pairs)) {
        return(annotation)
    }
    pairs <- unique(do.call(rbind, out_pairs))
    terms <- unique(do.call(rbind, out_terms))
    rownames(pairs) <- NULL
    rownames(terms) <- NULL
    universe <- lapply(c("BP", "MF", "CC"), function(o) {
        g <- terms$go_id[terms$ontology == o]
        sort(unique(pairs$transcript[pairs$go_id %in% g]))
    })
    names(universe) <- c("BP", "MF", "CC")
    out <- list(pairs = pairs, terms = terms, universe = universe,
                keytype = annotation$keytype, source = annotation$source,
                propagated = TRUE)
    class(out) <- "go_annotation"
    if (verbose) {
        message(sprintf("Propagated annotation: %d pairs, %d terms",
                        nrow(pairs), nrow(terms)))
    }
    out
}
