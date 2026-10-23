.dag_graph <- function(dag) {
    if (is.null(dag$edges) || !nrow(dag$edges)) {
        return(igraph::graph_from_data_frame(
            data.frame(parent = character(), child = character(),
                       stringsAsFactors = FALSE),
            directed = TRUE,
            vertices = data.frame(name = dag$nodes$id, stringsAsFactors = FALSE)))
    }
    igraph::graph_from_data_frame(
        as.data.frame(dag$edges[, c("parent", "child")]),
        directed = TRUE,
        vertices = data.frame(name = dag$nodes$id, stringsAsFactors = FALSE))
}

#' Ancestors of GO terms
#'
#' @param dag A `go_dag` or `pruned_go_dag` object.
#' @param ids GO identifiers.
#' @param include_self Include the query terms themselves.
#' @return A character vector of GO identifiers.
#' @export
dag_ancestors <- function(dag, ids, include_self = TRUE) {
    g <- .dag_graph(dag)
    idx <- match(unique(ids), igraph::V(g)$name)
    idx <- idx[!is.na(idx)]
    if (!length(idx)) return(character())
    v <- unique(unlist(lapply(idx, function(i) {
        as.integer(igraph::subcomponent(g, i, mode = "in"))
    })))
    out <- igraph::V(g)$name[v]
    if (!include_self) out <- setdiff(out, ids)
    out
}

#' Descendants of GO terms
#'
#' @inheritParams dag_ancestors
#' @return A character vector of GO identifiers.
#' @export
dag_descendants <- function(dag, ids, include_self = TRUE) {
    g <- .dag_graph(dag)
    idx <- match(unique(ids), igraph::V(g)$name)
    idx <- idx[!is.na(idx)]
    if (!length(idx)) return(character())
    v <- unique(unlist(lapply(idx, function(i) {
        as.integer(igraph::subcomponent(g, i, mode = "out"))
    })))
    out <- igraph::V(g)$name[v]
    if (!include_self) out <- setdiff(out, ids)
    out
}

.dag_levels <- function(dag) {
    g <- .dag_graph(dag)
    n <- igraph::vcount(g)
    lev <- rep(0L, n)
    if (n == 0L) return(setNames(integer(), character()))
    if (igraph::ecount(g) > 0L && !igraph::is_dag(g)) {
        warning("GO graph contains a cycle; levels may be unreliable.",
                call. = FALSE)
        ord <- seq_len(n)
    } else {
        ord <- as.integer(igraph::topo_sort(g, mode = "out"))
    }
    ed <- igraph::as_edgelist(g, names = FALSE)
    if (nrow(ed)) {
        rank <- integer(n)
        rank[ord] <- seq_along(ord)
        o <- order(rank[ed[, 1]])
        for (k in o) {
            p <- ed[k, 1]
            ch <- ed[k, 2]
            if (lev[p] + 1L > lev[ch]) lev[ch] <- lev[p] + 1L
        }
    }
    setNames(lev, igraph::V(g)$name)
}

.dag_ancestor_distance <- function(dag, ids) {
    g <- .dag_graph(dag)
    idx <- match(unique(ids), igraph::V(g)$name)
    idx <- idx[!is.na(idx)]
    if (!length(idx)) return(setNames(numeric(), character()))
    d <- igraph::distances(g, v = idx, mode = "in")
    out <- apply(d, 2, min)
    out[!is.finite(out)] <- Inf
    setNames(out, igraph::V(g)$name)
}

.dag_subset <- function(dag, ids) {
    ids <- unique(ids)
    list(nodes = dag$nodes[dag$nodes$id %in% ids, , drop = FALSE],
         edges = dag$edges[dag$edges$parent %in% ids & dag$edges$child %in% ids, ,
                           drop = FALSE],
         roots = intersect(dag$roots, ids),
         source = dag$source, version = dag$version,
         file = dag$file)
}
