#' Prune and annotate a GO DAG for display
#'
#' Selects enriched terms, adds the ancestors required for connectivity, and
#' removes or collapses the parts of the graph that make the figure unreadable.
#' Pruning never changes enrichment statistics: p-values and ratios are copied
#' verbatim from the `go_enrich` object.
#'
#' @param x A `go_enrich` object from [enrich_go()] or [as_go_enrich()].
#' @param dag A `go_dag` object.  Defaults to [go_dag()] for `x$ontology`.
#' @param top_n Maximum number of significant terms to display (by `stat`).
#'   Use `Inf` to keep all terms passing `p_cutoff`.
#' @param p_cutoff Significance threshold applied to `stat`.
#' @param stat Column used for ranking and colouring: `"p.adjust"` or `"pvalue"`.
#' @param keep GO identifiers forced into the figure (shown as selected).
#' @param drop GO identifiers removed from the figure.
#' @param drop_subtree GO identifiers whose entire subtree is removed.
#' @param keep_branches GO identifiers whose subtrees are the only branches kept.
#' @param exclude_regex Regular expression; matching term names are removed.
#' @param include_regex Regular expression; only matching terms (plus selected
#'   terms and their ancestors) are kept.
#' @param keep_ancestors Add ancestors of selected terms.
#' @param max_depth Maximum number of steps above a selected term.  `Inf` keeps
#'   the full ancestry up to the ontology root.
#' @param expand_depth Number of levels below selected terms to keep.
#' @param collapse_chains Collapse linear, non-selected ancestor chains.
#' @param collapse_min_length Minimum chain length (in terms) to collapse.
#' @param cluster_jaccard Merge sibling leaves whose gene sets have a Jaccard
#'   index of at least this value.  `NULL` disables merging.
#' @param max_nodes Hard cap on the number of displayed terms.
#' @param wrap Number of characters per line when wrapping term names.
#' @param counts `"propagated"` (default) counts input and background genes
#'   annotated to the term or any descendant (true-path rule, as in agriGO);
#'   `"direct"` counts only genes annotated to the term itself.
#' @param verbose Print a summary.
#'
#' @return A `pruned_go_dag` object with `nodes`, `edges`, `dropped`, `params`
#'   and `ontology`.
#'
#' @details
#' Deleting a term removes only descendants that lose every parent; children
#' shared with another parent are kept.
#'
#' The third label line follows the true-path rule by default: it counts input
#' and background genes annotated to the term or to any of its descendants.
#' This is the convention used by agriGO and avoids the misleading zeros that
#' direct annotation produces for high-level terms.  The enrichment table keeps
#' the direct counts (and statistics) returned by `enrich_go()`; use
#' `counts = "direct"` to show the same direct counts in the figure.
#'
#' @examples
#' \dontrun{
#' dag <- prune_go_dag(ego, top_n = 8, drop_subtree = "GO:0006950")
#' }
#' @export
prune_go_dag <- function(x, dag = NULL, top_n = 10, p_cutoff = 0.05,
                         stat = c("p.adjust", "pvalue"), keep = NULL,
                         drop = NULL, drop_subtree = NULL,
                         keep_branches = NULL, exclude_regex = NULL,
                         include_regex = NULL, keep_ancestors = TRUE,
                         max_depth = Inf, expand_depth = 0,
                         collapse_chains = FALSE, collapse_min_length = 2,
                         cluster_jaccard = NULL, max_nodes = 80,
                         wrap = 25, counts = c("propagated", "direct"),
                         verbose = TRUE) {
    if (!inherits(x, "go_enrich")) {
        stop("'x' must be a go_enrich object.", call. = FALSE)
    }
    stat <- match.arg(stat)
    counts <- match.arg(counts)
    ont <- x$ontology
    dag <- dag %||% go_dag(ont, verbose = verbose)
    res <- x$result
    if (!stat %in% names(res)) {
        stop("Column '", stat, "' not found in the enrichment result.",
             call. = FALSE)
    }
    res <- res[!is.na(res[[stat]]), , drop = FALSE]
    if (!nrow(res)) stop("The enrichment result is empty.", call. = FALSE)
    res <- res[order(res[[stat]]), , drop = FALSE]
    stat_map <- stats::setNames(res[[stat]], res$ID)

    sel_df <- res[res[[stat]] <= p_cutoff, , drop = FALSE]
    if (is.finite(top_n)) sel_df <- utils::head(sel_df, top_n)
    selected <- unique(c(sel_df$ID, keep))
    unknown <- setdiff(selected, dag$nodes$id)
    if (length(unknown)) {
        warning("Ignoring ", length(unknown), " identifier(s) absent from the ontology: ",
                paste(utils::head(unknown, 3), collapse = ", "), call. = FALSE)
    }
    selected <- intersect(selected, dag$nodes$id)
    if (!length(selected)) {
        stop("No term passes p_cutoff = ", p_cutoff, " (", stat, ").",
             call. = FALSE)
    }

    reason <- stats::setNames(rep("off_path", nrow(dag$nodes)), dag$nodes$id)
    reason[selected] <- "kept"

    branch_roots <- character()
    if (!is.null(keep_branches)) {
        branch_roots <- intersect(keep_branches, dag$nodes$id)
        if (!length(branch_roots)) {
            stop("None of 'keep_branches' exists in the ontology.", call. = FALSE)
        }
        allowed <- union(branch_roots, dag_descendants(dag, branch_roots))
        outside <- setdiff(selected, allowed)
        if (length(outside)) {
            warning("Dropping ", length(outside),
                    " selected term(s) outside 'keep_branches'.", call. = FALSE)
            reason[outside] <- "outside_branch"
        }
        selected <- intersect(selected, allowed)
        if (!length(selected)) {
            stop("No selected term remains inside 'keep_branches'.", call. = FALSE)
        }
    }

    cand <- selected
    if (isTRUE(keep_ancestors)) {
        anc <- dag_ancestors(dag, selected)
        if (is.finite(max_depth)) {
            dist <- .dag_ancestor_distance(dag, selected)
            anc <- anc[is.finite(dist[anc]) & dist[anc] <= max_depth]
        }
        cand <- union(cand, anc)
    }
    if (is.finite(expand_depth) && expand_depth > 0) {
        lev <- .dag_levels(dag)
        desc <- dag_descendants(dag, selected)
        max_lev <- max(lev[selected], na.rm = TRUE) + expand_depth
        cand <- union(cand, desc[lev[desc] <= max_lev])
    }
    cand <- union(cand, branch_roots)

    if (!is.null(drop_subtree)) {
        rem <- .cascade_orphans(dag, drop_subtree)
        hit <- intersect(rem, selected)
        if (length(hit)) {
            warning("'drop_subtree' removes selected term(s): ",
                    paste(utils::head(hit, 3), collapse = ", "),
                    call. = FALSE)
        }
        cand <- setdiff(cand, rem)
        reason[intersect(rem, names(reason))] <- "drop_subtree"
        selected <- setdiff(selected, rem)
    }
    if (!length(selected)) {
        stop("All selected terms were removed by 'drop_subtree'.", call. = FALSE)
    }
    if (!is.null(drop)) {
        rem <- intersect(drop, dag$nodes$id)
        hit <- intersect(rem, selected)
        if (length(hit)) {
            warning("'drop' removes selected term(s): ",
                    paste(utils::head(hit, 3), collapse = ", "), call. = FALSE)
        }
        cand <- setdiff(cand, rem)
        selected <- setdiff(selected, rem)
        reason[rem] <- "explicit"
    }
    if (!is.null(exclude_regex)) {
        hit <- dag$nodes$id[grepl(exclude_regex, dag$nodes$name, ignore.case = TRUE)]
        cand <- setdiff(cand, hit)
        selected <- setdiff(selected, hit)
        reason[intersect(hit, names(reason))] <- "exclude_regex"
    }
    if (!is.null(include_regex)) {
        must <- union(selected, dag_ancestors(dag, selected))
        hit <- dag$nodes$id[grepl(include_regex, dag$nodes$name, ignore.case = TRUE)]
        keep_re <- intersect(cand, union(must, hit))
        reason[setdiff(cand, keep_re)] <- "exclude_regex"
        cand <- keep_re
    }
    if (!length(selected)) {
        stop("All selected terms were removed by the filters.", call. = FALSE)
    }

    sub <- .dag_subset(dag, cand)
    needed <- union(selected, dag_ancestors(sub, selected))
    keep_final <- intersect(cand, union(needed, branch_roots))
    reason[setdiff(cand, keep_final)] <- "off_path"

    nodes <- dag$nodes[dag$nodes$id %in% keep_final, , drop = FALSE]
    edges <- dag$edges[dag$edges$parent %in% keep_final &
                           dag$edges$child %in% keep_final, , drop = FALSE]
    edges$n_skipped <- 0L
    n_collapsed <- stats::setNames(integer(nrow(nodes)), nodes$id)
    collapsed_all <- list()

    if (isTRUE(collapse_chains) && nrow(nodes) > 2) {
        res_c <- .collapse_chains(nodes, edges, selected, branch_roots,
                                  collapse_min_length)
        gone <- setdiff(nodes$id, res_c$nodes$id)
        reason[gone] <- "chain_collapsed"
        nodes <- res_c$nodes
        edges <- res_c$edges
        n_collapsed <- res_c$n_collapsed
        for (nm in names(res_c$collapsed)) {
            collapsed_all[[nm]] <- c(collapsed_all[[nm]], res_c$collapsed[[nm]])
        }
    }

    if (!is.null(cluster_jaccard) && is.finite(cluster_jaccard) &&
        nrow(nodes) > 2) {
        res_j <- .cluster_jaccard(nodes, edges, selected, x$annotation,
                                  stat_map, cluster_jaccard)
        nodes <- res_j$nodes
        edges <- res_j$edges
        n_collapsed <- n_collapsed[names(n_collapsed) %in% nodes$id]
        for (nm in names(res_j$n_collapsed)) {
            n_collapsed[nm] <- (if (nm %in% names(n_collapsed)) n_collapsed[[nm]] else 0L) +
                res_j$n_collapsed[[nm]]
        }
        reason[names(res_j$merged)] <- res_j$merged
        for (nm in names(res_j$collapsed)) {
            collapsed_all[[nm]] <- c(collapsed_all[[nm]], res_j$collapsed[[nm]])
        }
    }

    if (is.finite(max_nodes) && nrow(nodes) > max_nodes) {
        res_m <- .cap_nodes(nodes, edges, selected, stat_map, max_nodes)
        nodes <- res_m$nodes
        edges <- res_m$edges
        reason[names(res_m$dropped)] <- "max_nodes"
    }

    lev <- .dag_levels(list(nodes = nodes, edges = edges))
    nodes$level <- as.integer(lev[nodes$id])
    nodes$selected <- nodes$id %in% selected
    nodes$n_collapsed <- as.integer(n_collapsed[nodes$id])
    nodes$n_collapsed[is.na(nodes$n_collapsed)] <- 0L
    nodes$collapsed_ids <- vapply(nodes$id, function(id) {
        paste(unique(collapsed_all[[id]]), collapse = ";")
    }, character(1))
    nodes$collapsed_terms <- vapply(nodes$collapsed_ids, function(ids) {
        ids <- strsplit(ids, ";", fixed = TRUE)[[1]]
        ids <- ids[nzchar(ids)]
        if (!length(ids)) return("")
        paste(dag$nodes$name[match(ids, dag$nodes$id)], collapse = "; ")
    }, character(1))

    sets <- split(x$annotation$pairs$transcript, x$annotation$pairs$go_id)
    set_size <- function(id, pool) {
        g <- sets[[id]]
        if (is.null(g)) return(0L)
        length(intersect(g, pool))
    }
    if (counts == "propagated") {
        for (i in seq_len(nrow(nodes))) {
            ids <- dag_descendants(dag, nodes$id[i])
            pool <- unique(unlist(sets[ids], use.names = FALSE))
            nodes$n_in_term[i] <- length(intersect(pool, x$gene))
            nodes$n_bg_term[i] <- length(intersect(pool, x$universe))
        }
    } else {
        nodes$n_in_term <- vapply(nodes$id, set_size, integer(1), pool = x$gene)
        nodes$n_bg_term <- vapply(nodes$id, set_size, integer(1), pool = x$universe)
    }
    nodes$n_in_universe <- length(x$gene)
    nodes$n_bg_universe <- length(x$universe)
    nodes$pvalue <- unname(x$result$pvalue[match(nodes$id, x$result$ID)])
    nodes$p.adjust <- unname(x$result$p.adjust[match(nodes$id, x$result$ID)])

    show_stat <- nodes[[stat]]
    nodes$label_l1 <- nodes$id
    has <- nodes$selected & !is.na(show_stat)
    nodes$label_l1[has] <- paste0(nodes$id[has], " (", .fmt_p(show_stat[has]), ")")
    nodes$label_l2 <- .wrap_label(nodes$name, width = wrap, max_lines = Inf)
    add <- nodes$n_collapsed > 0
    nodes$label_l2[add] <- .wrap_label(
        paste0(nodes$name[add], " (+", nodes$n_collapsed[add], " hidden)"),
        width = wrap, max_lines = Inf)
    nodes$label_l3 <- paste0(nodes$n_in_term, "/", nodes$n_in_universe, " | ",
                             nodes$n_bg_term, "/", nodes$n_bg_universe)

    ordered <- match(dag$nodes$id, nodes$id, nomatch = 0L)
    miss <- dag$nodes$id[ordered == 0L]
    dropped <- data.frame(
        go_id = miss,
        name = dag$nodes$name[ordered == 0L],
        p.adjust = unname(x$result$p.adjust[match(miss, x$result$ID)]),
        action = unname(reason[miss]),
        reason = unname(reason[miss]),
        stringsAsFactors = FALSE)
    dropped <- dropped[order(dropped$action, dropped$go_id), , drop = FALSE]
    rownames(dropped) <- NULL

    out <- list(nodes = nodes, edges = edges, dropped = dropped,
                annotation = x$annotation, gene = x$gene,
                universe = x$universe, ontology = ont, stat = stat,
                params = list(top_n = top_n, p_cutoff = p_cutoff, stat = stat,
                              keep = keep, drop = drop,
                              drop_subtree = drop_subtree,
                              keep_branches = keep_branches,
                              exclude_regex = exclude_regex,
                              include_regex = include_regex,
                              max_depth = max_depth,
                              expand_depth = expand_depth,
                              collapse_chains = collapse_chains,
                              cluster_jaccard = cluster_jaccard,
                              max_nodes = max_nodes, wrap = wrap,
                              counts = counts))
    class(out) <- "pruned_go_dag"
    if (verbose) {
        message(sprintf("Pruned DAG: %d shown (%d selected), %d dropped",
                        nrow(nodes), sum(nodes$selected), nrow(dropped)))
    }
    out
}

.collapse_chains <- function(nodes, edges, selected, branch_roots, min_length) {
    n_collapsed <- stats::setNames(integer(nrow(nodes)), nodes$id)
    collapsed <- list()
    if (!nrow(edges)) {
        return(list(nodes = nodes, edges = edges, n_collapsed = n_collapsed,
                    collapsed = collapsed))
    }
    repeat {
        indeg <- table(factor(edges$child, levels = nodes$id))
        outdeg <- table(factor(edges$parent, levels = nodes$id))
        elig <- nodes$id[!(nodes$id %in% selected) &
                             indeg[nodes$id] == 1 & outdeg[nodes$id] == 1 &
                             !(nodes$id %in% branch_roots)]
        if (!length(elig)) break
        parent_of <- stats::setNames(edges$parent[match(elig, edges$child)], elig)
        tops <- elig[!(parent_of[elig] %in% elig)]
        if (!length(tops)) break
        changed <- FALSE
        for (v in tops) {
            if (!(v %in% nodes$id)) next
            chain <- v
            cur <- v
            repeat {
                ch <- edges$child[edges$parent == cur]
                if (length(ch) != 1L) break
                if (!(ch %in% elig)) break
                chain <- c(chain, ch)
                cur <- ch
            }
            if (length(chain) < max(2L, min_length)) next
            removed <- chain[-1]
            last <- chain[length(chain)]
            last_child <- edges$child[edges$parent == last]
            p <- edges$parent[edges$child == v]
            etype <- edges$type[edges$parent == p & edges$child == v]
            if (!length(p) || !length(last_child) || !length(etype)) next
            nodes <- nodes[!nodes$id %in% removed, , drop = FALSE]
            edges <- edges[!edges$parent %in% removed & !edges$child %in% removed, ,
                           drop = FALSE]
            n_collapsed[v] <- length(removed)
            collapsed[[v]] <- c(collapsed[[v]], removed)
            if (!any(edges$parent == v & edges$child == last_child)) {
                edges <- rbind(edges, data.frame(
                    parent = v, child = last_child, type = etype[[1]],
                    n_skipped = as.integer(length(removed)), stringsAsFactors = FALSE))
            }
            changed <- TRUE
        }
        if (!changed) break
    }
    list(nodes = nodes, edges = edges, n_collapsed = n_collapsed,
         collapsed = collapsed)
}

.cascade_orphans <- function(dag, roots) {
    roots <- intersect(unique(roots), dag$nodes$id)
    if (!length(roots)) return(character())
    kids <- split(dag$edges$child, dag$edges$parent)
    parents <- split(dag$edges$parent, dag$edges$child)
    removed <- roots
    repeat {
        children <- unique(unlist(kids[removed], use.names = FALSE))
        new <- character()
        for (v in children) {
            if (v %in% removed) next
            p <- parents[[v]]
            if (length(p) && all(p %in% removed)) new <- c(new, v)
        }
        if (!length(new)) break
        removed <- unique(c(removed, new))
    }
    removed
}

.jaccard <- function(a, b) {
    if (!length(a) || !length(b)) return(NA_real_)
    length(intersect(a, b)) / length(union(a, b))
}

.cluster_jaccard <- function(nodes, edges, selected, annotation, stat_map,
                             threshold) {
    merged <- character()
    n_collapsed <- integer()
    collapsed <- list()
    sets <- split(annotation$pairs$transcript, annotation$pairs$go_id)
    leaves <- nodes$id[!nodes$id %in% edges$parent]
    for (p in unique(edges$parent)) {
        ch <- intersect(edges$child[edges$parent == p], leaves)
        if (length(ch) < 2L) next
        ch <- ch[order(stat_map[ch], na.last = TRUE)]
        removed <- character()
        for (i in seq_along(ch)) {
            v <- ch[[i]]
            if (v %in% removed || v %in% selected) next
            for (j in seq_along(ch)) {
                if (i == j) next
                r <- ch[[j]]
                if (r %in% removed) next
                jac <- .jaccard(sets[[v]], sets[[r]])
                if (is.finite(jac) && jac >= threshold) {
                    removed <- c(removed, v)
                    merged[[v]] <- paste0("jaccard_merge:", r)
                    n_collapsed[r] <- (if (r %in% names(n_collapsed)) n_collapsed[[r]] else 0L) + 1L
                    collapsed[[r]] <- c(collapsed[[r]], v)
                    break
                }
            }
        }
        if (length(removed)) {
            nodes <- nodes[!nodes$id %in% removed, , drop = FALSE]
            edges <- edges[!edges$parent %in% removed & !edges$child %in% removed, ,
                           drop = FALSE]
        }
    }
    list(nodes = nodes, edges = edges, merged = merged,
         n_collapsed = n_collapsed, collapsed = collapsed)
}

.cap_nodes <- function(nodes, edges, selected, stat_map, max_nodes) {
    dropped <- character()
    lev <- .dag_levels(list(nodes = nodes, edges = edges))
    while (nrow(nodes) > max_nodes) {
        leaves <- nodes$id[!nodes$id %in% edges$parent]
        cand <- setdiff(leaves, selected)
        if (!length(cand)) break
        sc <- stat_map[cand]
        sc[is.na(sc)] <- 1
        v <- cand[order(-sc, -lev[cand])][[1]]
        nodes <- nodes[nodes$id != v, , drop = FALSE]
        edges <- edges[edges$parent != v & edges$child != v, , drop = FALSE]
        dropped[[v]] <- v
    }
    list(nodes = nodes, edges = edges, dropped = dropped)
}

#' Write the pruning decision table
#'
#' @param x A `pruned_go_dag` object.
#' @param file Output CSV path.
#' @return The path, invisibly.
#' @export
write_prune_table <- function(x, file) {
    if (!inherits(x, "pruned_go_dag")) {
        stop("'x' must be a pruned_go_dag object.", call. = FALSE)
    }
    kept <- data.frame(
        go_id = x$nodes$id, name = x$nodes$name, p.adjust = x$nodes$p.adjust,
        action = ifelse(x$nodes$selected, "kept_selected", "kept_ancestor"),
        reason = ifelse(x$nodes$n_collapsed > 0,
                        paste0("collapsed_", x$nodes$n_collapsed), ""),
        stringsAsFactors = FALSE)
    out <- rbind(kept, x$dropped)
    utils::write.csv(out, file, row.names = FALSE)
    invisible(file)
}

#' @export
print.pruned_go_dag <- function(x, ...) {
    cat("<pruned_go_dag> ontology:", x$ontology, "\n")
    cat("  shown   :", nrow(x$nodes), "terms (",
        sum(x$nodes$selected), "selected )\n")
    cat("  dropped :", nrow(x$dropped), "terms\n")
    cat("  edges   :", nrow(x$edges), "\n")
    invisible(x)
}

#' @export
as.data.frame.pruned_go_dag <- function(x, ...) {
    x$nodes[, c("id", "name", "level", "selected", "n_collapsed", "pvalue",
                "p.adjust", "n_in_term", "n_in_universe", "n_bg_term",
                "n_bg_universe")]
}
