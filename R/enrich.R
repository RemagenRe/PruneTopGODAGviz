#' GO over-representation analysis
#'
#' A self-contained over-representation test for transcript-level GO annotation
#' tables produced by [read_go_annotation()].  Terms are tested with the
#' hypergeometric distribution, p-values are adjusted with `p.adjust()`, and
#' effect sizes are reported in the same columns used by the rest of the
#' package.  No external enrichment package is required.
#'
#' @param gene Character vector of input identifiers (transcripts, genes, ...).
#' @param annotation A `go_annotation` object from [read_go_annotation()].
#' @param ont One of `"BP"`, `"MF"` or `"CC"`.
#' @param universe Background identifiers.  Defaults to all identifiers
#'   annotated in the selected ontology.
#' @param pvalueCutoff,qvalueCutoff Stored in the result metadata for backwards
#'   compatibility.  The returned table contains all tested terms; filtering is
#'   performed by [prune_go_dag()].
#' @param pAdjustMethod Method passed to [stats::p.adjust()].
#' @param minGSSize,maxGSSize Minimum and maximum gene-set size.
#' @param keytype Optional label stored in the result (e.g. `"transcript"`).
#' @param verbose Print a short summary.
#'
#' @return A `go_enrich` object: a list with `result` (data frame), `enrich`
#'   (internal result metadata), `gene`, `universe`, `ontology`, `annotation`
#'   and `params`.
#'
#' @examples
#' \dontrun{
#' ego <- enrich_go(my_genes, anno, ont = "BP")
#' head(ego$result)
#' }
#' @export
enrich_go <- function(gene, annotation, ont = "BP", universe = NULL,
                      pvalueCutoff = 0.05, pAdjustMethod = "BH",
                      qvalueCutoff = 0.05, minGSSize = 10, maxGSSize = 500,
                      keytype = NULL, verbose = TRUE) {
    if (!inherits(annotation, "go_annotation")) {
        stop("'annotation' must be a go_annotation object.", call. = FALSE)
    }
    ont <- .check_ont(ont)
    gene <- unique(trimws(as.character(gene)))
    gene <- gene[!is.na(gene) & nzchar(gene)]
    if (!length(gene)) stop("'gene' is empty.", call. = FALSE)

    uni <- universe %||% annotation$universe[[ont]]
    uni <- unique(trimws(as.character(uni)))
    uni <- uni[!is.na(uni) & nzchar(uni)]
    if (!length(uni)) {
        stop("No background identifiers annotated in ontology '", ont, "'.",
             call. = FALSE)
    }

    terms <- annotation$terms[annotation$terms$ontology == ont, , drop = FALSE]
    if (!nrow(terms)) {
        stop("No GO terms annotated in ontology '", ont, "'.")
    }
    pairs <- annotation$pairs[annotation$pairs$go_id %in% terms$go_id, ,
                              drop = FALSE]

    hit <- intersect(gene, uni)
    if (!length(hit)) {
        stop("None of the input identifiers is annotated in ontology '", ont,
             "' (checked against ", length(uni), " background identifiers).",
             call. = FALSE)
    }
    if (length(hit) < length(gene)) {
        warning(length(gene) - length(hit), " of ", length(gene),
                " input identifiers are not annotated and are ignored.",
                call. = FALSE)
    }

    t2g <- unique(pairs[, c("go_id", "transcript")])
    t2n <- unique(terms[, c("go_id", "name")])
    if (verbose) {
        message(sprintf("Enrichment: %d input / %d background IDs, %d terms (%s)",
                        length(hit), length(uni), nrow(t2g), ont))
    }

    ego <- .ora_enrich(hit = hit, universe = uni, t2g = t2g, t2n = t2n,
                       pAdjustMethod = pAdjustMethod,
                       minGSSize = minGSSize, maxGSSize = maxGSSize)
    if (is.null(ego)) {
        warning("No GO term passed the size or significance filter.",
                call. = FALSE)
        res <- .empty_enrich_result()
        obj <- NULL
    } else {
        res <- .augment_enrich_result(ego$result)
        res$qvalue <- .storey_qvalue(res$pvalue)
        ego$result <- res
        ego$ontology <- ont
        ego$keytype <- keytype %||% annotation$keytype %||% "transcript"
        obj <- ego
    }

    gene_out <- if (!is.null(ego)) ego$gene else hit
    uni_out <- if (!is.null(ego)) ego$universe else uni
    out <- list(result = res, enrich = obj, gene = gene_out, universe = uni_out,
                ontology = ont, annotation = annotation,
                params = list(pvalueCutoff = pvalueCutoff,
                              pAdjustMethod = pAdjustMethod,
                              qvalueCutoff = qvalueCutoff,
                              minGSSize = minGSSize, maxGSSize = maxGSSize,
                              keytype = keytype %||% annotation$keytype))
    class(out) <- "go_enrich"
    out
}

.storey_qvalue <- function(p, lambda = 0.5) {
    p <- as.numeric(p)
    n <- length(p)
    if (!n) return(numeric())
    pi0 <- mean(p >= lambda, na.rm = TRUE) / (1 - lambda)
    pi0 <- min(1, pi0)
    if (!is.finite(pi0) || pi0 <= 0) pi0 <- 1
    ord <- order(p)
    q_sorted <- pi0 * n * p[ord] / seq_len(n)
    q_sorted <- rev(cummin(rev(q_sorted)))
    q <- numeric(n)
    q[ord] <- q_sorted
    pmin(pmax(q, 0), 1)
}
.ora_enrich <- function(hit, universe, t2g, t2n, pAdjustMethod,
                        minGSSize, maxGSSize) {
    all_ids <- unique(as.character(t2g$transcript))
    ext_id <- intersect(all_ids, as.character(universe))
    if (!length(ext_id)) return(NULL)

    q <- unique(t2g[t2g$transcript %in% hit, , drop = FALSE])
    q_terms <- unique(q$go_id)
    if (!length(q_terms)) return(NULL)
    n_mapped <- length(unique(as.character(q$transcript)))

    sets <- split(as.character(t2g$transcript), as.character(t2g$go_id))
    sets <- sets[q_terms]
    sets <- lapply(sets, intersect, ext_id)
    sizes <- lengths(sets)
    keep <- sizes >= minGSSize & sizes <= maxGSSize
    if (!any(keep)) return(NULL)
    sets <- sets[keep]
    q_terms <- names(sets)
    q_sets <- lapply(q_terms, function(g) {
        unique(intersect(q$transcript[q$go_id == g], ext_id))
    })
    names(q_sets) <- q_terms
    k <- lengths(q_sets)
    M <- lengths(sets)[q_terms]
    N <- rep(length(ext_id), length(M))
    n <- rep(n_mapped, length(M))

    pvalue <- stats::phyper(k - 1, M, N - M, n, lower.tail = FALSE)
    gene_ratio <- sprintf("%s/%s", k, n)
    bg_ratio <- sprintf("%s/%s", M, N)
    rich_factor <- k / M
    fold <- rich_factor * N / n
    mu <- M * n / N
    sigma <- mu * (N - n) * (N - M) / N / (N - 1)
    zscore <- (k - mu) / sqrt(sigma)
    gene_id <- vapply(q_sets, paste, character(1), collapse = "/")
    desc <- t2n$name[match(q_terms, t2n$go_id)]
    desc[is.na(desc)] <- q_terms[is.na(desc)]

    res <- data.frame(
        ID = as.character(q_terms),
        Description = as.character(desc),
        GeneRatio = gene_ratio,
        BgRatio = bg_ratio,
        RichFactor = as.numeric(rich_factor),
        FoldEnrichment = as.numeric(fold),
        zScore = as.numeric(zscore),
        pvalue = as.numeric(pvalue),
        p.adjust = stats::p.adjust(pvalue, method = pAdjustMethod),
        qvalue = .storey_qvalue(pvalue),
        geneID = as.character(gene_id),
        Count = as.integer(k),
        stringsAsFactors = FALSE
    )
    res <- res[order(res$pvalue), , drop = FALSE]
    rownames(res) <- res$ID

    list(result = res, gene = unique(as.character(q$transcript)),
         universe = ext_id,
         pvalueCutoff = NA_real_, pAdjustMethod = pAdjustMethod,
         qvalueCutoff = NA_real_, geneSets = sets,
         keytype = "UNKNOWN", ontology = "UNKNOWN")
}

#' Re-use an existing over-representation result
#'
#' Wraps a result data frame or a list containing a `result` data frame so that
#' it can be drawn with [prune_go_dag()] without re-running the enrichment.
#'
#' @param x A data frame with `ID`, `Description`, `pvalue`, `p.adjust`,
#'   `GeneRatio` and `BgRatio`, or a list containing such a data frame in
#'   `x$result`.
#' @param annotation A `go_annotation` object used for gene counts.
#' @param ont Ontology; defaults to the supplied result or `"BP"`.
#' @param universe Background identifiers; defaults to the annotation universe.
#' @param gene Input identifiers.
#' @return A `go_enrich` object.
#' @export
as_go_enrich <- function(x, annotation, ont = NULL, universe = NULL, gene = NULL) {
    if (is.data.frame(x)) {
        res <- x
        ont0 <- NULL
        gene0 <- NULL
    } else if (is.list(x) && !is.null(x$result)) {
        res <- x$result
        ont0 <- x$ontology
        gene0 <- x$gene
    } else {
        stop("'x' must be a result data frame or a list with a 'result' data frame.",
             call. = FALSE)
    }
    ont <- .check_ont(ont %||% ont0 %||% "BP")
    uni <- universe %||% annotation$universe[[ont]]
    res <- .augment_enrich_result(res)
    out <- list(result = res, enrich = x,
                gene = gene %||% gene0 %||% character(),
                universe = uni, ontology = ont,
                annotation = annotation,
                params = list(keytype = annotation$keytype))
    class(out) <- "go_enrich"
    out
}

.augment_enrich_result <- function(res) {
    if (!nrow(res)) {
        for (col in c("FoldEnrichment", "RichFactor", "BgCount")) {
            if (!col %in% names(res)) res[[col]] <- numeric()
        }
        return(res)
    }
    gr <- .ratio_parts(res$GeneRatio)
    br <- .ratio_parts(res$BgRatio)
    fe <- (gr$num / gr$den) / (br$num / br$den)
    res$FoldEnrichment <- as.numeric(fe)
    res$RichFactor <- as.numeric(gr$num / br$num)
    res$BgCount <- as.integer(br$num)
    if ("pvalue" %in% names(res)) {
        q_all <- .storey_qvalue(res$pvalue)
        if (!"qvalue" %in% names(res)) {
            res$qvalue <- q_all
        } else {
            fill <- is.na(res$qvalue) & !is.na(res$pvalue)
            res$qvalue[fill] <- q_all[fill]
        }
    }
    res
}

#' @export
print.go_enrich <- function(x, n = 5L, ...) {
    cat("<go_enrich> ontology:", x$ontology, "\n")
    cat("  input:", length(x$gene), " background:", length(x$universe), "\n")
    cat("  significant terms:", nrow(x$result), "\n")
    if (nrow(x$result)) {
        show <- utils::head(x$result[, c("ID", "Description", "GeneRatio",
                                         "pvalue", "p.adjust")], n)
        print(show, row.names = FALSE)
    }
    invisible(x)
}

#' @export
as.data.frame.go_enrich <- function(x, ...) x$result
