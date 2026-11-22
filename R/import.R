#' Maximum accepted size for an imported external result file (50 MB).
#'
#' @keywords internal
.max_external_bytes <- 50 * 1024^2

#' Read an external GO over-representation result
#'
#' Import results produced by another enrichment tool so they can be edited,
#' pruned and exported without rerunning the analysis. The first supported
#' format is the PANTHER over-representation JSON download. A generic CSV/TSV
#' interface is also provided.
#'
#' External results are not re-tested. Their p-values, adjusted p-values and
#' counts are kept as supplied; the `qvalue` column is the supplied FDR, not a
#' newly calculated Storey q-value.
#'
#' @param file Path to a PANTHER JSON, CSV or TSV file.
#' @param ontology One of `"BP"`, `"MF"`, `"CC"` or `NULL` to detect it.
#' @param dag Optional `go_dag` object. If omitted, it is built from GO OBO.
#' @param obo Optional local GO OBO path.
#' @param max_bytes Maximum accepted file size in bytes. Files above this
#'   limit are rejected before parsing (default 50 MB).
#' @param verbose Print a short summary.
#' @return A `go_enrich` object with an `external_counts` component.
#' @export
import_go_enrich <- function(file, ontology = NULL, dag = NULL, obo = NULL,
                             max_bytes = .max_external_bytes, verbose = TRUE) {
    if (!file.exists(file)) stop("File not found: ", file, call. = FALSE)
    size <- file.info(file)$size
    if (is.finite(size) && size > max_bytes) {
        stop(sprintf(
            "External result file is %.1f MB; the limit is %s MB.",
            size / 1024^2, format(round(max_bytes / 1024^2), trim = TRUE)),
            call. = FALSE)
    }
    ext <- tolower(tools::file_ext(file))
    parsed <- if (identical(ext, "json")) {
        parse_panther_json(file)
    } else {
        parse_external_table(file, ontology = ontology)
    }
    ont <- .check_ont(ontology %||% parsed$ontology %||% "BP")
    if (is.null(dag)) dag <- go_dag(ont, obo = obo, verbose = FALSE)
    if (!identical(parsed$ontology, ont) && verbose) {
        warning("Using ontology '", ont, "' instead of detected '",
                parsed$ontology, "'.", call. = FALSE)
    }
    res <- parsed$result
    if (!nrow(res)) stop("No usable enrichment terms found.", call. = FALSE)
    missing <- setdiff(res$ID, dag$nodes$id)
    if (length(missing) == length(res$ID)) {
        stop("None of the imported terms exist in the ", ont,
             " graph. Choose the ontology of this result explicitly ",
             "(BP, MF or CC).", call. = FALSE)
    }
    if (length(missing)) {
        warning(length(missing), " imported term(s) are absent from the GO DAG: ",
                paste(utils::head(missing, 3), collapse = ", "), call. = FALSE)
    }
    counts <- parsed$counts
    out <- list(
        result = res,
        enrich = NULL,
        gene = parsed$input_ids %||% character(),
        universe = character(),
        ontology = ont,
        annotation = NULL,
        params = list(
            keytype = "external",
            source = parsed$source,
            source_file = basename(file),
            input_size = parsed$n_input,
            reference_size = parsed$n_ref,
            adjusted_source = parsed$adjusted_source %||% "source FDR",
            pvalue_kind = parsed$pvalue_kind %||% "source_adjusted",
            qvalue_is_alias = isTRUE(parsed$qvalue_is_alias),
            n_missing_adjusted = parsed$n_missing_adjusted %||% 0L,
            stat_default = if (all(is.na(res$p.adjust))) "pvalue" else "p.adjust",
            package_version = as.character(utils::packageVersion("PruneTopGODAGviz")),
            r_version = R.version.string,
            imported_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
        ),
        external_counts = counts
    )
    class(out) <- "go_enrich"
    if (verbose) {
        message(sprintf("Imported %d external terms from %s (%s)",
                        nrow(res), parsed$source, ont))
    }
    out
}

#' Parse a PANTHER over-representation JSON file
#'
#' @param file Path to the JSON file.
#' @return A parsed list with `result`, `counts`, `input_ids`, `n_input`,
#'   `n_ref`, `ontology`, `source` and `adjusted_source`.
#' @export
parse_panther_json <- function(file) {
    raw <- jsonlite::fromJSON(file, simplifyVector = FALSE)
    over <- raw$overrepresentation %||% raw
    groups <- over$group
    if (is.null(groups)) stop("No 'overrepresentation$group' found.", call. = FALSE)
    if (!is.list(groups)) groups <- list(groups)

    as_results <- function(x) {
        if (is.null(x)) return(list())
        if (is.list(x) && !is.null(x$term)) return(list(x))
        if (is.list(x)) return(x)
        list()
    }
    value <- function(x, key, default = NA) {
        v <- x[[key]]
        if (is.null(v) || !length(v)) default else v[[1]]
    }
    rows <- list()
    for (g in groups) {
        if (!is.list(g)) next
        for (r in as_results(g$result)) {
            if (!is.list(r) || is.null(r$term)) next
            term <- r$term
            input <- r$input_list %||% list()
            mids <- input$mapped_id_list$mapped_id %||% list()
            if (is.list(mids)) mids <- unlist(mids, use.names = FALSE)
            ids <- unique(trimws(as.character(mids)))
            ids <- ids[!is.na(ids) & nzchar(ids)]
            rows[[length(rows) + 1L]] <- data.frame(
                ID = as.character(value(term, "id", NA_character_)),
                Description = as.character(value(term, "label", NA_character_)),
                term_level = as.numeric(value(term, "level", NA_real_)),
                Count = as.numeric(value(input, "number_in_list",
                                         value(r, "number_in_list", NA_real_))),
                BgCount = as.numeric(value(r, "number_in_reference", NA_real_)),
                pvalue = as.numeric(value(input, "pValue",
                                          value(r, "pValue", NA_real_))),
                p.adjust = as.numeric(value(input, "fdr",
                                            value(r, "fdr", NA_real_))),
                FoldEnrichment = as.numeric(value(input, "fold_enrichment",
                                                  value(r, "fold_enrichment", NA_real_))),
                expected = as.numeric(value(input, "expected",
                                            value(r, "expected", NA_real_))),
                geneID = paste(ids, collapse = "/"),
                stringsAsFactors = FALSE
            )
        }
    }
    if (!length(rows)) stop("No PANTHER result rows found.", call. = FALSE)
    res <- do.call(rbind, rows)
    res <- res[nzchar(res$ID) & !is.na(res$ID), , drop = FALSE]
    if (anyDuplicated(res$ID)) {
        n_dup <- sum(duplicated(res$ID))
        warning(n_dup, " duplicate term identifier(s) ignored.", call. = FALSE)
        res <- res[!duplicated(res$ID), , drop = FALSE]
    }
    # The JSON occasionally stores the result as one object instead of an array;
    # the flattening above handles both shapes.
    input_ids <- unique(unlist(strsplit(res$geneID, "/", fixed = TRUE),
                               use.names = FALSE))
    input_ids <- input_ids[nzchar(input_ids)]
    n_input <- if (length(input_ids)) length(input_ids) else NA_real_
    n_ref <- .infer_reference_size(res, n_input)
    res$GeneRatio <- if (is.finite(n_input)) {
        paste0(format(res$Count, trim = TRUE, scientific = FALSE), "/",
               format(n_input, trim = TRUE, scientific = FALSE))
    } else NA_character_
    res$BgRatio <- if (is.finite(n_ref)) {
        paste0(format(res$BgCount, trim = TRUE, scientific = FALSE), "/",
               format(round(n_ref), trim = TRUE, scientific = FALSE))
    } else NA_character_
    res$qvalue <- res$p.adjust
    res$RichFactor <- ifelse(!is.na(res$BgCount) & res$BgCount > 0,
                             res$Count / res$BgCount, NA_real_)
    res$source <- "PANTHER"
    res$source_file <- basename(file)
    # Keep the count line deliberately raw: Count | BgCount.
    res$count_label <- ifelse(
        !is.na(res$Count) & !is.na(res$BgCount),
        paste0(format(res$Count, trim = TRUE, scientific = FALSE), " | ",
               format(res$BgCount, trim = TRUE, scientific = FALSE)),
        ""
    )
    counts <- data.frame(
        ID = res$ID,
        n_in_term = res$Count,
        n_bg_term = res$BgCount,
        n_in_universe = n_input,
        n_bg_universe = n_ref,
        count_label = res$count_label,
        stringsAsFactors = FALSE
    )
    ontology <- .ontology_from_text(
        over$annotation_type %||% over$annotationType %||% "")
    pvalue_kind <- if (all(is.na(res$p.adjust))) {
        "raw_only"
    } else if (anyNA(res$p.adjust)) {
        "source_adjusted_with_missing"
    } else {
        "source_adjusted"
    }
    list(result = res, counts = counts, input_ids = input_ids,
         n_input = n_input, n_ref = n_ref, ontology = ontology,
         source = "PANTHER", adjusted_source = "source FDR",
         pvalue_kind = pvalue_kind, qvalue_is_alias = TRUE,
         n_missing_adjusted = sum(is.na(res$p.adjust)))
}

# Numerator of "numerator/denominator" ratio cells.
#
# Returns the numerator. A cell is valid only when it holds two numbers, the
# denominator is positive and 0 <= numerator <= denominator; anything else
# becomes NA and is reported once, so a malformed ratio never silently turns
# into a wrong count.
.ratio_numerator <- function(x, field) {
    txt <- trimws(as.character(x))
    parts <- strsplit(txt, "/", fixed = TRUE)
    ok <- lengths(parts) == 2L
    num <- rep(NA_real_, length(txt))
    den <- rep(NA_real_, length(txt))
    if (any(ok)) {
        num[ok] <- suppressWarnings(as.numeric(vapply(parts[ok], `[[`, "", 1L)))
        den[ok] <- suppressWarnings(as.numeric(vapply(parts[ok], `[[`, "", 2L)))
    }
    good <- ok & is.finite(num) & is.finite(den) & num >= 0 & den > 0 &
        num <= den
    bad <- !good
    if (any(bad)) {
        warning(sum(bad), " invalid ", field,
                " value(s) ignored; expected 'numerator/denominator'.",
                call. = FALSE)
    }
    num[bad] <- NA_real_
    num
}

#' Parse a generic external enrichment table
#'
#' @param file CSV or TSV file.
#' @param ontology Optional ontology label.
#' @return A parsed list in the same shape as [parse_panther_json()].
#' @export
parse_external_table <- function(file, ontology = NULL) {
    ext <- tolower(tools::file_ext(file))
    tab <- if (ext %in% c("tsv", "txt")) {
        utils::read.delim(file, check.names = FALSE, stringsAsFactors = FALSE)
    } else {
        utils::read.csv(file, check.names = FALSE, stringsAsFactors = FALSE)
    }
    norm <- function(x) gsub("[^a-z0-9]", "", tolower(x))
    pick <- function(candidates) {
        nm <- norm(names(tab))
        hit <- match(norm(candidates), nm)
        hit <- hit[!is.na(hit)]
        if (!length(hit)) NULL else names(tab)[hit[[1]]]
    }
    id_col <- pick(c("ID", "TermID", "term_id", "go_id", "GO", "Term"))
    p_col <- pick(c("pvalue", "PValue", "p", "P"))
    if (is.null(id_col) || is.null(p_col)) {
        stop("Generic table needs an ID/term column and a p-value column.",
             call. = FALSE)
    }
    padj_col <- pick(c("p.adjust", "padj", "FDR", "fdr", "qvalue",
                       "q.value", "adjusted_pvalue"))
    desc_col <- pick(c("Description", "Term", "TermName", "name", "label"))
    count_col <- pick(c("Count", "number_in_list", "GeneCount", "n"))
    bg_col <- pick(c("BgCount", "number_in_reference", "background_count"))
    gene_col <- pick(c("geneID", "genes", "mapped_id", "mapped_ids"))
    fold_col <- pick(c("FoldEnrichment", "fold_enrichment", "fold"))
    gr_col <- pick(c("GeneRatio", "gene_ratio"))
    br_col <- pick(c("BgRatio", "bg_ratio"))
    ont_col <- pick(c("Namespace", "Ontology", "ont", "GO_namespace"))
    val <- function(col) if (is.null(col)) rep(NA, nrow(tab)) else tab[[col]]
    res <- data.frame(
        ID = trimws(as.character(val(id_col))),
        Description = if (is.null(desc_col)) as.character(val(id_col)) else
            as.character(val(desc_col)),
        pvalue = suppressWarnings(as.numeric(val(p_col))),
        p.adjust = suppressWarnings(as.numeric(val(padj_col))),
        Count = suppressWarnings(as.numeric(val(count_col))),
        BgCount = suppressWarnings(as.numeric(val(bg_col))),
        geneID = as.character(val(gene_col)),
        FoldEnrichment = suppressWarnings(as.numeric(val(fold_col))),
        stringsAsFactors = FALSE
    )
    if (is.null(padj_col)) {
        res$p.adjust <- rep(NA_real_, nrow(res))
        adjusted_source <- "raw p-value only"
    } else {
        adjusted_source <- "source FDR"
    }
    if (is.null(count_col) && !is.null(gr_col)) {
        res$Count <- .ratio_numerator(tab[[gr_col]], "GeneRatio")
    } else if (!is.null(count_col) && !is.null(gr_col)) {
        num <- .ratio_numerator(tab[[gr_col]], "GeneRatio")
        bad <- !is.na(num) & !is.na(res$Count) & num != res$Count
        if (any(bad)) {
            warning(sum(bad), " row(s) have Count that differs from the ",
                    "GeneRatio numerator; the Count column was kept.",
                    call. = FALSE)
        }
    }
    if (is.null(bg_col) && !is.null(br_col)) {
        res$BgCount <- .ratio_numerator(tab[[br_col]], "BgRatio")
    } else if (!is.null(bg_col) && !is.null(br_col)) {
        num <- .ratio_numerator(tab[[br_col]], "BgRatio")
        bad <- !is.na(num) & !is.na(res$BgCount) & num != res$BgCount
        if (any(bad)) {
            warning(sum(bad), " row(s) have BgCount that differs from the ",
                    "BgRatio numerator; the BgCount column was kept.",
                    call. = FALSE)
        }
    }
    if (is.null(gene_col)) res$geneID <- ""
    # A raw p-value must never impersonate an adjusted one: missing stays NA.
    res$qvalue <- res$p.adjust
    pvalue_kind <- if (!nrow(res) || all(is.na(res$p.adjust))) {
        "raw_only"
    } else if (anyNA(res$p.adjust)) {
        "source_adjusted_with_missing"
    } else {
        "source_adjusted"
    }
    res$RichFactor <- ifelse(!is.na(res$BgCount) & res$BgCount > 0,
                             res$Count / res$BgCount, NA_real_)
    res$source <- "external table"
    res$source_file <- basename(file)
    res$count_label <- ifelse(
        !is.na(res$Count) & !is.na(res$BgCount),
        paste0(format(res$Count, trim = TRUE, scientific = FALSE), " | ",
               format(res$BgCount, trim = TRUE, scientific = FALSE)), ""
    )
    res <- res[nzchar(res$ID) & !is.na(res$ID), , drop = FALSE]
    if (!nrow(res)) stop("No usable rows in external table.", call. = FALSE)
    input_ids <- unique(unlist(strsplit(res$geneID, "/", fixed = TRUE),
                               use.names = FALSE))
    input_ids <- input_ids[nzchar(input_ids)]
    n_input <- if (length(input_ids)) length(input_ids) else NA_real_
    n_ref <- NA_real_
    counts <- data.frame(
        ID = res$ID, n_in_term = res$Count, n_bg_term = res$BgCount,
        n_in_universe = n_input, n_bg_universe = n_ref,
        count_label = res$count_label, stringsAsFactors = FALSE
    )
    detected <- NULL
    if (!is.null(ont_col)) {
        vals <- trimws(as.character(tab[[ont_col]]))
        vals <- vals[nzchar(vals) & !is.na(vals)]
        if (length(vals)) {
            detected <- .ontology_from_text(names(sort(table(vals),
                                                       decreasing = TRUE))[1])
        }
    }
    list(result = res, counts = counts, input_ids = input_ids,
         n_input = n_input, n_ref = n_ref,
         ontology = .check_ont(ontology %||% detected %||% "BP"),
         source = "external table", adjusted_source = adjusted_source,
         pvalue_kind = pvalue_kind, qvalue_is_alias = TRUE,
         n_missing_adjusted = sum(is.na(res$p.adjust)))
}

.infer_reference_size <- function(res, n_input) {
    if (!is.finite(n_input) || !n_input) return(NA_real_)
    ok <- !is.na(res$Count) & !is.na(res$BgCount) & !is.na(res$expected) &
        res$Count > 0 & res$expected > 0
    if (!any(ok)) return(NA_real_)
    vals <- n_input * res$BgCount[ok] / res$expected[ok]
    vals <- vals[is.finite(vals) & vals > 0]
    if (!length(vals)) return(NA_real_)
    med <- stats::median(vals)
    if (stats::median(abs(vals - med)) > med * 1e-5) return(NA_real_)
    med
}

.ontology_from_text <- function(text) {
    text <- tolower(trimws(as.character(text)))
    text <- gsub("[_-]+", " ", text)
    if (!length(text) || !nzchar(text)) return("BP")
    if (text %in% c("bp", "biological process")) return("BP")
    if (text %in% c("mf", "molecular function")) return("MF")
    if (text %in% c("cc", "cellular component")) return("CC")
    if (grepl("biological process", text, fixed = TRUE)) return("BP")
    if (grepl("molecular function", text, fixed = TRUE)) return("MF")
    if (grepl("cellular component", text, fixed = TRUE)) return("CC")
    "BP"
}
