#' Read a transcript-level GO annotation table
#'
#' Reads a tabular file in which every row maps one transcript (or gene) to one
#' GO term, and converts it into a compact annotation object that can be used by
#' [enrich_go()].  The function is species agnostic: the GO hierarchy is not
#' required here, only the identifier column, the GO identifier and (optionally)
#' the term name, ontology and definition.
#'
#' @param file Path to a CSV or TSV file.
#' @param col_transcript,col_go,col_name,col_ontology,col_def Column names.
#'   `NULL` means auto-detection by regular expression.
#' @param sep Field separator, or `"auto"` to guess from the first line.
#' @param ont Restrict the annotation to one ontology (`"BP"`, `"MF"`, `"CC"`).
#'   `NULL` keeps all ontologies present in the file.
#' @param keytype Use `"transcript"` or `"gene"` identifiers.
#' @param include_def Read the definition column as well.  Definitions are long
#'   and are only needed for interactive tooltips, so they are skipped by
#'   default.
#' @param verbose Print progress messages.
#'
#' @return An object of class `go_annotation`, a list with elements `pairs`
#'   (unique transcript--GO pairs), `terms` (one row per GO term), `universe`
#'   (annotated identifiers per ontology), `keytype` and `source`.
#'
#' @examples
#' \dontrun{
#' anno <- read_go_annotation("At_GO_info.csv", ont = "BP")
#' head(anno$terms)
#' }
#' @export
read_go_annotation <- function(file, col_transcript = NULL, col_go = NULL,
                               col_name = NULL, col_ontology = NULL,
                               col_def = NULL, sep = "auto", ont = NULL, keytype = c("transcript", "gene"),
                               include_def = FALSE, verbose = TRUE) {
    if (!file.exists(file)) {
        stop("Annotation file not found: ", file, call. = FALSE)
    }
    if (identical(sep, "auto")) sep <- .guess_sep(file)
    header <- names(data.table::fread(file, nrows = 0L, sep = sep,
                                      showProgress = FALSE, data.table = TRUE))
    first_line <- .read_first_line(file)
    manual_header <- scan(text = first_line, what = "", sep = sep,
                          quiet = TRUE, quote = "\"")
    manual_header[manual_header == ""] <- "row_id"
    manual <- all(grepl("^V[0-9]+$", header)) &&
        length(manual_header) >= length(header)
    if (manual) header <- manual_header
    header_clean <- trimws(header)

    pick <- function(spec, patterns, what) {
        if (!is.null(spec)) {
            if (!spec %in% header) {
                stop("Column '", spec, "' (", what, ") not found. Available: ",
                     paste(header, collapse = ", "), call. = FALSE)
            }
            return(spec)
        }
        for (pat in patterns) {
            hit <- grep(pat, tolower(header_clean), value = FALSE, perl = TRUE)
            if (length(hit)) return(header[hit[[1]]])
        }
        stop("Could not detect the ", what, " column. Supply it explicitly. Available: ",
             paste(header, collapse = ", "), call. = FALSE)
    }

    keytype <- match.arg(keytype)
    key_patterns <- if (keytype == "gene") {
        c("gene_id$", "ensembl_gene", "^gene$")
    } else {
        c("transcript", "tgid", "mrna")
    }
    col_transcript <- pick(col_transcript, key_patterns, keytype)
    col_go <- pick(col_go, c("^go_?id$", "^go$", "goid"), "GO id")
    col_name <- pick(col_name, c("name_1006", "^term_?name$", "^term$", "^name$"), "term name")
    col_ontology <- pick(col_ontology, c("namespace", "ontology", "^ont$"), "ontology")
    if (!is.null(col_def) || isTRUE(include_def)) {
        col_def <- pick(col_def, c("^description$", "^def$"), "definition")
    } else {
        col_def <- NULL
    }

    select <- unique(c(col_transcript, col_go, col_name, col_ontology, col_def))
    if (manual) {
        idx <- match(select, manual_header)
        dt <- data.table::fread(file, sep = sep, skip = 1L, header = FALSE,
                                select = idx, col.names = manual_header[idx],
                                colClasses = "character", showProgress = FALSE,
                                na.strings = c("", "NA", "N/A", "na"))
    } else {
        dt <- data.table::fread(file, sep = sep, select = select,
                                colClasses = "character", showProgress = FALSE,
                                na.strings = c("", "NA", "N/A", "na"))
    }

    tr <- trimws(dt[[col_transcript]])
    go <- toupper(trimws(dt[[col_go]]))
    nm <- trimws(dt[[col_name]])
    on <- .norm_ontology(dt[[col_ontology]])
    df <- if (is.null(col_def)) {
        rep(NA_character_, length(tr))
    } else {
        trimws(dt[[col_def]])
    }

    ok <- !is.na(tr) & nzchar(tr) & !is.na(go) & grepl("^GO:[0-9]{7}$", go)
    n_bad <- sum(!ok)
    tr <- tr[ok]; go <- go[ok]; nm <- nm[ok]; on <- on[ok]; df <- df[ok]

    if (!length(go)) {
        stop("No valid 'GO:0000000'-style identifiers found in column '", col_go, "'.",
             call. = FALSE)
    }
    if (all(is.na(on))) {
        if (is.null(ont)) {
            stop("No ontology information found; supply 'ont' or a namespace column.",
                 call. = FALSE)
        }
        on[] <- .check_ont(ont)
    }
    if (anyNA(on)) {
        bad <- unique(go[is.na(on)])[1:min(3, length(unique(go[is.na(on)])))]
        warning("Dropping rows with unrecognised ontology (e.g. ",
                paste(bad, collapse = ", "), ").", call. = FALSE)
        keep <- !is.na(on)
        tr <- tr[keep]; go <- go[keep]; nm <- nm[keep]; on <- on[keep]; df <- df[keep]
    }

    dfp <- data.frame(transcript = tr, go_id = go, name = nm,
                      ontology = on, def = df, stringsAsFactors = FALSE)
    dfp <- unique(dfp[, c("transcript", "go_id", "name", "ontology", "def")])

    ord <- order(dfp$go_id, is.na(dfp$name))
    terms <- dfp[ord, c("go_id", "name", "ontology", "def")]
    terms <- terms[!duplicated(terms$go_id), , drop = FALSE]
    terms$name[is.na(terms$name) | !nzchar(terms$name)] <-
        terms$go_id[is.na(terms$name) | !nzchar(terms$name)]
    rownames(terms) <- NULL

    pairs <- unique(dfp[, c("transcript", "go_id")])
    rownames(pairs) <- NULL

    if (!is.null(ont)) {
        ont <- .check_ont(ont)
        keep_go <- terms$go_id[terms$ontology == ont]
        pairs <- pairs[pairs$go_id %in% keep_go, , drop = FALSE]
        terms <- terms[terms$go_id %in% keep_go, , drop = FALSE]
    }

    ont_levels <- c("BP", "MF", "CC")
    universe <- lapply(ont_levels, function(o) {
        g <- terms$go_id[terms$ontology == o]
        sort(unique(pairs$transcript[pairs$go_id %in% g]))
    })
    names(universe) <- ont_levels

    out <- list(pairs = pairs, terms = terms, universe = universe,
                keytype = keytype,
                source = list(file = normalizePath(file, mustWork = FALSE),
                              size = file.size(file),
                              mtime = as.numeric(file.mtime(file)),
                              n_rows = nrow(pairs),
                              n_bad = n_bad,
                              sep = sep))
    class(out) <- "go_annotation"
    if (verbose) {
        message(sprintf("GO annotation: %d pairs, %d terms (BP %d / MF %d / CC %d)",
                        nrow(pairs), nrow(terms),
                        sum(terms$ontology == "BP"), sum(terms$ontology == "MF"),
                        sum(terms$ontology == "CC")))
    }
    out
}

.read_first_line <- function(file) {
    con <- if (grepl("\\.gz$", file, ignore.case = TRUE)) {
        gzfile(file, "rt")
    } else {
        file(file, "rt")
    }
    on.exit(close(con), add = TRUE)
    readLines(con, n = 1L, warn = FALSE)
}

.guess_sep <- function(file) {
    l <- .read_first_line(file)
    n_comma <- lengths(regmatches(l, gregexpr(",", l, fixed = TRUE)))
    n_tab <- lengths(regmatches(l, gregexpr("\t", l, fixed = TRUE)))
    if (n_tab > n_comma) "\t" else ","
}

#' @export
print.go_annotation <- function(x, ...) {
    cat("<go_annotation>\n")
    cat("  pairs :", nrow(x$pairs), "\n")
    cat("  terms :", nrow(x$terms), "\n")
    for (o in c("BP", "MF", "CC")) {
        cat(sprintf("  %s    : %d transcripts\n", o, length(x$universe[[o]])))
    }
    cat("  source:", basename(x$source$file), "\n")
    invisible(x)
}

#' @export
as.data.frame.go_annotation <- function(x, ...) {
    x$pairs
}
