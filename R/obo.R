.pkg_cache <- new.env(parent = emptyenv())

#' Download the Gene Ontology basic OBO file
#'
#' The file is cached in the user cache directory and downloaded only once.
#' Proxy settings are taken from `proxy` or from the `https_proxy` /
#' `http_proxy` environment variables.
#'
#' @param dest Destination path.  Defaults to the package cache directory.
#' @param url Source URL.
#' @param proxy Optional proxy URL, e.g. `"http://127.0.0.1:7890"`.
#' @param overwrite Re-download even if a cached file exists.
#' @param verbose Print progress messages.
#'
#' @return The path to the OBO file, invisibly.
#' @export
download_go_obo <- function(dest = NULL,
                            url = "https://purl.obolibrary.org/obo/go/go-basic.obo",
                            proxy = NULL, overwrite = FALSE, verbose = TRUE) {
    if (is.null(dest)) {
        dest <- file.path(tools::R_user_dir("PruneTopGODAGviz", "cache"),
                          "go-basic.obo")
    }
    dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
    if (file.exists(dest) && !overwrite && file.size(dest) > 1e6) {
        return(invisible(dest))
    }
    old <- Sys.getenv(c("http_proxy", "https_proxy"), unset = NA_character_)
    if (!is.null(proxy)) {
        Sys.setenv(http_proxy = proxy, https_proxy = proxy)
    }
    on.exit({
        for (nm in names(old)) {
            if (is.na(old[[nm]])) Sys.unsetenv(nm) else
                do.call(Sys.setenv, stats::setNames(list(old[[nm]]), nm))
        }
    }, add = TRUE)
    if (verbose) message("Downloading GO OBO from ", url)
    ok <- tryCatch({
        utils::download.file(url, dest, method = "libcurl", quiet = !verbose,
                             mode = "wb")
        TRUE
    }, error = function(e) {
        if (verbose) message("Download failed: ", conditionMessage(e))
        FALSE
    })
    if (!ok || !file.exists(dest) || file.size(dest) < 1e6) {
        unlink(dest)
        stop("Could not download the GO OBO file. Supply a local 'obo' path ",
             "or set 'source = \"GO.db\"'.", call. = FALSE)
    }
    invisible(dest)
}

#' Parse a Gene Ontology OBO file
#'
#' A minimal, fast parser for `go-basic.obo` that keeps the fields needed for
#' plotting: identifiers, names, namespace, `is_a` and `part_of` relations,
#' obsolete flags and alternate identifiers.
#'
#' @param file Path to an OBO file.
#' @param ont Restrict to one ontology (`"BP"`, `"MF"`, `"CC"`); `NULL` keeps all.
#' @param include_obsolete Keep obsolete terms.
#' @param verbose Print progress messages.
#'
#' @return A `go_dag` object with `nodes`, `edges`, `roots`, `source` and
#'   `version`.
#' @export
read_obo <- function(file, ont = NULL, include_obsolete = FALSE, verbose = TRUE) {
    if (!file.exists(file)) stop("OBO file not found: ", file, call. = FALSE)
    if (verbose) message("Parsing OBO: ", file)
    lines <- readLines(file, warn = FALSE)
    hdr <- grep("^\\[[A-Za-z_]+\\]$", lines)
    if (!length(hdr)) stop("Not an OBO file: no stanza headers found.", call. = FALSE)
    header <- lines[hdr]
    pos <- findInterval(seq_along(lines), hdr)
    type <- rep(NA_character_, length(lines))
    has <- pos > 0L
    type[has] <- header[pos[has]]
    is_term <- !is.na(type) & type == "[Term]"

    ln <- lines[is_term]
    sid <- pos[is_term]
    field <- sub(":.*$", "", ln)
    value <- trimws(sub("^[^:]+:\\s*", "", ln))
    d <- data.table::data.table(sid = sid, field = field, value = value)

    sids <- sort(unique(d$sid))
    get1 <- function(f) {
        s <- d[d$field == f, c("sid", "value"), with = FALSE]
        if (!nrow(s)) return(rep(NA_character_, length(sids)))
        s <- s[!duplicated(s$sid)]
        s$value[match(sids, s$sid)]
    }
    terms <- data.table::data.table(
        sid = sids,
        id = get1("id"),
        name = get1("name"),
        namespace = .norm_ontology(get1("namespace")),
        definition = get1("def"),
        obsolete = tolower(get1("is_obsolete")) %in% c("true", "1", "yes")
    )
    terms <- terms[!is.na(terms$id) & nzchar(terms$id), ]

    isa <- d[d$field == "is_a", c("sid", "value"), with = FALSE]
    isa[, parent := sub("\\s*!.*$", "", value)]
    isa[, type := "is_a"]
    rel <- d[d$field == "relationship" & grepl("^part_of\\s", value),
             c("sid", "value"), with = FALSE]
    rel[, parent := sub("^part_of\\s+", "", value)]
    rel[, parent := sub("\\s*!.*$", "", parent)]
    rel[, type := "part_of"]
    ed <- data.table::rbindlist(list(isa[, c("sid", "parent", "type")],
                                     rel[, c("sid", "parent", "type")]))
    if (nrow(ed)) {
        ed <- merge(ed, terms[, c("sid", "id"), with = FALSE], by = "sid", all.x = FALSE)
        ed <- ed[, c("parent", "id", "type"), with = FALSE]
        data.table::setnames(ed, c("parent", "child", "type"))
    } else {
        ed <- data.table::data.table(parent = character(), child = character(),
                                     type = character())
    }

    alt <- d[d$field == "alt_id", c("sid", "value"), with = FALSE]
    if (nrow(alt)) {
        alt <- merge(alt, terms[, c("sid", "id"), with = FALSE], by = "sid")
        remap <- function(x) {
            m <- alt$id[match(x, alt$value)]
            ifelse(is.na(m) | !nzchar(m), x, m)
        }
        if (nrow(ed)) {
            ed[, parent := remap(parent)]
            ed[, child := remap(child)]
        }
    }

    if (!include_obsolete) {
        keep_id <- terms$id[!terms$obsolete]
        terms <- terms[terms$id %in% keep_id, ]
        ed <- ed[ed$parent %in% keep_id & ed$child %in% keep_id, ]
    }
    if (!is.null(ont)) {
        ont <- .check_ont(ont)
        keep_id <- terms$id[!is.na(terms$namespace) & terms$namespace == ont]
        terms <- terms[terms$id %in% keep_id, ]
        ed <- ed[ed$parent %in% keep_id & ed$child %in% keep_id, ]
    }
    ed <- as.data.frame(unique(ed))
    terms <- as.data.frame(unique(terms))

    nodes <- data.frame(
        id = terms$id, name = terms$name,
        namespace = terms$namespace, definition = terms$definition,
        obsolete = terms$obsolete, stringsAsFactors = FALSE
    )
    nodes$name[is.na(nodes$name) | !nzchar(nodes$name)] <- nodes$id[is.na(nodes$name) | !nzchar(nodes$name)]

    roots <- setdiff(nodes$id, ed$child)
    version <- grep("^(data-version|format-version|date):", lines[seq_len(min(50, length(lines)))],
                    value = TRUE)
    out <- list(nodes = nodes, edges = ed, roots = roots, source = "obo",
                version = paste(version, collapse = "; "),
                file = normalizePath(file, mustWork = FALSE))
    class(out) <- "go_dag"
    if (verbose) {
        message(sprintf("OBO graph: %d terms, %d edges, %d roots",
                        nrow(nodes), nrow(ed), length(roots)))
    }
    out
}

#' Build a GO directed acyclic graph
#'
#' @param ont One of `"BP"`, `"MF"` or `"CC"`.
#' @param obo Optional path to a local `go-basic.obo`.  When `NULL` the value
#'   of `getOption("PruneTopGODAGviz.obo")` is used; otherwise the cached file
#'   is used, and if no cache exists the file is downloaded.
#' @param source `"auto"` (default), `"obo"` or `"GO.db"`.  The `"GO.db"`
#'   fallback does not know relation types, so all edges are labelled `is_a`.
#' @param cache_dir Cache directory for the OBO file.
#' @param verbose Print progress messages.
#'
#' @return A `go_dag` object.
#' @export
go_dag <- function(ont = "BP", obo = NULL,
                   source = c("auto", "obo", "GO.db"), cache_dir = NULL,
                   verbose = TRUE) {
    ont <- .check_ont(ont)
    source <- match.arg(source)
    if (is.null(obo)) obo <- getOption("PruneTopGODAGviz.obo", NULL)
    if (!is.null(obo) && !file.exists(obo)) {
        warning("Ignoring option 'PruneTopGODAGviz.obo': file not found: ",
                obo, call. = FALSE)
        obo <- NULL
    }
    if (source == "auto") source <- "obo"

    if (source == "obo") {
        path <- obo
        if (is.null(path)) {
            path <- file.path(cache_dir %||%
                                  file.path(tools::R_user_dir("PruneTopGODAGviz", "cache")),
                              "go-basic.obo")
            if (!file.exists(path)) {
                path <- tryCatch(download_go_obo(dest = path, verbose = verbose),
                                 error = function(e) NULL)
            }
        }
        if (is.null(path) || !file.exists(path)) {
            warning("GO OBO file unavailable; falling back to GO.db. ",
                    "Relation types will all be labelled 'is_a'.",
                    call. = FALSE)
            return(.go_dag_from_godb(ont, verbose = verbose))
        }
        key <- paste("obo", ont, normalizePath(path, mustWork = FALSE),
                     file.mtime(path), sep = "|")
        if (!is.null(.pkg_cache[[key]])) return(.pkg_cache[[key]])
        dag <- read_obo(path, ont = ont, verbose = verbose)
        .pkg_cache[[key]] <- dag
        return(dag)
    }

    .go_dag_from_godb(ont, verbose = verbose)
}

.go_dag_from_godb <- function(ont, verbose = TRUE) {
    .check_pkg("GO.db", "the GO.db fallback")
    key <- paste("godb", ont, sep = "|")
    if (!is.null(.pkg_cache[[key]])) return(.pkg_cache[[key]])
    kids <- get(paste0("GO", ont, "CHILDREN"), envir = asNamespace("GO.db"))
    lst <- as.list(kids)
    lens <- lengths(lst)
    ed <- data.frame(parent = rep(names(lst), lens),
                     child = unlist(lst, use.names = FALSE),
                     type = "is_a", stringsAsFactors = FALSE)
    ids <- unique(c(ed$parent, ed$child))
    ann <- suppressMessages(AnnotationDbi::select(
        GO.db::GO.db, keys = ids, columns = c("TERM", "ONTOLOGY"),
        keytype = "GOID"))
    nodes <- data.frame(
        id = ann$GOID, name = ann$TERM, namespace = .norm_ontology(ann$ONTOLOGY),
        definition = NA_character_, obsolete = FALSE, stringsAsFactors = FALSE)
    nodes <- unique(nodes)
    roots <- setdiff(nodes$id, ed$child)
    out <- list(nodes = nodes, edges = ed, roots = roots, source = "GO.db",
                version = paste("GO.db", utils::packageVersion("GO.db")),
                file = NA_character_)
    class(out) <- "go_dag"
    if (verbose) {
        message(sprintf("GO.db graph: %d terms, %d edges (edge type unknown)",
                        nrow(nodes), nrow(ed)))
    }
    .pkg_cache[[key]] <- out
    out
}

#' @export
print.go_dag <- function(x, ...) {
    cat("<go_dag> source:", x$source, "\n")
    cat("  terms :", nrow(x$nodes), "\n")
    cat("  edges :", nrow(x$edges), "\n")
    if (nrow(x$edges)) {
        cat("  types :", paste(names(table(x$edges$type)),
                               table(x$edges$type), sep = "=", collapse = ", "), "\n")
    }
    cat("  roots :", paste(utils::head(x$roots, 3), collapse = ", "), "\n")
    invisible(x)
}
