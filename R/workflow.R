#' Run the full GO workflow for one or more ontologies
#'
#' Convenience wrapper that runs enrichment, pruning, layout and export for each
#' ontology, writing one figure per ontology plus result and decision tables.
#'
#' @param genes Input identifiers.
#' @param annotation A `go_annotation` object.
#' @param outdir Output directory.
#' @param ont Ontologies to run.
#' @param prefix File name prefix.
#' @param prune Named list of arguments for [prune_go_dag()].
#' @param plot Named list of arguments for [plot_go_dag()].
#' @param interactive Also write an interactive HTML file per ontology.
#' @param obo Optional path to a local OBO file.
#' @param overwrite Overwrite existing figure files.
#' @param verbose Print progress messages.
#'
#' @return A named list (one element per ontology) with `enrich`, `dag` and
#'   `files`, invisibly.
#' @export
run_go_workflow <- function(genes, annotation, outdir, ont = c("BP", "MF", "CC"),
                            prefix = "GO", prune = list(), plot = list(),
                            interactive = FALSE, obo = NULL, overwrite = TRUE,
                            verbose = TRUE) {
    if (!inherits(annotation, "go_annotation")) {
        stop("'annotation' must be a go_annotation object.", call. = FALSE)
    }
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
    out <- list()
    for (o in ont) {
        o <- .check_ont(o)
        if (verbose) message("== ", o, " ==")
        ego <- enrich_go(genes, annotation, ont = o, verbose = verbose)
        prune_args <- c(list(x = ego, dag = go_dag(o, obo = obo,
                                                            verbose = verbose)),
                            prune)
        if (is.null(prune_args$verbose)) prune_args$verbose <- verbose
        dag <- do.call(prune_go_dag, prune_args)
        p <- do.call(plot_go_dag, c(list(x = dag), plot))
        files <- save_go_dag(dag, file.path(outdir, paste0(prefix, "_", o)),
                             plot = p, overwrite = overwrite)
        utils::write.csv(ego$result,
                         file.path(outdir, paste0(prefix, "_", o, "_enrichment.csv")),
                         row.names = FALSE)
        write_prune_table(dag,
                          file.path(outdir, paste0(prefix, "_", o, "_prune_table.csv")))
        if (isTRUE(interactive)) {
            html <- file.path(outdir, paste0(prefix, "_", o, ".html"))
            plot_go_dag_interactive(dag, file = html)
            files <- c(files, html)
        }
        out[[o]] <- list(enrich = ego, dag = dag, files = files)
    }
    si <- file.path(outdir, paste0(prefix, "_sessionInfo.txt"))
    writeLines(utils::capture.output(utils::sessionInfo()), si)
    invisible(out)
}
