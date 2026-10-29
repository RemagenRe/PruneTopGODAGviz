.max_upload_bytes <- 2 * 1024^3 - 1

.builtin_species <- c(
    arabidopsis_thaliana = "At_GO_info_panther.csv",
    zea_mays = "Maize_GO_info_panther.csv",
    oryza_sativa = "Rice_GO_info_panther.csv",
    glycine_max = "Soybean_GO_info_panther.csv",
    homo_sapiens = "Homo_sapiens_GO_info_panther.csv",
    mus_musculus = "Mus_musculus_GO_info_panther.csv"
)

.builtin_dirs <- function() {
    d <- getOption("PruneTopGODAGviz.data", NULL)
    unique(c(d,
             file.path(getwd(), "data"),
             file.path(getwd(), "..", "data"),
             system.file("extdata", package = "PruneTopGODAGviz")))
}

.builtin_file <- function(name) {
    dirs <- .builtin_dirs()
    for (d in dirs) {
        if (!is.null(d) && nzchar(d) && dir.exists(d)) {
            for (candidate in c(name, paste0(name, ".gz"))) {
                p <- file.path(d, candidate)
                if (file.exists(p)) return(p)
            }
        }
    }
    file.path(dirs[[1]] %||% file.path(getwd(), "data"), name)
}

.builtin_available <- function() {
    names(.builtin_species)[vapply(.builtin_species, function(f) {
        file.exists(.builtin_file(f))
    }, logical(1))]
}

.builtin_download_ids <- c(
    arabidopsis_thaliana = "download_arabidopsis",
    zea_mays = "download_zea_mays",
    oryza_sativa = "download_rice",
    glycine_max = "download_soybean",
    homo_sapiens = "download_human",
    mus_musculus = "download_mouse"
)

.set_upload_limit <- function() {
    options(shiny.maxRequestSize = .max_upload_bytes)
    invisible(.max_upload_bytes)
}

.id_check_error <- function(message) {
    stop(structure(list(message = message),
                   class = c("id_check_error", "error", "condition")))
}

.check_id_match <- function(ids, universe, min_ids = 10L,
                            min_matched = 10L, min_rate = 0.2) {
    ids <- trimws(as.character(ids))
    ids <- ids[!is.na(ids) & nzchar(ids)]
    ids <- sort(unique(ids))
    n <- length(ids)
    if (n < min_ids) {
        .id_check_error(sprintf(
            "Only %d unique, non-empty IDs were found; at least %d are required.",
            n, min_ids))
    }
    universe <- unique(as.character(universe))
    universe <- universe[!is.na(universe) & nzchar(universe)]
    matched <- intersect(ids, universe)
    n_matched <- length(matched)
    rate <- if (n) n_matched / n else 0
    if (!length(universe)) {
        .id_check_error("The selected GO annotation has no IDs for this ontology.")
    }
    if (n_matched < min_matched) {
        .id_check_error(sprintf(
            "Only %d of %d IDs match the selected GO annotation; at least %d matches are required.",
            n_matched, n, min_matched))
    }
    if (rate < min_rate) {
        .id_check_error(sprintf(
            "Only %.1f%% of IDs match the selected GO annotation; at least %.0f%% are required.",
            100 * rate, 100 * min_rate))
    }
    list(ids = ids, matched = matched, n = n, n_matched = n_matched,
         rate = rate)
}

#' Launch the interactive GO DAG editor
#'
#' Opens a Shiny application in which enriched GO terms can be explored,
#' rearranged by dragging (with an optional force-directed layout), pruned by
#' deleting a term together with its descendants, recoloured with arbitrary
#' hexadecimal colours, and exported as PDF, PNG, SVG or TIFF using the current
#' arrangement.
#'
#' @param port Port for the Shiny server; `NULL` picks a free port.
#' @param launch.browser Open the default browser.
#' @param obo Optional path to a local go-basic.obo used by the example and
#'   by [go_dag()] when no explicit obo is supplied.
#' @param ... Passed to [shiny::runApp()].
#'
#' @return Called for its side effect; returns the value of [shiny::runApp()].
#' @export
run_app <- function(port = NULL, launch.browser = TRUE, obo = NULL, ...) {
    .check_pkg("shiny", "the interactive editor")
    .check_pkg("visNetwork", "the interactive editor")
    .check_pkg("colourpicker", "the colour controls")
    .check_pkg("future", "background enrichment")
    .check_pkg("promises", "background enrichment")
    old_libs <- Sys.getenv("R_LIBS")
    libs <- unique(c(.libPaths(),
                     strsplit(old_libs, .Platform$path.sep)[[1]]))
    Sys.setenv(R_LIBS = paste(libs[nzchar(libs)],
                              collapse = .Platform$path.sep))
    future::plan(future::multisession)
    .set_upload_limit()
    if (!is.null(obo)) options(PruneTopGODAGviz.obo = obo)
    if (is.null(getOption("PruneTopGODAGviz.data"))) {
        d <- file.path(getwd(), "data")
        if (dir.exists(d)) {
            options(PruneTopGODAGviz.data = normalizePath(d))
        }
    }
    appdir <- system.file("app", package = "PruneTopGODAGviz")
    if (!nzchar(appdir)) {
        appdir <- file.path(getwd(), "PruneTopGODAGviz", "inst", "app")
    }
    if (dir.exists(appdir)) {
        try(shiny::addResourcePath("appwww", appdir), silent = TRUE)
    }
    shiny::runApp(shiny::shinyApp(ui = app_ui(), server = app_server),
                  port = port, launch.browser = launch.browser, ...)
}

app_ui <- function() {
    avail <- .builtin_available()
    shiny::fluidPage(
        shiny::tags$head(shiny::tags$style(shiny::HTML(
            "body{background:#f4f5f7;font-family:Arial, Helvetica, sans-serif;}
             .well{padding:10px 12px;margin-bottom:10px;background:#fff;border:1px solid #e2e4e8;border-radius:6px;}
             .well h5{margin:2px 0 8px 0;font-weight:600;color:#2c3e50;font-size:12px;text-transform:uppercase;letter-spacing:.05em;}
             .form-group{margin-bottom:7px;}
             .control-label{font-size:12px;margin-bottom:2px;color:#333;}
             .form-control,input[type=number]{height:30px;padding:2px 8px;font-size:12px;}
             select.form-control{height:30px;padding:2px 8px;font-size:12px;}
             .btn{font-size:12px;padding:5px 10px;}
             .graph-wrap{position:relative;background:#fff;border:1px solid #d9dde3;border-radius:6px;overflow:hidden;}
             #graph{background:#fff;}
             .welcome{position:absolute;top:36%;left:0;right:0;text-align:center;color:#7a8797;font-size:15px;line-height:1.9;}
             .legend-note{font-size:11px;color:#6b7785;line-height:1.5;margin-top:2px;}
             .shiny-input-container .progress{height:0;margin-bottom:0;overflow:hidden;}
             .cbox{margin-top:2px;margin-bottom:2px;font-size:12px;}
             .cbox label{font-size:12px;}
             h2{font-size:20px;margin:8px 0 10px 0;font-weight:600;color:#1f2d3d;}
             .right-panel{position:sticky;top:10px;}
             .sel-box{font-size:12px;color:#333;background:#fafbfc;border:1px solid #e2e4e8;border-radius:4px;padding:6px;min-height:58px;white-space:pre-wrap;}
             .status{font-size:12px;color:#56606b;white-space:pre-wrap;}
             .colourpicker-input{height:30px !important;font-size:12px;}
             .anno-file-input{position:absolute;left:-9999px;width:1px;height:1px;overflow:hidden;}
            ")),
            shiny::tags$script(shiny::HTML(
                "Shiny.addCustomMessageHandler('anno_picker_delay', function(msg) {",
                "  if (window.__annoPickerTimer) clearTimeout(window.__annoPickerTimer);",
                "  window.__annoPickerTimer = setTimeout(function() {",
                "    window.__annoPickerTimer = null;",
                "    var el = document.getElementById('anno_file');",
                "    if (el) el.click();",
                "  }, (msg && msg.delay) ? msg.delay : 6000);",
                "});",
                "$(document).on('change', '#anno_file', function() {",
                "  if (window.__annoPickerTimer) {",
                "    clearTimeout(window.__annoPickerTimer);",
                "    window.__annoPickerTimer = null;",
                "  }",
                "});",
                "$(document).on('hidden.bs.modal', '.modal', function() {",
                "  if (window.__annoPickerTimer) {",
                "    clearTimeout(window.__annoPickerTimer);",
                "    window.__annoPickerTimer = null;",
                "  }",
                "});",
                "function __ptgTryReconnect() {",
                "  if (!document.getElementById('shiny-disconnected-overlay')) return;",
                "  if (window.Shiny && Shiny.shinyapp && !Shiny.shinyapp.isConnected()) {",
                "    try { Shiny.shinyapp.reconnect(); } catch(e) {}",
                "  }",
                "}",
                "document.addEventListener('visibilitychange', function() {",
                "  if (!document.hidden) setTimeout(__ptgTryReconnect, 800);",
                "});",
                "$(document).on('shiny:disconnected', function() {",
                "  setTimeout(__ptgTryReconnect, 2500);",
                "});"
            ))),
        shiny::titlePanel("PruneTopGODAGviz Studio"),
        shiny::fluidRow(
            shiny::column(3,
                shiny::div(class = "well",
                    shiny::h5("Data"),
                    shiny::actionButton("load_example", "Load Example Data",
                                        class = "btn-default", width = "100%"),
                    shiny::div(class = "anno-file-input",
                        shiny::fileInput("anno_file", "GO Annotation File",
                                         accept = c(".csv", ".tsv", ".txt"))),
                    shiny::actionButton("browse_anno", "Browse...", width = "100%"),
                    shiny::div(class = "legend-note", shiny::textOutput("anno_file_name")),
                    shiny::div(class = "legend-note",
                        "One transcript-GO pair per row: transcript/gene ID, GO ID, term name, namespace. Column names are auto-detected."),
                    shiny::selectInput("builtin_anno", "Or Use a Built-in Annotation",
                        c("None (Upload My Own)" = "",
                          stats::setNames(avail, avail))),
                    if (!length(avail)) shiny::div(class = "legend-note",
                        "Built-in annotation data directory not found; upload an annotation file or set PruneTopGODAGviz.data."),
                    shiny::fileInput("gene_file", "Gene / Transcript List",
                                     accept = c(".txt", ".csv", ".tsv")),
                    shiny::div(class = "legend-note",
                        "One ID per line. IDs are deduplicated and sorted; at least 10 IDs, 10 matches and 20% match rate are required."),
                    shiny::checkboxInput("propagate",
                        "Propagate Annotations to Ancestors", TRUE),
                    shiny::selectInput("keytype", "ID Type",
                                       c("Transcript" = "transcript",
                                         "Gene" = "gene")),
                    shiny::selectInput("ont", "Ontology",
                                       c("BP" = "BP", "MF" = "MF", "CC" = "CC")),
                    shiny::fluidRow(
                        shiny::column(6, shiny::numericInput("min_gs", "MinGSSize", 10, 1, 1000, 1)),
                        shiny::column(6, shiny::numericInput("max_gs", "MaxGSSize", 5000, 10, 50000, 10))
                    ),
                    shiny::actionButton("run", "Run Enrichment",
                                        class = "btn-primary btn-lg",
                                        style = "width:100%;font-size:15px;padding:9px;background:#2878B5;border-color:#2878B5;")
                ),
                shiny::div(class = "well",
                    shiny::h5("Display and Pruning"),
                    shiny::fluidRow(
                        shiny::column(6, shiny::numericInput("top_n", "Top Significant Terms", 5, 1, 200, 1)),
                        shiny::column(6, shiny::numericInput("p_cutoff", "P / FDR Cutoff", 0.05, 0, 1, 0.01))
                    ),
                    shiny::sliderInput("wrap", "Term Name Width", 15, 60, 25, 1),
                    shiny::numericInput("fontsize", "Font Size (pt)", 12, 4, 48, 0.5),
                    shiny::numericInput("lineheight", "Line Spacing", 0.9, 0.6, 1.5, 0.05),
                    shiny::checkboxGroupInput("show", "Label Lines",
                        c("Term Name" = "name", "Statistic" = "p"),
                        selected = c("name", "p"), inline = TRUE),
                    shiny::checkboxInput("show_counts",
                        "Show Enrichment Counts (GeneRatio | BgRatio)", TRUE),
                    shiny::checkboxInput("physics", "Elastic Layout (Physics)", FALSE),
                    shiny::checkboxInput("rearrange", "Re-arrange After Delete", TRUE),
                    shiny::checkboxInput("collapse", "Collapse Single-child Chains", FALSE),
                    shiny::fluidRow(
                        shiny::column(6, shiny::sliderInput("spring", "Spring Length", 50, 400, 150, 10)),
                        shiny::column(6, shiny::sliderInput("gravity", "Gravity", -250, -10, -60, 5))
                    ),
                    shiny::actionButton("relayout", "Auto-Arrange", width = "100%"),
                    shiny::div(style = "height:4px"),
                    shiny::fluidRow(
                        shiny::column(6, shiny::actionButton("undo", "Undo Delete", width = "100%")),
                        shiny::column(6, shiny::actionButton("reset", "Reset All", width = "100%"))
                    )
                ),
                shiny::div(class = "well",
                    shiny::h5("Colours"),
                    shiny::selectInput("preset", "Preset",
                        c("White-gold-rose (default)" = "warm", "Blue-white-red" = "bwr",
                          "AgriGO Style" = "agrigo", "Blue-Teal" = "teal",
                          "Custom" = "custom"), selected = "warm"),
                    shiny::fluidRow(
                        shiny::column(6, colourpicker::colourInput("col_low", "Low", "#FFFFFF")),
                        shiny::column(6, colourpicker::colourInput("col_mid", "Mid", "#E9C46A"))
                    ),
                    colourpicker::colourInput("col_high", "High", "#EA8379"),
                    shiny::fluidRow(
                        shiny::column(6, colourpicker::colourInput("col_ancestor", "Ancestor Fill", "#FFFFFF")),
                        shiny::column(6, colourpicker::colourInput("col_border", "Border", "#000000"))
                    ),
                    shiny::fluidRow(
                        shiny::column(6, colourpicker::colourInput("col_edge", "Is_a Edge", "#000000")),
                        shiny::column(6, colourpicker::colourInput("col_partof", "Part_of Edge", "#2878B5"))
                    ),
                    colourpicker::colourInput("col_bg", "Canvas", "#FFFFFF")
                )
            ),
            shiny::column(6,
                shiny::div(class = "graph-wrap",
                    visNetwork::visNetworkOutput("graph", height = "84vh"),
                    shiny::conditionalPanel(
                        condition = "output.has_data == false",
                        shiny::div(class = "welcome",
                            "Load the example data, or upload an annotation table and a gene list, then run enrichment.",
                            shiny::br(),
                            "Click to select, double-click to delete a branch, drag to arrange."
                        )
                    )
                ),
                shiny::div(class = "legend-note",
                    "Click to select | double-click to delete a branch | drag to arrange | export on the right",
                    shiny::textOutput("status", inline = TRUE),
                    shiny::textOutput("task_info", inline = TRUE)
                ),
                shiny::div(class = "well",
                    shiny::h5("Download Built-in GO Info"),
                    if (length(avail)) {
                        idx <- seq_along(avail)
                        chunks <- split(idx, ceiling(idx / 2))
                        shiny::tagList(lapply(chunks, function(ii) {
                            shiny::fluidRow(lapply(ii, function(i) {
                                shiny::column(6, shiny::downloadButton(
                                    .builtin_download_ids[[avail[[i]]]],
                                    avail[[i]], style = "width:100%"))
                            }))
                        }))
                    } else {
                        shiny::div(class = "legend-note",
                            "No built-in GO info files available.")
                    }
                )
            ),
            shiny::column(3,
                shiny::div(class = "right-panel",
                    shiny::div(class = "well",
                        shiny::h5("Selected Term"),
                        shiny::div(class = "sel-box", shiny::textOutput("sel_info")),
                        shiny::actionButton("delete_sel", "Delete Term and Orphaned Descendants",
                                            class = "btn-danger", width = "100%"),
                        shiny::div(style = "height:4px"),
                        shiny::actionButton("undo_last", "Undo Last Delete",
                                            class = "btn-warning", width = "100%")
                    ),
                    shiny::div(class = "well",
                        shiny::h5("FDR Colour Scale"),
                        shiny::fluidRow(
                            shiny::column(4, shiny::numericInput(
                                "fdr_blue", "Blue", 0.05, 1e-6, 1, 0.001)),
                            shiny::column(4, shiny::numericInput(
                                "fdr_white", "White", 0.01, 1e-6, 1, 0.001)),
                            shiny::column(4, shiny::numericInput(
                                "fdr_red", "Red", 0.0001, 1e-8, 1, 0.0001))
                        ),
                        shiny::uiOutput("legend")
                    ),
                    shiny::div(class = "well",
                        shiny::h5("Export Figure"),
                        shiny::selectInput("format", "Format",
                            c("PDF" = "pdf", "PNG" = "png", "SVG" = "svg",
                              "TIFF" = "tiff")),
                        shiny::checkboxInput("match_screen",
                            "Match on-screen layout (curved edges, screen spacing)",
                            TRUE),
                        shiny::checkboxInput("auto_size", "Auto Size", TRUE),
                        shiny::conditionalPanel(
                            condition = "input.auto_size == false",
                            shiny::fluidRow(
                                shiny::column(6, shiny::numericInput("width", "Width (in)", 12, 2, 100, 0.5)),
                                shiny::column(6, shiny::numericInput("height", "Height (in)", 8, 2, 100, 0.5))
                            )
                        ),
                        shiny::checkboxInput("export_legend",
                            "Include colour legend in the figure", TRUE),
                        shiny::numericInput("dpi", "DPI (PNG/TIFF)", 800, 72, 1200, 50),
                        shiny::downloadButton("download", "Download Current Layout",
                                              class = "btn-primary", style = "width:100%"),
                        shiny::div(class = "legend-note",
                            "The current arrangement is exported. A small collision-relaxation step prevents text overlap.")
                    ),
                    shiny::div(class = "well",
                        shiny::h5("Export Tables"),
                        shiny::downloadButton("download_csv",
                            "Enrichment Results (CSV)", style = "width:100%"),
                        shiny::div(style = "height:4px"),
                        shiny::downloadButton("download_prune",
                            "Pruning Decisions (CSV)", style = "width:100%")
                    ),
                    shiny::div(class = "well",
                        shiny::h5("Deleted Branch Roots"),
                        shiny::div(class = "status", shiny::textOutput("dropped_info"))
                    )
                )
            )
        )
    )
}

app_server <- function(input, output, session) {
    session$allowReconnect("force")
    .set_upload_limit()
    rv <- shiny::reactiveValues(
        anno = NULL, ego = NULL, dag = NULL,
        drop_roots = character(), positions = NULL, selected = NULL,
        relayout = 0L
    )

    output$has_data <- shiny::reactive(!is.null(rv$ego))
    shiny::outputOptions(output, "has_data", suspendWhenHidden = FALSE)

    show_error <- function(e) {
        shiny::showNotification(conditionMessage(e), type = "error", duration = 12)
    }

    show_anno_example <- function() {
        shiny::showModal(shiny::modalDialog(
            title = "GO Annotation Table Example",
            shiny::tags$img(src = "appwww/GO_info_example.png",
                            style = "max-width:100%;height:auto;"),
            easyClose = TRUE,
            footer = shiny::tagList(
                shiny::actionButton("cancel_browse", "Cancel"),
                shiny::tags$button(
                    "Choose File", type = "button",
                    class = "btn btn-primary",
                    onclick = "if(window.__annoPickerTimer){clearTimeout(window.__annoPickerTimer);window.__annoPickerTimer=null;} document.getElementById('anno_file').click();")
            )
        ))
    }
    shiny::observeEvent(input$browse_anno, {
        show_anno_example()
        session$sendCustomMessage("anno_picker_delay", list(delay = 6000))
    }, ignoreInit = TRUE)
    shiny::observeEvent(input$cancel_browse, shiny::removeModal(),
                        ignoreInit = TRUE)
    shiny::observeEvent(input$anno_file, shiny::removeModal(),
                        ignoreInit = TRUE)
    output$anno_file_name <- shiny::renderText({
        if (is.null(input$anno_file)) return("No File Selected")
        sprintf("%s (%.1f MB)", input$anno_file$name,
                input$anno_file$size / 1024^2)
    })

    load_example_data <- function() {
        tryCatch({
            res <- shiny::withProgress(message = "Loading Example Data",
                                       value = 0.2, {
                mini <- function() {
                    a <- read_go_annotation(
                        system.file("extdata", "example_mini_annotation.csv",
                                    package = "PruneTopGODAGviz"),
                        ont = "BP", verbose = FALSE)
                    g <- readLines(system.file("extdata", "example_mini_genes.txt",
                                               package = "PruneTopGODAGviz"))
                    d <- read_obo(system.file("extdata", "mini.obo",
                                              package = "PruneTopGODAGviz"),
                                  ont = "BP", verbose = FALSE)
                    if (!identical(input$propagate, FALSE)) {
                        a <- propagate_go_annotation(a, d, ont = "BP",
                                                     verbose = FALSE)
                    }
                    e <- enrich_go(g, a, ont = "BP", minGSSize = 10,
                                   maxGSSize = 50000, verbose = FALSE)
                    list(anno = a, genes = g, dag = d, ego = e, top = 5)
                }
                find_data <- function(name) {
                    p <- .builtin_file(name)
                    if (file.exists(p)) return(p)
                    for (d in c(file.path(getwd(), "data"),
                                file.path(getwd(), "..", "data"))) {
                        q <- file.path(d, name)
                        if (file.exists(q)) return(q)
                    }
                    ""
                }
                res <- tryCatch({
                    at_anno <- find_data("At_GO_info_panther.csv")
                    if (!nzchar(at_anno)) {
                        stop("Built-in At annotation not found; install the package data or set PruneTopGODAGviz.data.")
                    }
                    at_genes <- find_data("At_diffpeak_genes.GID")
                    if (!nzchar(at_genes)) {
                        at_genes <- system.file("extdata",
                            "example_at_genes.txt",
                            package = "PruneTopGODAGviz")
                    }
                    a <- read_go_annotation(at_anno, keytype = "gene",
                                            ont = "BP", verbose = FALSE)
                    g <- readLines(at_genes)
                    d <- go_dag("BP", verbose = FALSE)
                    if (!identical(input$propagate, FALSE)) {
                        a <- propagate_go_annotation(a, d, ont = "BP",
                                                     verbose = FALSE)
                    }
                    e <- enrich_go(g, a, ont = "BP",
                                   maxGSSize = 50000, verbose = FALSE)
                    list(anno = a, genes = g, dag = d, ego = e, top = 15)
                }, error = function(e) mini())
                shiny::incProgress(0.8)
                res
            })
            rv$anno <- res$anno
            rv$ego <- res$ego
            rv$dag <- res$dag
            rv$drop_roots <- character()
            rv$positions <- NULL
            rv$selected <- NULL
            rv$relayout <- rv$relayout + 1L
            shiny::updateNumericInput(session, "top_n", value = res$top)
            shiny::updateCheckboxInput(session, "physics", value = FALSE)
            shiny::updateSelectInput(session, "keytype", selected = "gene")
            shiny::updateSelectInput(session, "ont", selected = "BP")
            shiny::showNotification(
                sprintf("Example Data Loaded: %d genes", length(res$genes)),
                type = "message")
        }, error = show_error)
    }
    shiny::observeEvent(input$load_example, load_example_data())
    shiny::observe({
        query <- shiny::parseQueryString(session$clientData$url_search)
        if (identical(query$demo, "1") && is.null(rv$ego)) {
            load_example_data()
        }
    })

    enrich_job <- function(args) {
        tryCatch({
            dag <- PruneTopGODAGviz::go_dag(args$ont, obo = args$obo,
                                            verbose = FALSE)
            anno <- PruneTopGODAGviz::read_go_annotation(
                args$path, keytype = args$keytype, ont = args$ont,
                verbose = FALSE)
            check_ids <- utils::getFromNamespace(
                ".check_id_match", "PruneTopGODAGviz")
            chk <- check_ids(args$genes, anno$universe[[args$ont]])
            genes <- chk$ids
            if (isTRUE(args$propagate)) {
                anno <- PruneTopGODAGviz::propagate_go_annotation(
                    anno, dag, ont = args$ont, verbose = FALSE)
            }
            ego <- PruneTopGODAGviz::enrich_go(
                genes, anno, ont = args$ont, pvalueCutoff = 1,
                qvalueCutoff = 1, pAdjustMethod = "BH",
                minGSSize = args$min_gs, maxGSSize = args$max_gs,
                verbose = FALSE)
            list(ok = TRUE, anno = anno, ego = ego, dag = dag, chk = chk)
        }, error = function(e) {
            list(ok = FALSE, error = conditionMessage(e),
                 id_error = inherits(e, "id_check_error"))
        })
    }
    enrich_task <- shiny::ExtendedTask$new(function(args) {
        future::future(enrich_job(args), seed = TRUE)
    })
    task_done <- shiny::reactiveVal(FALSE)

    shiny::observeEvent(input$run, {
        tryCatch({
            .check_fdr_breaks(c(input$fdr_blue, input$fdr_white,
                                 input$fdr_red))
            shiny::req(input$gene_file)
            genes <- readLines(input$gene_file$datapath)
            keytype <- input$keytype
            if (nzchar(input$builtin_anno)) {
                anno_name <- switch(input$builtin_anno,
                    arabidopsis_thaliana = "At_GO_info_panther.csv",
                    zea_mays = "Maize_GO_info_panther.csv",
                    oryza_sativa = "Rice_GO_info_panther.csv",
                    glycine_max = "Soybean_GO_info_panther.csv",
                    homo_sapiens = "Homo_sapiens_GO_info_panther.csv",
                    mus_musculus = "Mus_musculus_GO_info_panther.csv",
                    stop("Unknown built-in annotation: ", input$builtin_anno)
                )
                path <- .builtin_file(anno_name)
                if (!file.exists(path)) {
                    stop("Built-in annotation not found at ", path)
                }
            } else {
                shiny::req(input$anno_file)
                path <- input$anno_file$datapath
            }
            task_done(FALSE)
            enrich_task$invoke(list(
                path = path, genes = genes, keytype = keytype,
                ont = input$ont, obo = getOption("PruneTopGODAGviz.obo", NULL),
                propagate = input$propagate, min_gs = input$min_gs,
                max_gs = input$max_gs))
            session$sendCustomMessage("portal_enrichment_state",
                                      list(state = "running"))
        }, error = show_error)
    })

    shiny::observe({
        if (!identical(enrich_task$status(), "success") || isTRUE(task_done())) {
            return()
        }
        task_done(TRUE)
        res <- enrich_task$result()
        if (!isTRUE(res$ok)) {
            if (isTRUE(res$id_error)) {
                shiny::showModal(shiny::modalDialog(
                    title = "ID Check Failed",
                    shiny::HTML(paste0(
                        res$error, "<br><br>Selected annotation: <b>",
                        if (nzchar(input$builtin_anno)) input$builtin_anno else "uploaded",
                        "</b><br>ID Type: <b>", input$keytype,
                        "</b><br>Required: at least 10 unique IDs, 10 matches and 20% match rate.")),
                    easyClose = TRUE,
                    footer = shiny::modalButton("Close")))
            } else {
                show_error(res$error)
            }
            session$sendCustomMessage("portal_enrichment_state",
                                      list(state = "failed"))
            return()
        }
        rv$dag <- res$dag
        rv$anno <- res$anno
        rv$ego <- res$ego
        rv$drop_roots <- character()
        rv$positions <- NULL
        rv$selected <- NULL
        rv$relayout <- rv$relayout + 1L
        shiny::showNotification(
            sprintf("Enrichment Finished: %d terms tested",
                    nrow(rv$ego$result)), type = "message")
        session$sendCustomMessage("portal_enrichment_state",
                                  list(state = "finished",
                                       terms = nrow(rv$ego$result)))
    })
    pruned <- shiny::reactive({
        shiny::req(rv$ego, rv$dag)
        run <- function(cut) suppressWarnings(prune_go_dag(
            rv$ego, dag = rv$dag,
            top_n = input$top_n, p_cutoff = cut,
            drop_subtree = if (length(rv$drop_roots)) rv$drop_roots else NULL,
            collapse_chains = isTRUE(input$collapse), cluster_jaccard = NULL,
            max_nodes = 300, wrap = input$wrap,
            counts = "propagated", verbose = FALSE))
        tryCatch(run(input$p_cutoff), error = function(e) {
            if (grepl("No term passes", conditionMessage(e), fixed = TRUE)) {
                run(1)
            } else {
                stop(e)
            }
        })
    })

    shiny::observeEvent(input$node_click, {
        pos <- .positions_from_js(input$node_click$positions)
        if (!is.null(pos)) rv$positions <- pos
        rv$selected <- input$node_click$id
    })

    delete_node <- function(id) {
        if (is.null(id) || !length(id) || is.na(id) || !nzchar(id)) {
            shiny::showNotification("Select a term first.", type = "message",
                                    duration = 3)
            return()
        }
        d <- pruned()
        if (!id %in% d$nodes$id) return()
        removed <- .cascade_orphans(rv$dag, id)
        selected_displayed <- d$nodes$id[d$nodes$selected]
        if (!length(selected_displayed) ||
            all(selected_displayed %in% removed)) {
            shiny::showNotification(
                "Cannot delete: this would remove every enriched term.",
                type = "error", duration = 6)
            return()
        }
        rv$drop_roots <- unique(c(rv$drop_roots, id))
        if (isTRUE(input$rearrange)) rv$positions <- NULL
        rv$selected <- NULL
        n_desc <- length(setdiff(intersect(removed, d$nodes$id), id))
        shiny::showNotification(
            if (n_desc > 0L) {
                sprintf("Deleted %s and %d orphaned descendant(s) (undo available)",
                        id, n_desc)
            } else {
                sprintf("Deleted %s (undo available)", id)
            },
            type = "warning", duration = 5)
    }
    shiny::observeEvent(input$node_dblclick, {
        pos <- .positions_from_js(input$node_dblclick$positions)
        if (!is.null(pos)) rv$positions <- pos
        delete_node(input$node_dblclick$id)
    })
    shiny::observeEvent(input$delete_sel, delete_node(rv$selected))

    undo_delete <- function() {
        if (length(rv$drop_roots)) {
            rv$drop_roots <- utils::head(rv$drop_roots, -1)
            if (isTRUE(input$rearrange)) rv$positions <- NULL
            shiny::showNotification("Last Deletion Undone", type = "message")
        }
    }
    shiny::observeEvent(input$undo, undo_delete())
    shiny::observeEvent(input$undo_last, undo_delete())

    shiny::observeEvent(input$reset, {
        rv$drop_roots <- character()
        rv$positions <- NULL
        rv$selected <- NULL
        shiny::showNotification("Reset", type = "message")
    })

    shiny::observeEvent(input$relayout, {
        rv$positions <- NULL
        rv$relayout <- rv$relayout + 1L
        shiny::showNotification("Re-arranging Layout", type = "message")
    })

    shiny::observeEvent(input$preset, {
        p <- switch(input$preset,
            bwr = list(fill = c("#2878B5", "#FFFFFF", "#c82423"),
                       ancestor = "#FFFFFF", part_of = "#2878B5"),
            agrigo = list(fill = c("#BB9727", "#F27970", "#c82423"),
                          ancestor = "#FFFFFF", part_of = "#54B345"),
            teal = list(fill = c("#CFEAF1", "#96CCCB", "#2878B5"),
                        ancestor = "#FFFFFF", part_of = "#32B897"),
            warm = list(fill = c("#FFFFFF", "#E9C46A", "#EA8379"),
                        ancestor = "#FFFFFF", part_of = "#E53528"),
            custom = NULL)
        if (!is.null(p)) {
            colourpicker::updateColourInput(session, "col_low", value = p$fill[1])
            colourpicker::updateColourInput(session, "col_mid", value = p$fill[2])
            colourpicker::updateColourInput(session, "col_high", value = p$fill[3])
            colourpicker::updateColourInput(session, "col_ancestor", value = p$ancestor)
            colourpicker::updateColourInput(session, "col_partof", value = p$part_of)
        }
    }, ignoreInit = TRUE)

    output$graph <- visNetwork::renderVisNetwork({
        if (is.null(rv$ego)) {
            return(visNetwork::visNetwork(
                data.frame(id = character()),
                data.frame(from = character(), to = character())))
        }
        shiny::req(rv$relayout >= 0L)
        d <- pruned()
        .build_vis(
            d,
            positions = shiny::isolate(rv$positions),
            physics = isTRUE(input$physics),
            spring = input$spring,
            gravity = input$gravity,
            fontsize = input$fontsize,
            show_name = "name" %in% input$show,
            show_counts = isTRUE(input$show_counts),
            show_p = "p" %in% input$show,
            colors = list(low = input$col_low, mid = input$col_mid,
                          high = input$col_high, ancestor = input$col_ancestor,
                          border = input$col_border, edge = input$col_edge,
                          part_of = input$col_partof, bg = input$col_bg,
                          fdr_breaks = c(input$fdr_blue, input$fdr_white,
                                         input$fdr_red))
        )
    })

    shiny::observeEvent(input$graph_drag, {
        visNetwork::visNetworkProxy("graph") |> visNetwork::visGetPositions()
    })
    shiny::observeEvent(input$graph_stab, {
        visNetwork::visNetworkProxy("graph") |> visNetwork::visGetPositions()
    })

    shiny::observeEvent(input$graph_positions, {
        pos <- .positions_from_js(input$graph_positions)
        if (!is.null(pos)) rv$positions <- pos
    })

    output$sel_info <- shiny::renderText({
        if (is.null(rv$selected)) return("Click a term in the graph to select it")
        d <- pruned()
        row <- d$nodes[d$nodes$id == rv$selected, , drop = FALSE]
        if (!nrow(row)) return(paste0(rv$selected, "\n(not in the current graph)"))
        txt <- c(row$id, row$name,
                 if (!is.na(row$p.adjust)) paste0("FDR: ", .fmt_p(row$p.adjust)),
                 row$label_l3,
                 if (isTRUE(row$selected)) "Significant Term" else "Ancestor Term")
        paste(txt, collapse = "\n")
    })

    output$task_info <- shiny::renderText({
        if (identical(enrich_task$status(), "running")) {
            " | enrichment running in background"
        } else {
            ""
        }
    })

    output$status <- shiny::renderText({
        if (is.null(rv$ego)) return("")
        d <- pruned()
        passes <- any(rv$ego$result$p.adjust < input$p_cutoff, na.rm = TRUE)
        note <- if (passes) "" else " | no term passes cutoff; showing top nominal terms"
        sprintf(" | %d terms shown (%d selected)%s", nrow(d$nodes),
                sum(d$nodes$selected), note)
    })

    output$legend <- shiny::renderUI({
        fdr_labels <- tryCatch(
            .fmt_break(.check_fdr_breaks(c(input$fdr_blue, input$fdr_white,
                                           input$fdr_red))),
            error = function(e) c("0.05", "0.01", "0.0001"))
        shiny::tagList(
            shiny::div(style = sprintf(
                "height:12px;border-radius:3px;border:1px solid #d9dde3;background:linear-gradient(to right,%s,%s,%s);",
                input$col_low, input$col_mid, input$col_high)),
            shiny::div(style = "display:flex;justify-content:space-between;font-size:11px;color:#56606b;margin-top:3px;",
                       shiny::span(fdr_labels[1]),
                       shiny::span(fdr_labels[2]),
                       shiny::span(fdr_labels[3])),
            shiny::HTML(sprintf(
                '<svg width="100%%" height="46" viewBox="0 0 170 46" style="margin-top:4px">
                   <defs>
                     <marker id="ar1" markerWidth="7" markerHeight="7" refX="5.5" refY="3.5" orient="auto"><path d="M0,0 L7,3.5 L0,7 z" fill="%s"/></marker>
                     <marker id="ar2" markerWidth="7" markerHeight="7" refX="5.5" refY="3.5" orient="auto"><path d="M0,0 L7,3.5 L0,7 z" fill="%s"/></marker>
                   </defs>
                   <line x1="2" y1="12" x2="42" y2="12" stroke="%s" stroke-width="1.4" marker-end="url(#ar1)"/>
                   <text x="52" y="16" font-family="Arial" font-size="11" fill="#333">is_a</text>
                   <line x1="2" y1="32" x2="42" y2="32" stroke="%s" stroke-width="1.4" stroke-dasharray="5,3" marker-end="url(#ar2)"/>
                   <text x="52" y="36" font-family="Arial" font-size="11" fill="#333">part_of</text>
                 </svg>', input$col_edge, input$col_partof, input$col_edge,
                input$col_partof)))
    })

    output$dropped_info <- shiny::renderText({
        if (!length(rv$drop_roots)) return("None")
        paste(rv$drop_roots, collapse = "\n")
    })

    output$download <- shiny::downloadHandler(
        filename = function() {
            # Include seconds so repeated exports in the same minute do not collide.
            paste0("GO_DAG_", input$ont, "_",
                   format(Sys.time(), "%Y_%m%d_%H%M%S"), ".", input$format)
        },
        content = function(file) {
            d <- pruned()
            pos <- rv$positions
            # Default export reproduces the on-screen visNetwork figure:
            # same spacing, line height and cubic-bezier edges.
            match_screen <- isTRUE(input$match_screen)
            gap_level <- if (match_screen) 0.10 * input$fontsize else 0.28
            gap_sibling <- if (match_screen) 0.028 * input$fontsize else 0.12
            lh_used <- if (match_screen) 1.05 else input$lineheight
            edge_style <- if (match_screen) "bezier" else "hybrid"
            lay <- if (!is.null(pos) && nrow(pos)) {
                layout_from_positions(d, pos, fontsize = input$fontsize,
                                      lineheight = lh_used,
                                      family = "Arial", rankdir = "TB",
                                      verbose = FALSE)
            } else {
                layout_go_dag(d, fontsize = input$fontsize,
                              lineheight = lh_used, family = "Arial",
                              level_gap = gap_level, sibling_gap = gap_sibling,
                              edge_style = edge_style, verbose = FALSE)
            }
            p <- plot_go_dag(
                d, layout = lay,
                fill = c(input$col_low, input$col_mid, input$col_high),
                ancestor_fill = input$col_ancestor,
                bg = input$col_bg, font_family = "Arial",
                edge_col = c(is_a = input$col_edge,
                             part_of = input$col_partof,
                             unknown = input$col_edge),
                fdr_breaks = c(input$fdr_blue, input$fdr_white,
                               input$fdr_red),
                fontsize = input$fontsize, lineheight = lh_used,
                wrap = input$wrap,
                show_name = "name" %in% input$show,
                show_counts = isTRUE(input$show_counts),
                show_p = "p" %in% input$show,
                legend = isTRUE(input$export_legend))
            sz <- plot_go_dag_size(d, layout = lay)
            if (isTRUE(input$auto_size)) {
                w <- sz$width
                h <- sz$height
            } else {
                # Honour the requested page size. When Match is off, adapt the
                # spacing and only reduce type when the page is genuinely tight;
                # when Match is on, preserve the panel geometry and scale uniformly.
                w <- as.numeric(input$width)
                h <- as.numeric(input$height)
                rd_use <- "TB"
                if (is.finite(w) && is.finite(h) && w > 0 && h > 0) {
                    if (match_screen) {
                        sz0 <- plot_go_dag_size(d, layout = lay)
                        fit0 <- min(w / sz0$width, h / sz0$height)
                        if (!is.finite(fit0) || fit0 <= 0) fit0 <- 1
                        lay <- scale_go_layout(lay, fit0)
                        fs_use <- input$fontsize * fit0
                        margin_use <- 0.08 * fit0
                        export_note <- sprintf(
                            paste0("Exported %.1f x %.1f in at %.0f%% of the on-screen layout",
                                   " (text %.1f pt). The figure matches the panel."),
                            w, h, 100 * fit0, fs_use)
                    } else if (is.null(pos) || !nrow(pos)) {
                        fitted <- fit_go_layout_to_page(
                            d, width_in = w, height_in = h,
                            fontsize = input$fontsize,
                            lineheight = lh_used, family = "Arial",
                            level_gap = gap_level, sibling_gap = gap_sibling,
                            edge_style = edge_style)
                        lay <- fitted$layout
                        fs_use <- fitted$fontsize
                        margin_use <- fitted$margin
                        rd_use <- fitted$rankdir
                        if (isTRUE(fitted$shrunk)) {
                            export_note <- sprintf(
                                paste0("Exported %.1f x %.1f in. The page is tight for %d terms:",
                                       " type reduced from %.1f to %.1f pt (direction %s). Enlarge Width/Height",
                                       " or reduce Top Significant Terms to keep the type size."),
                                w, h, nrow(d$nodes), input$fontsize, fs_use, rd_use)
                        } else {
                            export_note <- sprintf(
                                "Exported %.1f x %.1f in, text %.1f pt, direction %s.",
                                w, h, fs_use, rd_use)
                        }
                    } else {
                        sz2 <- plot_go_dag_size(d, layout = lay)
                        fit <- min(w / sz2$width, h / sz2$height)
                        if (!is.finite(fit) || fit <= 0) fit <- 1
                        lay <- scale_go_layout(lay, fit)
                        fs_use <- input$fontsize * fit
                        margin_use <- 0.08 * fit
                    }
                    p <- plot_go_dag(
                        d, layout = lay,
                        fill = c(input$col_low, input$col_mid, input$col_high),
                        ancestor_fill = input$col_ancestor,
                        bg = input$col_bg, font_family = "Arial",
                        edge_col = c(is_a = input$col_edge,
                                     part_of = input$col_partof,
                                     unknown = input$col_edge),
                        fdr_breaks = c(input$fdr_blue, input$fdr_white,
                                       input$fdr_red),
                        fontsize = fs_use,
                        margin = margin_use,
                        rankdir = rd_use,
                        lineheight = lh_used,
                        wrap = input$wrap,
                        show_name = "name" %in% input$show,
                        show_counts = isTRUE(input$show_counts),
                        show_p = "p" %in% input$show,
                        legend = isTRUE(input$export_legend))
                }
            }
            prefix <- sub("\\.(pdf|png|svg|tiff)$", "", file)
            save_go_dag(d, prefix = prefix, plot = p, width = w, height = h,
                        dpi = input$dpi, formats = input$format,
                        overwrite = TRUE, bg = input$col_bg)
            if (!isTRUE(input$auto_size) && exists("fs_use")) {
                msg <- if (exists("export_note")) {
                    export_note
                } else {
                    sprintf("Exported %.1f x %.1f in, text %.1f pt.",
                            as.numeric(w), as.numeric(h), fs_use)
                }
                shiny::showNotification(msg, type = "message", duration = 12)
            }
        }
    )
    output$download_csv <- shiny::downloadHandler(
        filename = function() {
            paste0("GO_enrichment_", input$ont, "_",
                   format(Sys.time(), "%Y_%m%d_%H%M"), ".csv")
        },
        content = function(file) {
            shiny::req(rv$ego)
            utils::write.csv(rv$ego$result, file, row.names = FALSE)
        }
    )

    output$download_prune <- shiny::downloadHandler(
        filename = function() {
            paste0("GO_pruning_", input$ont, "_",
                   format(Sys.time(), "%Y_%m%d_%H%M"), ".csv")
        },
        content = function(file) {
            write_prune_table(pruned(), file)
        }
    )

    output$download_arabidopsis <- shiny::downloadHandler(
        filename = function() basename(.builtin_file("At_GO_info_panther.csv")),
        content = function(file) {
            src <- .builtin_file("At_GO_info_panther.csv")
            if (!file.exists(src)) {
                stop("At_GO_info_panther.csv not found at ", src)
            }
            file.copy(src, file, overwrite = TRUE)
        }
    )

    output$download_zea_mays <- shiny::downloadHandler(
        filename = function() basename(.builtin_file("Maize_GO_info_panther.csv")),
        content = function(file) {
            src <- .builtin_file("Maize_GO_info_panther.csv")
            if (!file.exists(src)) {
                stop("Maize_GO_info_panther.csv not found at ", src)
            }
            file.copy(src, file, overwrite = TRUE)
        }
    )

    output$download_rice <- shiny::downloadHandler(
        filename = function() basename(.builtin_file("Rice_GO_info_panther.csv")),
        content = function(file) {
            src <- .builtin_file("Rice_GO_info_panther.csv")
            if (!file.exists(src)) {
                stop("Rice_GO_info_panther.csv not found at ", src)
            }
            file.copy(src, file, overwrite = TRUE)
        }
    )

    output$download_soybean <- shiny::downloadHandler(
        filename = function() basename(.builtin_file("Soybean_GO_info_panther.csv")),
        content = function(file) {
            src <- .builtin_file("Soybean_GO_info_panther.csv")
            if (!file.exists(src)) {
                stop("Soybean_GO_info_panther.csv not found at ", src)
            }
            file.copy(src, file, overwrite = TRUE)
        }
    )

    output$download_human <- shiny::downloadHandler(
        filename = function() basename(.builtin_file("Homo_sapiens_GO_info_panther.csv")),
        content = function(file) {
            src <- .builtin_file("Homo_sapiens_GO_info_panther.csv")
            if (!file.exists(src)) {
                stop("Homo_sapiens_GO_info_panther.csv not found at ", src)
            }
            file.copy(src, file, overwrite = TRUE)
        }
    )

    output$download_mouse <- shiny::downloadHandler(
        filename = function() basename(.builtin_file("Mus_musculus_GO_info_panther.csv")),
        content = function(file) {
            src <- .builtin_file("Mus_musculus_GO_info_panther.csv")
            if (!file.exists(src)) {
                stop("Mus_musculus_GO_info_panther.csv not found at ", src)
            }
            file.copy(src, file, overwrite = TRUE)
        }
    )
}

.positions_from_js <- function(pos) {
    if (is.null(pos) || !is.list(pos) || !length(pos)) return(NULL)
    ids <- names(pos)
    if (is.null(ids)) ids <- as.character(seq_along(pos))
    xy <- t(vapply(pos, function(p) {
        px <- if (!is.null(p$x)) p$x else p[[1]]
        py <- if (!is.null(p$y)) p$y else p[[2]]
        c(as.numeric(px), as.numeric(py))
    }, numeric(2)))
    data.frame(id = ids, x = xy[, 1], y = xy[, 2], stringsAsFactors = FALSE)
}

.vis_colours <- function(nodes, colors,
                            fdr_breaks = c(0.05, 0.01, 0.0001)) {
    sel <- nodes$selected & !is.na(nodes$p.adjust)
    fdr_breaks <- .check_fdr_breaks(fdr_breaks)
    pos <- .fdr_position(nodes$p.adjust, fdr_breaks[1], fdr_breaks[2],
                         fdr_breaks[3])
    pal <- grDevices::colorRampPalette(c(colors$low, colors$mid,
                                         colors$high))(101)
    idx <- ifelse(sel, 1 + round(100 * pos), NA_real_)
    fill <- ifelse(sel, pal[pmin(100, pmax(1, idx))], colors$ancestor)
    fill[is.na(fill)] <- colors$ancestor
    fill
}

.vis_labels <- function(nodes, fontsize, show_name, show_counts, show_p) {
    l1 <- nodes$id
    if (isTRUE(show_p)) {
        has <- nodes$selected & !is.na(nodes$p.adjust)
        l1[has] <- paste0(nodes$id[has], " (", .fmt_p(nodes$p.adjust[has]), ")")
    }
    l2 <- if (isTRUE(show_name)) nodes$label_l2 else rep("", nrow(nodes))
    l3 <- if (isTRUE(show_counts)) nodes$label_l3 else rep("", nrow(nodes))
    out <- paste(l1, l2, l3, sep = "\n")
    gsub("\n+", "\n", out)
}

.vis_nodes <- function(d, col_low, col_mid, col_high, col_ancestor,
                       col_border, fontsize, show_name, show_counts, show_p,
                       fdr_breaks = c(0.05, 0.01, 0.0001)) {
    nodes <- d$nodes
    colors <- list(low = col_low, mid = col_mid, high = col_high,
                   ancestor = col_ancestor)
    data.frame(
        id = nodes$id,
        label = .vis_labels(nodes, fontsize, show_name, show_counts, show_p),
        title = .vis_tooltips(d),
        color.background = .vis_colours(nodes, colors, fdr_breaks),
        color.border = col_border,
        color.highlight.background = "#F6CAE5",
        color.highlight.border = col_border,
        stringsAsFactors = FALSE)
}

.vis_edges <- function(d, col_edge, col_partof) {
    e <- d$edges
    data.frame(
        from = e$parent, to = e$child,
        dashes = e$type == "part_of",
        color = ifelse(e$type == "part_of", col_partof, col_edge),
        stringsAsFactors = FALSE)
}

.vis_tooltips <- function(d) {
    sets <- split(d$annotation$pairs$transcript, d$annotation$pairs$go_id)
    genes <- lapply(d$nodes$id, function(id) {
        intersect(sets[[id]] %||% character(), d$gene)
    })
    vapply(seq_len(nrow(d$nodes)), function(i) {
        g <- genes[[i]]
        shown <- utils::head(g, 12)
        more <- if (length(g) > length(shown)) {
            paste0(" +", length(g) - length(shown), " more")
        } else ""
        paste0("<div style='max-width:340px;white-space:normal'>",
               .html_escape(d$nodes$id[i]), "<br>",
               .html_escape(d$nodes$name[i]), "<br>",
               if (!is.na(d$nodes$p.adjust[i])) {
                   paste0("FDR: ", .html_escape(.fmt_p(d$nodes$p.adjust[i])),
                          "<br>")
               } else "",
               .html_escape(d$nodes$label_l3[i]),
               if (nzchar(d$nodes$collapsed_terms[i] %||% "")) {
                   paste0("<br>Hidden: ",
                          .html_escape(d$nodes$collapsed_terms[i]))
               } else "",
               if (length(g)) {
                   paste0("<br>Genes: ",
                          .html_escape(paste(shown, collapse = ", ")), more)
               } else "",
               "</div>")
    }, character(1))
}

.build_vis <- function(d, positions = NULL, physics = TRUE, spring = 150,
                       gravity = -60, fontsize = 7, show_name = TRUE,
                       show_counts = TRUE, show_p = TRUE, colors) {
    nodes <- .vis_nodes(d, colors$low, colors$mid, colors$high,
                        colors$ancestor, colors$border, fontsize,
                        show_name, show_counts, show_p,
                        fdr_breaks = colors$fdr_breaks %||%
                            c(0.05, 0.01, 0.0001))
    edges <- .vis_edges(d, colors$edge, colors$part_of)

    nodes$x <- NA_real_
    nodes$y <- NA_real_
    if (!is.null(positions) && nrow(positions)) {
        m <- match(nodes$id, positions$id)
        ok <- !is.na(m)
        nodes$x[ok] <- positions$x[m[ok]]
        nodes$y[ok] <- positions$y[m[ok]]
    }
    if (anyNA(nodes$x) || anyNA(nodes$y)) {
        base <- layout_go_dag(d, fontsize = fontsize, family = "Arial",
                              level_gap = 0.10 * fontsize,
                              sibling_gap = 0.028 * fontsize,
                              lineheight = 1.05, verbose = FALSE)
        bx <- stats::setNames(base$nodes$x * 96, base$nodes$id)
        by <- stats::setNames(base$nodes$y * 96, base$nodes$id)
        miss <- is.na(nodes$x) | is.na(nodes$y)
        nodes$x[miss] <- bx[nodes$id[miss]]
        nodes$y[miss] <- by[nodes$id[miss]]
    }

    visNetwork::visNetwork(nodes, edges, background = colors$bg,
                           width = "100%", height = "100%") |>
        visNetwork::visNodes(
            shape = "box",
            font = list(size = fontsize * 1.333, face = "Arial", vadjust = 0),
            margin = 6, borderWidth = 1, borderWidthSelected = 2,
            shapeProperties = list(borderRadius = 0),
            shadow = FALSE, chosen = TRUE) |>
        visNetwork::visEdges(
            arrows = list(to = list(enabled = TRUE, scaleFactor = 0.8)),
            smooth = list(enabled = TRUE, type = "cubicBezier",
                          forceDirection = "vertical", roundness = 0.5),
            arrowStrikethrough = FALSE, selectionWidth = 1.5) |>
        visNetwork::visPhysics(
            enabled = physics, solver = "forceAtlas2Based",
            forceAtlas2Based = list(
                gravitationalConstant = gravity, centralGravity = 0.005,
                springLength = spring, springConstant = 0.08,
                damping = 0.5, avoidOverlap = 1),
            stabilization = list(enabled = TRUE, iterations = 200,
                                 updateInterval = 25)) |>
        visNetwork::visInteraction(
            dragNodes = TRUE, dragView = TRUE, zoomView = TRUE, hover = TRUE,
            tooltipDelay = 120, navigationButtons = TRUE, keyboard = TRUE,
            multiselect = FALSE) |>
        visNetwork::visOptions(
            highlightNearest = list(enabled = TRUE, degree = 1, hover = FALSE),
            nodesIdSelection = FALSE) |>
        visNetwork::visEvents(
            click = paste0(
                "function(params){ if(params.nodes.length){",
                "var w = HTMLWidgets.find('#graph');",
                "var p = (w && w.network) ? w.network.getPositions() : null;",
                "Shiny.setInputValue('node_click',",
                "{id: params.nodes[0], positions: p, nonce: Math.random()},",
                "{priority: 'event'}); } }"),
            doubleClick = paste0(
                "function(params){ if(params.nodes.length){",
                "var w = HTMLWidgets.find('#graph');",
                "var p = (w && w.network) ? w.network.getPositions() : null;",
                "Shiny.setInputValue('node_dblclick',",
                "{id: params.nodes[0], positions: p, nonce: Math.random()},",
                "{priority: 'event'}); } }"),
            dragEnd = paste0(
                "function(){ Shiny.setInputValue('graph_drag',",
                "Math.random(), {priority: 'event'}); }"),
            stabilizationIterationsDone = paste0(
                "function(){ var w = HTMLWidgets.find('#graph');",
                "if (w && w.network) { w.network.setOptions({physics:{enabled:false}}); }",
                "Shiny.setInputValue('graph_stab', Math.random(),",
                "{priority: 'event'}); }"))
}
