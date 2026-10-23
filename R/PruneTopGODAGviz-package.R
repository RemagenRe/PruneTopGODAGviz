#' PruneTopGODAGviz: publication-quality pruned GO DAGs
#'
#' The package reads transcript-level GO annotation tables, runs
#' over-representation analysis with a self-contained hypergeometric test and draws an
#' agriGO-style directed acyclic graph that can be pruned, collapsed and
#' rendered without text overlap.
#'
#' @section Main functions:
#' * [read_go_annotation()] reads a transcript--GO table.
#' * [go_dag()] builds the GO hierarchy (OBO or GO.db).
#' * [enrich_go()] runs a self-contained hypergeometric over-representation test.
#' * [prune_go_dag()] selects terms and prunes branches for display.
#' * [plot_go_dag()] and [save_go_dag()] render PDF and PNG figures.
#' * [run_go_workflow()] does all of the above for BP/MF/CC.
#'
#' @importFrom rlang .data
#' @import data.table
#' @importFrom stats setNames
#' @importFrom utils head tail globalVariables
#' @keywords internal
"_PACKAGE"

utils::globalVariables(c("parent", "child"))
