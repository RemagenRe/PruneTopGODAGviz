# PruneTopGODAGviz 0.1.4

* Enable Shiny forced reconnection for local desktop use, so a backgrounded or
  temporarily disconnected browser tab recovers instead of remaining on the
  grey disconnected overlay.
# PruneTopGODAGviz 0.1.3

* Ensure background `future` workers inherit the clean project library path.
* Keep `getOption("PruneTopGODAGviz.obo")` and uploaded annotation paths
  visible to the worker.
# PruneTopGODAGviz 0.1.2

* Enrichment now runs in a background `ExtendedTask`/future worker, so the
  Shiny interface stays responsive during long calculations and after the
  browser tab is backgrounded.
* Added `future` and `promises` to the optional Shiny dependencies.
* Updated dependency instructions for GO.db/AnnotationDbi and background
  enrichment.
# PruneTopGODAGviz 0.1.1

* `GO.db` and `AnnotationDbi` are now mandatory dependencies, installed from
  Bioconductor.
* README now documents the exact dependency installation order and the
  additional Shiny/export packages required by `run_app()`.
* No functional or statistical changes.
# PruneTopGODAGviz 0.1.0

* First public release.
* Transcript- and gene-level GO over-representation analysis with a
  self-contained hypergeometric test.
* Pruned GO DAG layout, chain collapsing, Jaccard merging and
  publication-quality PDF/PNG/SVG/TIFF export.
* Shiny editor with built-in compressed annotation tables for Arabidopsis,
  maize, rice, soybean, human and mouse.
* Example fixtures, data provenance and citation metadata included.