# PruneTopGODAGviz

Publication-quality Gene Ontology DAGs from over-representation results, with
display pruning that leaves the statistics unchanged.

`PruneTopGODAGviz` reads transcript- or gene-level GO annotation tables, runs a
self-contained hypergeometric over-representation test, and draws an
agriGO-style DAG with user-controlled branch pruning, chain collapsing and
redundant-leaf merging. It also provides a Shiny editor and PDF/PNG/SVG/TIFF
export.

### Dependencies

`GO.db` and `AnnotationDbi` are required package dependencies. Install them
from Bioconductor before installing PruneTopGODAGviz:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("GO.db", "AnnotationDbi"))
```

Installing the tarball with `repos = NULL` still does not install all Suggests.
Before `run_app()`, install the Shiny/export packages explicitly:

```r
install.packages(c(
  "shiny", "visNetwork", "colourpicker",
  "ragg", "svglite", "ggiraph", "htmlwidgets", "htmltools"
))
```

Minimum for `run_app()`: `shiny`, `visNetwork`, `colourpicker`.
`ragg`/`svglite` are needed for PNG/SVG export; `ggiraph`/`htmlwidgets` for
interactive HTML. `GO.db`/`AnnotationDbi` are mandatory dependencies and are
installed via BiocManager above.

## Installation

Requires R >= 4.2.

```r
install.packages("pak")
pak::pkg_install("RemagenRe/PruneTopGODAGviz", dependencies = TRUE)

# reproducible release
pak::pkg_install("RemagenRe/PruneTopGODAGviz@v0.1.1")
```

Or install the release tarball:

```r
install.packages(
  "https://github.com/RemagenRe/PruneTopGODAGviz/releases/download/v0.1.1/PruneTopGODAGviz_0.1.1.tar.gz",
  repos = NULL, type = "source"
)
```

The package contains compressed PANTHER annotation tables for Arabidopsis,
maize, rice, soybean, human and mouse, plus example data. No external `data/`
directory is required after installation.

## Quick start

```r
library(PruneTopGODAGviz)

PruneTopGODAGviz::run_app()
```

Select a built-in species annotation, upload a gene or transcript ID list,
choose the ID type and ontology, then run enrichment. The app supports node
dragging, branch deletion, undo/reset, colour control, and PDF/PNG/SVG/TIFF
export.

The full GO hierarchy is downloaded and cached on first use, or can be supplied
locally:

```r
PruneTopGODAGviz::run_app(obo = "path/to/go-basic.obo")
```

## Programmatic use

```r
anno <- read_go_annotation(
  system.file("extdata", "At_GO_info_panther.csv.gz",
              package = "PruneTopGODAGviz"),
  keytype = "gene"
)

genes <- readLines("genes.txt")

ego <- enrich_go(genes, anno, ont = "BP")

dag <- prune_go_dag(
  ego,
  top_n = 10,
  p_cutoff = 0.05,
  collapse_chains = TRUE,
  cluster_jaccard = 0.85
)

p <- plot_go_dag(dag)
save_go_dag(dag, "BP_pruned", plot = p, dpi = 800)
plot_go_dag_interactive(dag, file = "BP_pruned.html")
```

Run BP, MF and CC together:

```r
run_go_workflow(genes, anno, outdir = "GO_figures",
                ont = c("BP", "MF", "CC"))
```

## Citation

Concept DOI: https://doi.org/10.5281/zenodo.23587324

Version DOI for this release: https://doi.org/10.5281/zenodo.23587325


```r
citation("PruneTopGODAGviz")
```

`CITATION.cff` supports GitHub's **Cite this repository** feature. Please cite
the specific package version/DOI used, plus the Gene Ontology Consortium,
PANTHER and the genome annotation resources listed in
`inst/extdata/DATA_PROVENANCE.md`.

## License

Package code: MIT. Built-in annotation tables retain the terms of their
upstream sources; see `inst/extdata/DATA_PROVENANCE.md`.

[中文介绍](README_zh.md)
