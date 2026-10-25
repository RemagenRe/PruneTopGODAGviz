# PruneTopGODAGviz

Publication-quality Gene Ontology DAGs from over-representation results, with
display pruning that leaves the underlying statistics unchanged.

`PruneTopGODAGviz` reads transcript- or gene-level GO annotation tables, runs a
self-contained hypergeometric over-representation test, and draws an
agriGO/topGO-style directed acyclic graph (DAG). The graph can be pruned,
rearranged, simplified and exported without changing the enrichment result.

## What it is for

TopGO-style GO DAGs are informative, but a DAG with many significant terms can
be too large for a paper figure: nodes overlap, labels collide, and branches
cross. `PruneTopGODAGviz` keeps the same GO hierarchy and statistics while
giving you control over what is displayed.

Typical uses:

- select the most informative enriched terms for a compact figure;
- remove branches that are not relevant to the story;
- collapse long single-child chains;
- merge highly similar sibling terms;
- adjust fonts, spacing, colours and edge styles;
- export a publication-ready figure.

## Features

- Reads transcript- or gene-level GO annotation tables.
- Self-contained hypergeometric ORA; `clusterProfiler` is not required.
- Supports BP, MF and CC ontologies.
- Deterministic DAG layout, branch pruning, chain collapsing and Jaccard-based
  merging of similar terms.
- Shiny editor with node dragging, branch deletion, undo/reset and colour
  controls.
- PDF, PNG, SVG and TIFF export at 800 dpi.
- Standalone interactive HTML figures.
- Compressed PANTHER annotation tables for six species are included in the
  package, together with example data.
- Background enrichment in the Shiny app, so the interface stays responsive
  during long calculations.

## Requirements

- R >= 4.2.
- Bioconductor packages `GO.db` and `AnnotationDbi` (mandatory).
- R packages for the Shiny app: `shiny`, `visNetwork`, `colourpicker`,
  `future` and `promises`.
- Export and interactive-output packages: `ragg`, `svglite`, `ggiraph`,
  `htmlwidgets` and `htmltools`.

## Installation

The package is distributed through GitHub Releases, not CRAN. Install the
dependencies **in the order below**.

### 1. Install the mandatory Bioconductor dependencies

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(
  c("GO.db", "AnnotationDbi"),
  update = FALSE,
  ask = FALSE
)
```

`GO.db` and `AnnotationDbi` are in `Imports`, so the package will not install
without them. Do not try `install.packages("GO.db")`; use `BiocManager`.

### 2. Install the Shiny and export packages

```r
install.packages(
  c(
    "shiny", "visNetwork", "colourpicker",
    "future", "promises",
    "ragg", "svglite",
    "ggiraph", "htmlwidgets", "htmltools"
  ),
  dependencies = TRUE
)
```

What each group is for:

| Package(s) | Needed for |
|---|---|
| `shiny`, `visNetwork`, `colourpicker` | The interactive editor |
| `future`, `promises` | Background enrichment and a responsive UI |
| `ragg` | PNG/TIFF export |
| `svglite` | SVG export |
| `ggiraph`, `htmlwidgets` | Standalone interactive HTML output |
| `htmltools` | Shiny/HTML support |

### 3. Install PruneTopGODAGviz

Install the latest release tarball:

```r
install.packages(
  "https://github.com/RemagenRe/PruneTopGODAGviz/releases/download/v0.1.4/PruneTopGODAGviz_0.1.4.tar.gz",
  repos = NULL,
  type = "source",
  dependencies = TRUE
)
```

Alternatively, with `pak`:

```r
install.packages("pak")
pak::pkg_install(
  "RemagenRe/PruneTopGODAGviz@v0.1.4",
  dependencies = TRUE
)
```

### 4. Verify the installation

Restart R/RStudio, then run:

```r
library(PruneTopGODAGviz)
packageVersion("PruneTopGODAGviz")
```

The expected version is `0.1.4`.

## Quick start: Shiny editor

```r
library(PruneTopGODAGviz)
PruneTopGODAGviz::run_app()
```

In the editor:

1. Select one of the built-in species annotations or upload your own table.
2. Upload a gene or transcript ID list.
3. Choose the ID type and ontology.
4. Run enrichment.
5. Drag, delete or collapse terms in the graph.
6. Export the figure.

The app reports progress in the status area and runs enrichment in a background
worker, so you can switch tabs or leave the page briefly while a calculation is
running.

## GO OBO file

The complete `go-basic.obo` file is **not bundled** in the R package. It is
about 30 MB, changes over time, and is therefore downloaded and cached on first
use.

`go_dag()` and `run_app()` resolve the OBO file in this order:

1. an explicitly supplied `obo` path;
2. `getOption("PruneTopGODAGviz.obo")`;
3. the user cache directory created by
   `tools::R_user_dir("PruneTopGODAGviz", "cache")`;
4. a download from `https://purl.obolibrary.org/obo/go/go-basic.obo`.

If the file is unavailable, the package falls back to `GO.db` and warns that
all relations are labelled `is_a`. This fallback works, but it does not preserve
`part_of` relation types and may produce a slightly different DAG.

To use a local file:

```r
PruneTopGODAGviz::run_app(
  obo = "D:/GO/go-basic.obo"
)
```

Or set a session-wide option:

```r
options(PruneTopGODAGviz.obo = "D:/GO/go-basic.obo")
PruneTopGODAGviz::run_app()
```

To download the OBO file explicitly:

```r
PruneTopGODAGviz::download_go_obo(dest = "D:/GO/go-basic.obo")
```

The clean project distributed with this package already includes
`data/go-basic.obo`; its launcher passes that file to `run_app()`.

## Built-in data

The package contains compressed PANTHER annotation tables for:

- Arabidopsis thaliana
- Zea mays
- Oryza sativa
- Glycine max
- Homo sapiens
- Mus musculus

It also includes a mini OBO file, example annotations and 5,167 Arabidopsis
example genes. No external `data/` directory is required after installation.
Source and licence information for the built-in tables is in
[`inst/extdata/DATA_PROVENANCE.md`](inst/extdata/DATA_PROVENANCE.md).

## Programmatic use

```r
library(PruneTopGODAGviz)

anno <- read_go_annotation(
  system.file(
    "extdata", "At_GO_info_panther.csv.gz",
    package = "PruneTopGODAGviz"
  ),
  keytype = "gene"
)

genes <- readLines("genes.txt")

ego <- enrich_go(genes, anno, ont = "BP")

dag <- go_dag("BP")
pruned <- prune_go_dag(
  ego,
  dag = dag,
  top_n = 10,
  p_cutoff = 0.05,
  collapse_chains = TRUE,
  cluster_jaccard = 0.85
)

p <- plot_go_dag(pruned)
save_go_dag(pruned, "BP_pruned", plot = p, dpi = 800)
```

Run BP, MF and CC together:

```r
run_go_workflow(
  genes,
  anno,
  outdir = "GO_figures",
  ont = c("BP", "MF", "CC")
)
```

Create an interactive HTML figure:

```r
plot_go_dag_interactive(pruned, file = "BP_pruned.html")
```

## Notes on export

- PDF is vector output and is the preferred format for journal submission.
- PNG uses `ragg`.
- SVG uses `svglite`.
- TIFF uses `ragg`.
- The default export resolution is 800 dpi.
- Pruning changes only the displayed DAG; the enrichment statistics and the
  pruned result table remain available.

## Citation

Concept DOI (all versions):

```text
https://doi.org/10.5281/zenodo.23587324
```

The archived v0.1.0 version DOI is:

```text
https://doi.org/10.5281/zenodo.23587325
```

Use the concept DOI when citing the software in general. Cite the exact
version you used in the Methods section. The package also provides a CITATION
file:

```r
citation("PruneTopGODAGviz")
```

`CITATION.cff` supports GitHub's **Cite this repository** feature. In addition
to the software citation, cite the Gene Ontology Consortium, PANTHER and the
genome annotation resources listed in
[`inst/extdata/DATA_PROVENANCE.md`](inst/extdata/DATA_PROVENANCE.md).

## License

Package code is released under the MIT license. Built-in annotation tables
retain the terms of their upstream sources; see
[`inst/extdata/DATA_PROVENANCE.md`](inst/extdata/DATA_PROVENANCE.md).

## Links

- Repository: <https://github.com/RemagenRe/PruneTopGODAGviz>
- Latest release: <https://github.com/RemagenRe/PruneTopGODAGviz/releases/tag/v0.1.4>
- Issue tracker: <https://github.com/RemagenRe/PruneTopGODAGviz/issues>
- Chinese README: [README_zh.md](README_zh.md)
