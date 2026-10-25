# PruneTopGODAGviz 中文介绍

## 灵感来源

这个 R 包的灵感来自两套常用的 GO 富集与可视化工具：

- [`topGO`](https://bioconductor.org/packages/topGO/)：在 R 中做 GO 富集分析，并输出 GO DAG 图。
- [agriGO v2 物种分析页面](https://systemsbiology.cau.edu.cn/agriGOv2/specises_analysis.php?SpeciseID=1&latin=Arabidopsis_thaliana)：提供拟南芥等物种的 GO 富集和 topGO 风格 DAG 展示。

这些工具生成的 topGO 图信息完整，但当显著 term 较多时，整张图会变得很大：节点密集、连线交叉、文字拥挤，难以直接放进论文或汇报中做紧凑展示。

## 为什么做这个包

`PruneTopGODAGviz` 重点解决“图太大、term 太多、布局不可控”的问题：

- 可以自由拖动节点，按自己的叙事逻辑重新排列；
- 可以删除不想突出展示的分支，只保留关键 term；
- 可以折叠单子节点链，减少冗余层级；
- 可以合并基因集高度相似的 sibling term；
- 可以统一调整字体、节点大小、颜色和边样式；
- 可以直接导出 PDF、PNG、SVG、TIFF，方便论文排版。

它不改变富集分析的统计值，只改变“如何展示”。这样既能保留 topGO/agriGO 风格的 GO DAG，又能得到更紧凑、更适合发表的图。

## 主要功能

- 读取 transcript 或 gene 级别的 GO 注释表；
- 包内自包含的超几何检验，不依赖 `clusterProfiler`；
- 支持 BP、MF、CC 三个本体；
- GO DAG 布局、剪枝、链折叠、Jaccard 合并；
- Shiny 交互式编辑器：拖动、双击删除分支、撤销、重置、颜色控制；
- 发表级 PDF/PNG/SVG/TIFF 导出；
- 可生成交互式 HTML 图；
- 内置压缩注释表和示例数据，安装后即可试用。

## 内置注释数据

安装包内已包含压缩的 PANTHER 补充注释表：

- Arabidopsis thaliana
- Zea mays
- Oryza sativa
- Glycine max
- Homo sapiens
- Mus musculus

同时包含 mini GO OBO、示例 annotation 和 5,167 个拟南芥示例基因。数据来源和上游许可说明见：

```text
inst/extdata/DATA_PROVENANCE.md
```

### 依赖安装

`GO.db` 和 `AnnotationDbi` 是必须依赖。请先从 Bioconductor 安装，再安装
PruneTopGODAGviz：

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("GO.db", "AnnotationDbi"))
```

直接安装 tarball 时，`repos = NULL, type = "source"` 仍不会自动安装全部
Suggests。运行 `run_app()` 前请显式安装 Shiny 和导出相关包：

```r
install.packages(c(
  "shiny", "visNetwork", "colourpicker",
  "ragg", "svglite", "ggiraph", "htmlwidgets", "htmltools"
))
```

`run_app()` 最小依赖是：`shiny`、`visNetwork`、`colourpicker`。
`ragg`/`svglite` 用于 PNG/SVG 导出；`ggiraph`/`htmlwidgets` 用于交互 HTML。
`GO.db`/`AnnotationDbi` 是必须依赖，已在上方通过 BiocManager 安装。

## 安装

要求 R >= 4.2。

```r
install.packages("pak")
pak::pkg_install("RemagenRe/PruneTopGODAGviz@v0.1.1",
                 dependencies = TRUE)
```

或安装 GitHub Release 中的 tarball：

```r
install.packages(
  "https://github.com/RemagenRe/PruneTopGODAGviz/releases/download/v0.1.1/PruneTopGODAGviz_0.1.1.tar.gz",
  repos = NULL, type = "source"
)
```

## 快速使用

```r
library(PruneTopGODAGviz)

# 打开交互式编辑器
PruneTopGODAGviz::run_app()
```

在应用中选择内置物种注释，上传自己的 gene/transcript ID，选择本体后运行富集；随后可以拖动节点、删除分支、调整颜色并导出图片。

编程方式：

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
```

## 与 topGO/agriGO 图的关系

本包借鉴的是它们“GO DAG + 富集结果”的展示思路，但不是复制其代码或界面。它更适合需要手动组织图形、突出少数关键分支、并把图压缩到论文版式的场景。

## 引用与链接

- 软件作者：Yuzhou Liu
- 仓库：<https://github.com/RemagenRe/PruneTopGODAGviz>
- Concept DOI：<https://doi.org/10.5281/zenodo.23587324>
- v0.1.1 Version DOI：<https://doi.org/10.5281/zenodo.23587325>

```text
Liu, Yuzhou. (2026). PruneTopGODAGviz (Version 0.1.0)
[Computer software]. Zenodo.
https://doi.org/10.5281/zenodo.23587325
```

包代码使用 MIT 许可；内置注释数据的上游条款见 `DATA_PROVENANCE.md`。
