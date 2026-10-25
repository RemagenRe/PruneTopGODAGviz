# PruneTopGODAGviz 中文说明

从 GO 富集结果生成适合发表使用的 Gene Ontology DAG，并且只剪掉“展示出来”的
分支，不改变富集统计结果。

`PruneTopGODAGviz` 可以读取 transcript 或 gene 级别的 GO 注释表，使用包内自带的
超几何检验完成 GO over-representation analysis，再绘制 agriGO/topGO 风格的
GO DAG。图形可以剪枝、重新排列、折叠和导出，而富集分析的统计值保持不变。

## 设计背景

这个包参考了 [`topGO`](https://bioconductor.org/packages/topGO/) 和
[agriGO v2](https://systemsbiology.cau.edu.cn/agriGOv2/specises_analysis.php?SpeciseID=1&latin=Arabidopsis_thaliana)
展示 GO DAG 的思路。

这些工具生成的图信息完整，但当显著 term 较多时，整张图容易变得很大：节点密集、
连线交叉、标签拥挤，很难直接放进论文或汇报中。`PruneTopGODAGviz` 重点解决
“图太大、term 太多、布局不可控”的问题。

## 主要功能

- 读取 transcript 或 gene 级别的 GO 注释表；
- 包内自包含的超几何检验，不依赖 `clusterProfiler`；
- 支持 BP、MF、CC 三个本体；
- GO DAG 布局、剪枝、单子节点链折叠和相似 term 合并；
- Shiny 交互编辑器：拖动节点、删除分支、撤销、重置、调整颜色；
- 导出 PDF、PNG、SVG 和 TIFF，默认 800 dpi；
- 生成独立交互式 HTML 图；
- 内置六个物种的压缩 PANTHER 注释表和示例数据；
- 富集计算在后台运行，界面在长时间计算时仍可响应。

## 环境要求

- R >= 4.2；
- Bioconductor 包：`GO.db`、`AnnotationDbi`（必须安装）；
- Shiny 界面：`shiny`、`visNetwork`、`colourpicker`、
  `future`、`promises`；
- 导出和交互图：`ragg`、`svglite`、`ggiraph`、`htmlwidgets`、
  `htmltools`。

## 安装

R 包发布在 GitHub Release，不在 CRAN。请**按下面的顺序安装**。

### 第一步：安装必须的 Bioconductor 依赖

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

`GO.db` 和 `AnnotationDbi` 位于 `Imports` 中，缺少它们时包无法安装。  
不要使用 `install.packages("GO.db")`，它不在 CRAN 上。

### 第二步：安装 Shiny、后台计算和导出相关包

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

各包用途：

| 包 | 用途 |
|---|---|
| `shiny`、`visNetwork`、`colourpicker` | 交互式编辑器 |
| `future`、`promises` | 后台富集计算，保持界面响应 |
| `ragg` | PNG、TIFF 导出 |
| `svglite` | SVG 导出 |
| `ggiraph`、`htmlwidgets` | 独立交互式 HTML 图 |
| `htmltools` | Shiny/HTML 支持 |

### 第三步：安装 PruneTopGODAGviz

安装最新版 tarball：

```r
install.packages(
  "https://github.com/RemagenRe/PruneTopGODAGviz/releases/download/v0.1.4/PruneTopGODAGviz_0.1.4.tar.gz",
  repos = NULL,
  type = "source",
  dependencies = TRUE
)
```

也可以使用 `pak`：

```r
install.packages("pak")
pak::pkg_install(
  "RemagenRe/PruneTopGODAGviz@v0.1.4",
  dependencies = TRUE
)
```

### 第四步：检查安装

重启 R/RStudio，然后运行：

```r
library(PruneTopGODAGviz)
packageVersion("PruneTopGODAGviz")
```

应显示 `0.1.4`。

## 快速开始：Shiny 界面

```r
library(PruneTopGODAGviz)
PruneTopGODAGviz::run_app()
```

操作顺序：

1. 选择内置物种注释，或上传自己的 GO 注释表；
2. 上传 gene/transcript ID 列表；
3. 选择 ID 类型和本体（BP、MF 或 CC）；
4. 运行富集；
5. 在图中拖动、删除或折叠 term；
6. 导出图片。

富集计算在后台 worker 中运行，状态区会显示进度；计算期间可以切换标签页或短时间
离开页面。

## GO OBO 文件说明

完整的 `go-basic.obo` **没有随 R 包一起打包**。它大约 30 MB，而且 GO 本体
会持续更新，因此程序采用“本地优先、需要时下载并缓存”的方式。

`go_dag()` 和 `run_app()` 按以下顺序查找 OBO：

1. 明确传入的 `obo` 路径；
2. `getOption("PruneTopGODAGviz.obo")` 设置的路径；
3. 用户缓存目录中的 `go-basic.obo`，目录位置由
   `tools::R_user_dir("PruneTopGODAGviz", "cache")` 给出；
4. 从 `https://purl.obolibrary.org/obo/go/go-basic.obo` 下载。

如果找不到且无法下载，程序会回退到 `GO.db`，并提示：

```text
GO OBO file unavailable; falling back to GO.db.
Relation types will all be labelled 'is_a'.
```

`GO.db` 可以作为备用层级，但它不保存 `part_of` 等关系类型，因此 DAG 的边类型和
部分拓扑可能与 OBO 版本不同。需要准确展示 `part_of` 关系时，建议使用 OBO。

### 使用本地 OBO

```r
PruneTopGODAGviz::run_app(
  obo = "D:/GO/go-basic.obo"
)
```

或者在当前 R 会话中设置选项：

```r
options(PruneTopGODAGviz.obo = "D:/GO/go-basic.obo")
PruneTopGODAGviz::run_app()
```

### 手动下载 OBO

```r
PruneTopGODAGviz::download_go_obo(
  dest = "D:/GO/go-basic.obo"
)
```

clean project 已经包含 `data/go-basic.obo`，它的启动脚本会自动把该文件传给
`run_app()`，因此不需要再次下载。

## 内置数据

包内包含以下六个物种的压缩 PANTHER 注释表：

- Arabidopsis thaliana（拟南芥）
- Zea mays（玉米）
- Oryza sativa（水稻）
- Glycine max（大豆）
- Homo sapiens（人）
- Mus musculus（小鼠）

此外还包含 mini OBO、示例 annotation 和 5,167 个拟南芥示例基因。安装包后不需要
额外的 `data/` 目录。内置数据来源和上游许可说明见：

```text
inst/extdata/DATA_PROVENANCE.md
```

## 编程方式

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

同时运行 BP、MF 和 CC：

```r
run_go_workflow(
  genes,
  anno,
  outdir = "GO_figures",
  ont = c("BP", "MF", "CC")
)
```

生成交互式 HTML：

```r
plot_go_dag_interactive(pruned, file = "BP_pruned.html")
```

## 导出格式

- PDF：矢量图，优先用于期刊投稿；
- PNG：使用 `ragg`；
- SVG：使用 `svglite`；
- TIFF：使用 `ragg`；
- 默认导出分辨率为 800 dpi。

剪枝只改变图形中显示的节点，不会改变富集统计值；被剪掉的 term 仍保留在结果表
和对象中，可以继续用于导出或后续分析。

## 引用

Concept DOI（代表整个软件，适合一般引用）：

```text
https://doi.org/10.5281/zenodo.23587324
```

已归档的 v0.1.0 Version DOI：

```text
https://doi.org/10.5281/zenodo.23587325
```

建议在论文中引用 Concept DOI，并在 Methods 中注明实际使用的版本号。包内也提供
CITATION 文件：

```r
citation("PruneTopGODAGviz")
```

`CITATION.cff` 支持 GitHub 的 **Cite this repository** 功能。除了软件本身，
还建议引用 Gene Ontology Consortium、PANTHER 以及
`inst/extdata/DATA_PROVENANCE.md` 中列出的基因组注释资源。

## 许可

包代码使用 MIT 许可。内置注释数据保留其上游来源的许可条款，具体见
`inst/extdata/DATA_PROVENANCE.md`。

## 链接

- GitHub 仓库：<https://github.com/RemagenRe/PruneTopGODAGviz>
- 最新 Release：<https://github.com/RemagenRe/PruneTopGODAGviz/releases/tag/v0.1.4>
- 问题反馈：<https://github.com/RemagenRe/PruneTopGODAGviz/issues>
- English README：[README.md](README.md)
