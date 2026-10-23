# Built-in annotation data provenance

These compressed tables are provided as convenience examples for
`PruneTopGODAGviz`. They contain public gene identifiers and Gene Ontology
identifiers only; no personal or patient data are included.

| File | Source annotation | Genome/transcript models | Enrichment supplement |
|---|---|---|---|
| `At_GO_info_panther.csv.gz` | TAIR/Ensembl/GO/orgDB-derived `At_GO_info.csv` | TAIR10.60 GTF | PANTHER gene-info API (organism 3702) |
| `Maize_GO_info_panther.csv.gz` | MaizeGDB/Ensembl/agriGO-Pannzer table | Zm-B73-REFERENCE-NAM-5.0.60 GTF | PANTHER gene-info API (organism 4577) |
| `Rice_GO_info_panther.csv.gz` | RAP-DB/Ensembl/GO table | IRGSP-1.0.60 GTF | PANTHER gene-info API (organism 39947) |
| `Soybean_GO_info_panther.csv.gz` | PlantRegMap/GO table | Glycine_max_v2.1.60 GTF | PANTHER gene-info API (organism 3847) |
| `Homo_sapiens_GO_info_panther.csv.gz` | PANTHER gene-info API | GRCh38.115 GTF | PANTHER (organism 9606) |
| `Mus_musculus_GO_info_panther.csv.gz` | PANTHER gene-info API | GRCm39.115 GTF | PANTHER (organism 10090) |

Processing notes:

- PANTHER annotations were retrieved from the public PANTHER gene-info
  service. Complete BP/MF/CC annotations were used where available.
- Obsolete/nonexistent GO terms were filtered against GO release
  **2026-07-26**.
- Each valid PANTHER gene-GO pair was expanded to every transcript of that
  gene in the corresponding GTF, so multi-isoform genes have the same gene-level
  annotation on all annotated transcripts.
- Tables are gzip-compressed for package size; `read_go_annotation()` reads
  `.csv.gz` directly.

Users redistributing or publishing analyses with these tables should follow the
terms of the upstream TAIR, Ensembl, MaizeGDB, RAP-DB, PlantRegMap, PANTHER and
Gene Ontology resources. The package's MIT license covers the package code, not
the upstream annotation data.
