# (PART) DEMO DATASET {-}

# Microbiome Demo Datasets {#demo-datasets}

In this section, we provide demo datasets for microbiome analysis. These datasets have been preprocessed and transformed in previous sessions, and now we are preparing to analyze them. We begin by importing the R and Phyloseq objects.






## Import preprepared data

```r
# Load required packages
library(phyloseq)
library(tidyverse)

cat("\nSaved RData objects\n\n")
```

```

Saved RData objects
```

```r
load("data/ps_transformed.rda", verbose = T)
```

```
Loading objects:
  ps_asin
  ps_identity
  ps_compositional
  ps_z_otu
  ps_z_sample
  ps_log10
  ps_log10p
  ps_clr
  ps_shift
  ps_scale
```

```r
load("data/bray_distances.rda", verbose = T)
```

```
Loading objects:
  ps_asin_bray_dist
  ps_compositional_bray_dist
  ps_z_otu_bray_dist
  ps_z_sample_bray_dist
  ps_log10_bray_dist
  ps_log10p_bray_dist
  ps_clr_bray_dist
  ps_shift_bray_dist
  ps_scale_bray_dist
```

```r
load("data/reduced_dimension.rda", verbose = T)
```

```
Loading objects:
  pca_asin_bray_metrics
  mds_asin_bray_metrics
  pcoa_asin_bray_metrics
  tsne_asin_bray_metrics
```

```r
load("data/phyloseq_raw_rel_psextra_df_objects.rda", verbose = T)
```

```
Loading objects:
  ps_raw
  ps_rel
  psextra_raw
  psextra_rel
  ps_df
```


## Demo data for QIIME2 analysis

The demo data for QIIME2 analysis is downloaded and stored in the `data/qiime2` directory. You can use these files to practice QIIME2 workflows and analysis techniques.

The data was sourced from the QIIME2 tutorials available at [https://docs.qiime2.org/2024.4/data/tutorials/moving-pictures/](https://docs.qiime2.org/2018.4/data/tutorials/moving-pictures/).

It contains the following files:

- Sample metadata
- Feature table
- Taxonomy
- Shannon vector
- Rooted tree

This data is essential for leveraging downstream analysis aimed at understanding data analyzed using the QIIME2 pipeline.


```R
if (!dir.exists("data")){dir.create("data")}
if (!dir.exists("data/qiime2")){dir.create("data/qiime2")}

download.file("https://docs.qiime2.org/2024.2/data/tutorials/moving-pictures/table.qza", "data/qiime2/feature_table.qza")
download.file("https://data.qiime2.org/2024.2/tutorials/moving-pictures/sample_metadata.tsv", "data/qiime2/sample_metadata.tsv")
download.file("https://docs.qiime2.org/2024.2/data/tutorials/moving-pictures/taxonomy.qza", "data/qiime2/taxonomy.qza")
download.file("https://docs.qiime2.org/2024.2/data/tutorials/moving-pictures/rooted-tree.qza", "data/qiime2/rooted_tree.qza")
download.file("https://docs.qiime2.org/2024.2/data/tutorials/moving-pictures/core-metrics-results/shannon_vector.qza", "data/qiime2/shannon_vector.qza")

```

