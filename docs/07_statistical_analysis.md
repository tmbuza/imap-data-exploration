# (PART) STATISTICAL ANALYSIS

# Statistical Significance Analysis

This section delves into statistical analysis methods tailored specifically for microbiome data. Statistical analysis plays a crucial role in understanding the complex relationships and patterns within microbiome datasets, helping researchers uncover significant findings and insights into microbial community dynamics, composition, and responses to various environmental factors or treatments.





```r
# Load required packages
library(phyloseq)
library(tidyverse)

cat("\nSaved RData objects\n\n")

Saved RData objects
load("data/dataframe_objects.rda", verbose = T)
Loading objects:
  df_GlobalPatterns
  df_dietswap
  df_caporaso
  df_kostic_crc
load("data/phyloseq_objects.rda", verbose = T)
Loading objects:
  ps_GlobalPatterns
  ps_dietswap
  ps_caporaso
  ps_kostic_crc
load("data/phyloseq_extra_objects.rda", verbose = T)
Loading objects:
  psextra_clr_dietswap
  psextra_id_dietswap
  psextra_log10p_dietswap
load("data/ps_transformed.rda", verbose = T)
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
load("data/bray_distances.rda", verbose = T)
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
load("data/psextra_distances.rda", verbose = T)
Loading objects:
  psextra_clr_asin_bray_dist
  psextra_id_asin_bray_dist
  psextra_log10p_asin_bray_dist
load("data/reduced_dimension.rda", verbose = T)
Loading objects:
  pca_asin_bray_metrics
  mds_asin_bray_metrics
  pcoa_asin_bray_metrics
  tsne_asin_bray_metrics
load("data/phyloseq_raw_rel_psextra_df_objects.rda", verbose = T)
Loading objects:
  ps_raw
  ps_rel
  psextra_raw
  psextra_rel
  ps_df
```

## PERMANOVA
PERMANOVA (Permutational Multivariate Analysis of Variance) is a statistical test used to assess the significance of differences between groups of microbial communities. Particularly suited for analyzing multivariate data like microbiome composition, PERMANOVA offers valuable insights into how different experimental conditions or treatments impact microbial community structure.


```r
library(microViz) 

bray_perm <- psextra_log10p_asin_bray_dist %>%
  dist_permanova(
    seed = 1234, # for set.seed to ensure reproducibility of random process
    n_processes = 1, n_perms = 99, # you should use at least 999!
    variables = "bmi_group"
  )

perm_get(bray_perm) %>% as.data.frame()
           Df  SumOfSqs         R2        F Pr(>F)
bmi_group   2 0.1233611 0.04544245 5.212833   0.01
Residual  219 2.5913059 0.95455755       NA     NA
Total     221 2.7146670 1.00000000       NA     NA


info_get(bray_perm)
psExtra info:
tax_agg = "Genus" tax_trans = "log10p" dist_method = "bray" 
```


