# (PART) START EXPLORING {-}

# Exploration of Visualization Types




## Import data

```r
# Load required packages
library(phyloseq)
library(tidyverse)

cat("\nSaved RData objects\n\n")

Saved RData objects
load("data/external_ps_objects.rda", verbose = T)
Loading objects:
  df_GlobalPatterns
  df_dietswap
  df_caporaso
  df_kostic_crc
  ps_GlobalPatterns
  ps_dietswap
  ps_caporaso
  ps_kostic_crc
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
load("data/phyloseq_extra_objects.rda", verbose = T)
Loading objects:
  psextra_clr_dietswap
  psextra_id_dietswap
  psextra_log10p_dietswap
load("data/phyloseq_raw_rel_psextra_df_objects.rda", verbose = T)
Loading objects:
  ps_raw
  ps_rel
  psextra_raw
  psextra_rel
  ps_df
```


## Major visualization R colors

In R, there are several built-in palettes that we can use for color schemes in plots. Some commonly used palettes include:


```
[1] "#66C2A5" "#FC8D62" "#8DA0CB" "#E78AC3" "#A6D854" "#FFD92F" "#E5C494"
[8] "#B3B3B3"
[1] "#B3B3B3" "#E5C494" "#FFD92F" "#A6D854" "#E78AC3" "#8DA0CB" "#FC8D62"
[8] "#66C2A5"
[1] "#440154FF" "#46337EFF" "#365C8DFF" "#277F8EFF" "#1FA187FF" "#4AC16DFF"
[7] "#9FDA3AFF" "#FDE725FF"
[1] "#FDE725FF" "#9FDA3AFF" "#4AC16DFF" "#1FA187FF" "#277F8EFF" "#365C8DFF"
[7] "#46337EFF" "#440154FF"
[1] "#FF0000" "#FFBF00" "#80FF00" "#00FF40" "#00FFFF" "#0040FF" "#8000FF"
[8] "#FF00BF"
[1] "#FF00BF" "#8000FF" "#0040FF" "#00FFFF" "#00FF40" "#80FF00" "#FFBF00"
[8] "#FF0000"
```

> - **viridis**: A perceptually uniform and colorblind-friendly palette.
> - **magma**: A palette with a dark-to-light color scheme.
> - **inferno**: A palette with a light-to-dark color scheme.
> - **plasma**: A palette with a dark-to-light color scheme.
> - **cool**: A palette with cool colors.
> - **hot**: A palette with hot colors.
> - **terrain.colors**: A palette with colors resembling a terrain map.
> - **rainbow**: A palette with colors of the rainbow.
> - **heat.colors**: A palette with colors ranging from dark red to yellow.
> 
> - In ggpubr: "npg", "aaas", "lancet", "jco", "ucscgb", "uchicago", "simpsons" and "rickandmorty".


## Major visualization techniques {#visual-types}

Below are some of the key visualization techniques used in microbiome research, along with their descriptions and the corresponding tools in R.


| Visual Type | Description |
|-------|-------------|
| Barplots | Display the relative abundances of different taxa across groups. |
| Heatmaps | Represent the abundance or presence/absence of taxa across samples. |
| Scatter plots | Useful for visualizing relationships between numerical variables |
| Box plots | Summarize the distribution of a variable |
| PCA plots | Dimensionality reduction technique for visualizing similarities or dissimilarities between samples based on their microbial composition |
| Alpha diversity plots | Measure the diversity within a sample, e.g. rarefaction plot |
| Beta diversity plots | Measure the dissimilarity between samples, e.g. PCoA ordination |
| Line plot | Visualize changes in the abundance of specific taxa over time or across different conditions. |
| Network plots | Depict interactions or associations between taxa|
| Volcano plots | Identify statistically significant differences in abundance between groups |
| Correlation plots | Visualize correlations between taxa or between taxa and metadata variables |
| UpSet plots | Display intersections of sets and their size in a matrix layout |
| Venn diagrams | Show overlap between taxa or groups |
| Differential abundance plots | Visualize differences in abundance between groups while controlling for confounding factors |
| Indicator species analysis plots | Identify taxa associated with specific groups or conditions |



# (PART) BAR CHARTS {-}

# Bar Plots 

Bar plots are commonly used to visualize the relative abundance of taxa or other features within different groups or conditions in microbiome data. In this section, we will demonstrate different methods for plotting bar plots.

## Explore categorical variables
The freq() function is in the funModeling package provides a concise summary of the frequency distribution of a categorical variable in a dataset. It displays the count and percentage of each category, allowing users to quickly understand the distribution and prevalence of different categories within the variable. This function is useful for exploratory data analysis and understanding the composition of categorical variables in the dataset.


```r
library(phyloseq)
library(dplyr)
library(microbiome)
library(funModeling)

df <-ps_df %>% 
  select(-sample_id) %>% 
  rename_all(tolower)

freq(df, 
     input = c("sex", "nationality", "grp", "bmi"), 
     plot = TRUE,
     na.rm = FALSE, 
     # path_out="figures"
     )
```

<img src="./figures/prepare_df-1.png" width="672" />

```
     sex frequency percentage cumulative_perc
1   male     49604      54.09           54.09
2 female     42096      45.91          100.00
```

<img src="./figures/prepare_df-2.png" width="672" />

```
  nationality frequency percentage cumulative_perc
1         AAM     50468      55.04           55.04
2         AFR     41232      44.96          100.00
```

<img src="./figures/prepare_df-3.png" width="672" />

```
  grp frequency percentage cumulative_perc
1  ED     31080      33.89           33.89
2  HE     31068      33.88           67.77
3  DI     29552      32.23          100.00
```

<img src="./figures/prepare_df-4.png" width="672" />

```
         bmi frequency percentage cumulative_perc
1      obese     37080      40.44           40.44
2 overweight     31320      34.15           74.59
3       lean     23300      25.41          100.00
[1] "Variables processed: sex, nationality, grp, bmi"
```


## Comparative results from Mothur and QIIME2


```r
# scripts/create_barplots.R
library(readr)
library(dplyr)
library(ggplot2)
library(svglite)

read_csv("../imap-data-preparation/data/mothur/mothur_composite.csv", show_col_types = FALSE) %>% 
  mutate(pipeline="MOTHUR", .before=2) %>% 
  select(-OTU, -3) %>% 
  ggplot(aes(x=Genus, y=rel_abund, fill=pipeline)) +
  facet_grid(~pipeline) +
  geom_col() +
  coord_flip() +
  labs(y="Relative Abundance (%)") +
  theme(legend.position="none")
```

<img src="./figures/unnamed-chunk-4-1.png" width="576" />

```r


read_csv("../imap-data-preparation/data/qiime2/qiime2_composite.csv", show_col_types = FALSE) %>% 
  mutate(pipeline="QIIME2", .before=2) %>% 
  select(-feature, -c(3:14)) %>% 
  ggplot(aes(x=Genus, y=rel_abund, fill=pipeline)) +
  facet_grid(~pipeline) +
  geom_col() +
  coord_flip() +
  labs(y="Relative Abundance (%)") + 
  theme(legend.position="none")
```

<img src="./figures/unnamed-chunk-4-2.png" width="576" />

```r



ggsave(file="figures/taxon_barplot.png", width=10, height=10)
ggsave(file="figures/taxon_barplot.svg", width=10, height=10)
```


## Compute sequence count per sample
Here we only show the head and tail count data


```r
seqcount_per_sample <- ps_df %>%
  group_by(sample_id) %>% 
  summarise(nseqs = sum(count), .groups = "drop") %>% 
  arrange(-nseqs)

head_tail <- rbind(head(seqcount_per_sample, 5), tail(seqcount_per_sample, 5))

max_y <- max(head_tail$nseqs)

head_tail %>% 
  mutate(sample_id = factor(sample_id),
         sample_id = fct_reorder(sample_id, nseqs, .desc = F),
         sample_id = fct_shift(sample_id, n = 0)) %>%
  ggplot(aes(x = reorder(sample_id, nseqs), y = nseqs)) +
  geom_col(position = position_stack(), fill = "steelblue") +
  labs(x = "sample ID", y = "Number of sequences", subtitle = "Phylum: Top and bottom sequence count \nwith no data labels") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
```

<img src="./figures/unnamed-chunk-5-1.png" width="672" />

```r

ggsave("figures/basic_barplot.png", width=5, height=4)


#---------------------------------------


head_tail %>% 
  mutate(sample_id = factor(sample_id),
         sample_id = fct_reorder(sample_id, nseqs, .desc = F),
         sample_id = fct_shift(sample_id, n = 0)) %>%
  ggplot(aes(x = reorder(sample_id, nseqs), y = nseqs)) +
  geom_col(position = position_stack(), fill = "steelblue") +
  labs(x = "sample ID", y = "Number of sequences", subtitle = "Phylum: Top and bottom sequence count \nwith data labels") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
  coord_cartesian(ylim = c(0, max_y)) +
  geom_text(aes(label = nseqs), vjust = -0.3, color = "#AAAAAA")
```

<img src="./figures/unnamed-chunk-5-2.png" width="672" />

```r

ggsave("figures/barplot_w_labels.png", width=5, height=4)
```




## Relative abundance of taxa

The `comp_barplot` function in the `microViz` R package is designed to create comparative bar plots for visualizing microbiome data. This function allows users to compare the relative abundance of taxa or other features across different groups or conditions. It supports customization options such as color palette selection and plot title specification to enhance the appearance of the bar plot. 


```r
library(microViz)

psextra_raw %>%
  comp_barplot(
    tax_level = "Genus", n_taxa = 15, other_name = "Other",
    taxon_renamer = function(x) stringr::str_remove(x, " [ae]t rel."),
    palette = distinct_palette(n = 15, add = "grey90"),
    merge_other = FALSE, bar_outline_colour = "darkgrey"
  ) +
  coord_flip() +
  facet_wrap("nationality", nrow = 1, scales = "free") +
  labs(x = NULL, y = NULL) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
```

<img src="./figures/barplotmicroViz-1.png" width="960" />


## Bar plot on ps_caporaso dataset


```r
library(microViz)

ps_caporaso %>%
  tax_fix() %>% 
  comp_barplot(
    tax_level = "Genus", n_taxa = 15, other_name = "Other",
    taxon_renamer = function(x) stringr::str_remove(x, " [ae]t rel."),
    palette = distinct_palette(n = 15, add = "grey90"),
    merge_other = FALSE, bar_outline_colour = "darkgrey"
  ) +
  coord_flip() +
  # facet_wrap("SampleType", nrow = 1, scales = "free") +
  labs(x = NULL, y = NULL) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
```

<img src="./figures/unnamed-chunk-6-1.png" width="960" />


## Bar plots for Taxa relative abundane


```r
set.seed(1234)
library(tidyverse)
library(readxl)
library(ggtext)
library(RColorBrewer)



otu_rel_abund <- ps_df %>%
  mutate(nationality = factor(nationality, 
                      levels = c("AAM", "AFR"),
                      labels = c("African American", "African")),
         bmi = factor(bmi,
                      levels = c("lean", "overweight", "obese"),
                      labels = c("Lean", "Overweight", "Obese")),
         sex = factor(sex,
                      levels =c("female", "male"),
                      labels = c("Female", "Male")))

save(otu_rel_abund, file = "data/otu_rel_abund.rda")

phylum_rel_abund <- otu_rel_abund %>%
  filter(level=="phylum",
         !grepl(".*unassigned.*|.*nclassified.*|.*ncultured.*",taxon)) %>%
  group_by(bmi, sample_id, taxon) %>%
  summarise(rel_abund = sum(rel_abund), .groups="drop") %>%
  group_by(bmi, taxon) %>%
  summarise(mean_rel_abund = 100*mean(rel_abund), .groups="drop") %>%
  mutate(taxon = str_replace(taxon, "(.*)_unclassified", "Unclassified<br>*\\1*"),
         taxon = str_replace(taxon, "^$", "*\\1*"),
         taxon = str_replace(taxon, "_", " "))

## Pool low abundance to "Others"
taxon_pool <- phylum_rel_abund %>%
  group_by(taxon) %>%
  summarise(pool = max(mean_rel_abund) < 1,
            mean = mean(mean_rel_abund), .groups="drop")


inner_join(phylum_rel_abund, taxon_pool, by="taxon") %>%
  mutate(taxon = if_else(pool, "Other", taxon)) %>%
  group_by(bmi, taxon) %>%
  summarise(mean_rel_abund = sum(mean_rel_abund),
            mean = min(mean), .groups="drop") %>%
  mutate(taxon = factor(taxon),
         taxon = fct_reorder(taxon, mean, .desc=TRUE),
         taxon = fct_shift(taxon, n=1)) %>%
  ggplot(aes(x=bmi, y=mean_rel_abund, fill=taxon)) +
  geom_col(position = position_stack()) +
  scale_fill_discrete(name=NULL) +
  scale_x_discrete(breaks = c("lean", "overweight", "obese"), labels = c("Lean", "Overweight", "Obese")) +
  scale_fill_manual(name=NULL,
                    breaks=c("*Actinobacteria*","*Bacteroidetes*", "*Cyanobacteria*", "*Firmicutes*","*Proteobacteria*", "Other"),
                    labels = c("Actinobacteria","Bacteroidetes", "Cyanobacteria", "Firmicutes", "Proteobacteria", "Other"),
                    values = c(brewer.pal(5, "Dark2"), "gray")) +
  labs(x = NULL, y = "Mean Relative Abundance", subtitle = "Phyla stacked bars filled by taxon") +
  theme_classic() +
  theme(axis.text.x = element_markdown(),
        legend.text = element_markdown(),
        legend.key.size = unit(10, "pt"))
```

<img src="./figures/unnamed-chunk-7-1.png" width="672" />

```r

ggsave("figures/stacked_phyla.png", width=5, height=4)

```



```r
#---------------------------------------

taxon_rel_abund <- otu_rel_abund %>%
  filter(level=="genus",
         !grepl(".*unassigned.*|.*nclassified.*|.*ncultured.*",taxon)) %>%
  group_by(bmi, sample_id, taxon) %>%
  summarise(rel_abund = sum(rel_abund), .groups="drop") %>%
  group_by(bmi, taxon) %>%
  summarise(mean_rel_abund = 100*mean(rel_abund), .groups="drop") %>%
  mutate(taxon = str_replace(taxon, "(.*)_unclassified", "Unclassified<br>*\\1*"),
         taxon = str_replace(taxon, "^$", "*\\1*"),
         taxon = str_replace(taxon, "_", " "))

#---------------------------------------

## Pool low abundance to "Others"

taxon_pool <- taxon_rel_abund %>%
  group_by(taxon) %>%
  summarise(pool = max(mean_rel_abund) < 2,
            mean = mean(mean_rel_abund), .groups="drop")


inner_join(taxon_rel_abund, taxon_pool, by="taxon") %>%
  mutate(taxon = if_else(pool, "Other", taxon)) %>%
  group_by(bmi, taxon) %>%
  summarise(mean_rel_abund = sum(mean_rel_abund),
            mean = min(mean), .groups="drop") %>%
  mutate(taxon = factor(taxon),
         taxon = fct_reorder(taxon, mean, .desc=TRUE),
         taxon = fct_shift(taxon, n=1)) %>%
  ggplot(aes(x=bmi, y=mean_rel_abund, fill=taxon)) +
  geom_col(position = position_stack()) +
  scale_fill_discrete(name=NULL) +
  scale_x_discrete(breaks = c("lean", "overweight", "obese"), labels = c("Lean", "Overweight", "Obese")) +
  labs(x = NULL, y = "Mean Relative Abundance", subtitle = "Taxa stacked bars filled by taxon") +
  theme_classic() +
  theme(axis.text.x = element_markdown(),
        legend.text = element_markdown(),
        legend.key.size = unit(10, "pt"))
```

<img src="./figures/unnamed-chunk-8-1.png" width="672" />



```r
ggsave("figures/stacked_taxa.png", width=5, height=4)

#---------------------------------------

phylum_rel_abund <- otu_rel_abund %>%
  filter(level=="phylum",
         !grepl(".*unassigned.*|.*nclassified.*|.*ncultured.*",taxon)) %>%
  group_by(bmi, sample_id, taxon) %>%
  summarise(rel_abund = sum(rel_abund), .groups="drop") %>%
  group_by(bmi, taxon) %>%
  summarise(mean_rel_abund = 100*mean(rel_abund), .groups="drop") %>%
  mutate(taxon = str_replace(taxon, "(.*)_unclassified", "Unclassified<br>*\\1*"),
         taxon = str_replace(taxon, "^$", "*\\1*"),
         taxon = str_replace(taxon, "_", " "))

## Pool low abundance to "Others"
taxon_pool <- phylum_rel_abund %>%
  group_by(taxon) %>%
  summarise(pool = max(mean_rel_abund) < 1,
            mean = mean(mean_rel_abund), .groups="drop")


inner_join(phylum_rel_abund, taxon_pool, by="taxon") %>%
  mutate(taxon = if_else(pool, "Other", taxon)) %>%
  group_by(bmi, taxon) %>%
  summarise(mean_rel_abund = sum(mean_rel_abund),
            mean = min(mean), .groups="drop") %>%
  mutate(taxon = factor(taxon),
         taxon = fct_reorder(taxon, mean, .desc=TRUE),
         taxon = fct_shift(taxon, n=1)) %>%
  ggplot(aes(y=taxon, x = mean_rel_abund, fill= bmi)) +
  geom_col(width=0.8, position = position_dodge()) +
  labs(y = NULL, x = "Mean Relative Abundance", subtitle = "Phyla grouped by BMI category", fill = NULL) +
  theme_classic() +
  theme(axis.text.x = element_markdown(angle = 0, hjust = 1, vjust = 1),
        axis.text.y = element_markdown(),
        legend.text = element_markdown(),
        legend.key.size = unit(12, "pt"),
        panel.background = element_blank(),
        # panel.grid.major.x =  element_line(colour = "lightgray", size = 0.1),
        panel.border = element_blank()) +
  guides(fill = guide_legend(ncol=1)) +
  scale_x_continuous(expand = c(0, 0))
```

<img src="./figures/unnamed-chunk-9-1.png" width="672" />

```r


ggsave("figures/grouped_phyla.png", width=5, height=4)
```



```r

#---------------------------------------

taxon_rel_abund <- otu_rel_abund %>%
  filter(level=="genus",
         !grepl(".*unassigned.*|.*nclassified.*|.*ncultured.*",taxon)) %>%
  group_by(bmi, sample_id, taxon) %>%
  summarise(rel_abund = sum(rel_abund), .groups="drop") %>%
  group_by(bmi, taxon) %>%
  summarise(mean_rel_abund = 100*mean(rel_abund), .groups="drop") %>%
  mutate(taxon = str_replace(taxon, "(.*)_unclassified", "Unclassified<br>*\\1*"),
         taxon = str_replace(taxon, "^$", "*\\1*"),
         taxon = str_replace(taxon, "_", " "))

## Pool low abundance to "Others"

taxon_pool <- taxon_rel_abund %>%
  group_by(taxon) %>%
  summarise(pool = max(mean_rel_abund) < 2,
            mean = mean(mean_rel_abund), .groups="drop")

## Plotting
inner_join(taxon_rel_abund, taxon_pool, by="taxon") %>%
  mutate(taxon = if_else(pool, "Other", taxon)) %>%
  group_by(bmi, taxon) %>%
  summarise(mean_rel_abund = sum(mean_rel_abund),
            mean = min(mean), .groups="drop") %>%
  mutate(taxon = factor(taxon),
         taxon = fct_reorder(taxon, mean, .desc=TRUE),
         taxon = fct_shift(taxon, n=1)) %>%
  ggplot(aes(y=taxon, x = mean_rel_abund, fill= bmi)) +
  geom_col(width=0.8, position = position_dodge()) +
  labs(y = NULL, x = "Mean Relative Abundance", subtitle = "Taxa grouped by BMI category", fill = NULL) +
  theme_classic() +
  theme(axis.text.x = element_markdown(angle = 0, hjust = 1, vjust = 1),
        axis.text.y = element_markdown(),
        legend.text = element_markdown(),
        legend.key.size = unit(12, "pt"),
        panel.background = element_blank(),
        # panel.grid.major.x =  element_line(colour = "lightgray", size = 0.1),
        panel.border = element_blank()) +
  guides(fill = guide_legend(ncol=1)) +
  scale_x_continuous(expand = c(0, 0))
```

<img src="./figures/unnamed-chunk-10-1.png" width="672" />

```r


ggsave("figures/grouped_taxa.png", width=5, height=4)

```



```r
#---------------------------------------

phylum <- ps_df %>%
  filter(level == "phylum") 

phylum_pool <- phylum %>%
  group_by(taxon) %>%
  summarise(pool = max(rel_abund) < 0.02,
            mean = mean(rel_abund), .groups="drop")  
  
inner_join(phylum, phylum_pool, by="taxon") %>%
  mutate(taxon = if_else(pool, "Other", taxon)) %>%
  ggplot(aes(x = sample_id, y = 100*rel_abund, fill = taxon)) +
  theme_bw() + 
  geom_bar(stat = "identity", position = position_stack()) +
  labs(x="Sample", y="Relative Abundance", subtitle = "Taxa Relative Abundance faceted by group", fill = NULL) +
  facet_grid(~ nationality, space = "free", scales = "free") +
  theme(axis.text.x = element_blank(),
        legend.text = element_markdown(),
        strip.background = element_rect(colour = "lightblue", fill = "lightblue")) +
  scale_fill_manual(name=NULL,
                    breaks=c("*Actinobacteria*","*Bacteroidetes*", "*Firmicutes*","*Proteobacteria*", "*Verrucomicrobia*", "Other"),
                    labels=c("*Actinobacteria*","*Bacteroidetes*", "*Firmicutes*","*Proteobacteria*", "*Verrucomicrobia*", "Other"),
                    values = c(brewer.pal(5, "Paired"), "gray"))
```

<img src="./figures/unnamed-chunk-11-1.png" width="672" />

```r

 
ggsave("figures/faceted_phyla_bar.png", width=5, height=4)

```



```r
#---------------------------------------

genus <-  ps_df %>%
  filter(level == "genus") 

genus_pool <- genus %>%
  group_by(taxon) %>%
  summarise(pool = max(rel_abund) < 0.15,
            mean = mean(rel_abund), .groups="drop")  
  
inner_join(genus, genus_pool, by="taxon") %>%
  mutate(taxon = if_else(pool, "Other", taxon)) %>%
  ggplot(aes(x = sample_id, y = 100*rel_abund, fill = taxon)) +
  theme_bw() + 
  geom_bar(stat = "identity", position = position_stack()) +
  labs(x="Sample", y="Relative Abundance", subtitle = "Taxa Relative Abundance faceted by group", fill = NULL) +
  facet_grid(~ nationality, space = "free", scales = "free") +
  theme(axis.text.x = element_blank(),
        legend.text = element_markdown(),
        strip.background = element_rect(colour = "lightblue", fill = "lightblue"))
```

<img src="./figures/unnamed-chunk-12-1.png" width="672" />

```r


ggsave("figures/faceted_taxa_bar.png", width=5, height=4)
```



```r
#---------------------------------------

ps_rel %>% 
microbial::plotbar( level="Phylum", group = "nationality", top = 10) +
  theme(axis.text.x = element_text(angle = 0)) + 
  labs(x="Sample", 
       y="Relative Abundance", 
       subtitle = "Phyla Relative Abundance by plotbar()\nin microbial package", 
       fill = NULL)
```

<img src="./figures/unnamed-chunk-13-1.png" width="672" />

```r

ggsave("figures/plotbar_phyla_bar.png", width=5, height=4)

#---------------------------------------

ps_rel %>% 
microbial::plotbar(level="Genus", group = "nationality", top = 10) +
  theme(axis.text.x = element_text(angle = 0)) + 
  labs(x="Sample", 
       y="Relative Abundance", 
       subtitle = "Taxa Relative Abundance by plotbar()\nin microbial package", 
       fill = NULL)
```

<img src="./figures/unnamed-chunk-13-2.png" width="672" />

```r

ggsave("figures/plotbar_taxa_bar.png", width=5, height=4)

```



```r
#---------------------------------------
# Using ggpubr package
#---------------------------------------

library(ggpubr)
library(RColorBrewer)
library(viridis)

psmelt(ps_rel) %>% 
  select(Sample, nationality, Phylum, Abundance) %>% 
  group_by(Sample, nationality, Phylum) %>% 
  summarise(Abundance = sum(Abundance), .groups = "drop") %>% 
  filter(Abundance > 0.01) %>% 
  mutate(Phylum = fct_reorder(Phylum, -Abundance)) %>% 
  ggbarplot(x = "Phylum", y = "Abundance", 
            fill = "nationality",
            color = "nationality",
            palette = "jco",
            sort.by.groups = FALSE,     
            rotate = TRUE,
            position = position_dodge()) +
  theme_classic() +
  theme(axis.text.x = element_markdown(angle = 0, hjust = 1, vjust = 1),
        axis.text.y = element_markdown(),
        legend.text = element_markdown(),
        legend.key.size = unit(12, "pt"))  +
  labs(x=NULL, 
       y="Relative Abundance (%)",
       subtitle = "Phyla Relative Abundance by ggbarplot()\nin ggpubr package grouped by nationality", 
       fill = NULL) +
  guides(color = "none")
```

<img src="./figures/unnamed-chunk-14-1.png" width="672" />

```r

ggsave("figures/ggpubr_phyla_bar.png", width=5, height=4)
```



```r
#------------------------------
# Plot top n
#------------------------------
library(phyloseq)
library(tidyverse)
library(ggpubr)

df <- psmelt(ps_rel) %>% 
  select(Sample, nationality, Genus, Abundance) %>% 
  group_by(Sample, nationality, Genus) %>% 
  summarise(Abundance = sum(Abundance), .groups = "drop")

# Function
top_n_unique <- function(df, n, group_var, value_var) {
  df %>%
    group_by({{ group_var }}) %>%
    mutate(rank = dense_rank(desc({{ value_var }}))) %>%
    filter(rank <= n) %>%
    ungroup() %>%
    select(-rank)
}

# Plot
n <- 10  # Number of top unique taxa to select
top_n_unique(df, n, Genus, Abundance) %>% 
  filter(Abundance > 0.2) %>% 
  mutate(Genus = fct_reorder(Genus, -Abundance)) %>% 
  ggbarplot(x = "Genus", y = "Abundance", 
            fill = "nationality",
            color = "nationality",
            palette = "jco",
            sort.by.groups = FALSE,     
            rotate = TRUE,
            position = position_dodge()) +
  theme_classic() +
  theme(axis.text.x = element_markdown(angle = 0, hjust = 1, vjust = 1),
        axis.text.y = element_markdown(),
        legend.text = element_markdown(),
        legend.key.size = unit(12, "pt"))  +
  labs(x=NULL, 
       y="Relative Abundance (%)",
       subtitle = "Genera Relative Abundance by ggbarplot()\nin ggpubr package grouped by nationality", 
       fill = NULL) +
  guides(color = "none")
```

<img src="./figures/unnamed-chunk-15-1.png" width="672" />

```r

ggsave("figures/ggpubr_top_n_bar.png", width=5, height=4)
```




# (PART) DIFFERENTIAL ABUNDANCE {-}

# Differential Analysis of Microbiome Data 

Differential abundance analysis is a fundamental technique in microbiome research used to identify features that exhibit significant differences in abundance across different biological conditions or treatments. This analysis is crucial for understanding the dynamics of microbial communities and their response to various environmental factors. In this section, we will explore different methods for conducting differential abundance analysis using microbiome data.

## Using run_lefse() in microbiomeMaker Package

The `run_lefse()` function in the `microbiomeMaker` R package provides a convenient method for performing Differential Abundance Analysis on microbiome data. By leveraging the Linear Discriminant Analysis Effect Size (LEfSe) algorithm, `run_lefse()` enables researchers to interpret the results in the context of biological class labels, uncovering important insights into microbial community dynamics.

**Dataset: dietswap**

```r
library(phyloseq)
library(microbiomeMarker)

# Run LEfSe analysis
run_lefse(
  psextra_raw,
  wilcoxon_cutoff = 0.0001,
  group = "nationality",
  taxa_rank = "Genus",
  transform = "log10p",
  kw_cutoff = 0.01,
  multigrp_strat = TRUE,
  lda_cutoff = 2
) %>% 
plot_heatmap(group = "nationality", color = "rainbow")
```

<img src="./figures/lefse_microbiomeMarker1-1.png" width="960" />


**Dataset: caporaso**

```r
library(microbiomeMarker)

ps_caporaso <- ps_caporaso %>%  tax_fix()

run_lefse(
    ps_caporaso, 
    wilcoxon_cutoff = 0.001,
    group = "SampleType",
    taxa_rank = "Genus",
    transform = "log10p",
    kw_cutoff = 0.01,
    multigrp_strat = TRUE,
    lda_cutoff = 2) %>%
  plot_heatmap(group = "SampleType", color = "rainbow")
```

<img src="./figures/lefse_microbiomeMarker2-1.png" width="960" />



## Cladogram on kostic_crc dataset
Using plot_cladogram() from microbiomeMarker package

```r
library(microbiomeMarker)
kostic_crc_small <- phyloseq::subset_taxa(
    ps_kostic_crc,
    Phylum %in% c("Firmicutes")
)
mm_lefse <- run_lefse(
    kostic_crc_small,
    wilcoxon_cutoff = 0.01,
    group = "DIAGNOSIS",
    kw_cutoff = 0.01,
    multigrp_strat = TRUE,
    lda_cutoff = 4
)
plot_cladogram(mm_lefse, color = c("darkgreen", "red"))
```

<img src="./figures/ps_kostic_crc_cladogram-1.png" width="960" />



# (PART) ORDINATION {-}

# Ordination Plots

## Ordination with auto caption


```r
library(microViz)

# perform ordination
psextra_clr_dietswap %>%
  ord_calc(method = "PCA") %>%
  ord_plot(
    plot_taxa = 1:6, colour = "bmi_group", size = 1.5,
    tax_vec_length = 0.325,
    tax_lab_style = tax_lab_style(max_angle = 90, aspect_ratio = 0.5),
    auto_caption = 8
  )
```

<img src="./figures/unconstrained_aitchison_pca-1.png" width="672" />

## Customizing an ordination plot

```r
library(ggplot2)

# customise plot
psextra_clr_dietswap %>%
  ord_calc(method = "PCA") %>%
  ord_plot(
    plot_taxa = 1:6, colour = "bmi_group", size = 1.5,
    tax_vec_length = 0.325,
    tax_lab_style = tax_lab_style(max_angle = 90, aspect_ratio = 0.5),
    auto_caption = 8) +
  stat_ellipse(aes(linetype = bmi_group, colour = bmi_group), linewidth = 0.3) + # linewidth not size, since ggplot 3.4.0
  scale_colour_brewer(palette = "Set1") +
  theme(legend.position = "right") +
  coord_fixed(ratio = 0.5, clip = "off") # makes rotated labels align correctly
```

<img src="./figures/unnamed-chunk-16-1.png" width="672" />


# (PART) RELATIONSHIP ANALYSIS {-}

# Patterns and Associations Analysis

In this section, we will explore methods for analyzing patterns and associations within datasets. These methods are instrumental in uncovering relationships between variables, identifying underlying patterns, and elucidating complex structures within the data. Specifically, we will discuss correlation analysis, heatmap visualization, hierarchical clustering, and other techniques for exploring patterns and associations in various types of datasets.

## Correlation analysis and visualization

Correlation analysis quantifies the strength and direction of the linear relationship between two continuous variables. It helps in understanding how variables are related to each other and can identify potential associations or dependencies within the data. Heatmap visualization is commonly integrated with correlation analysis to visually represent correlation matrices, allowing for the identification of patterns and relationships between variables.


### Pearson Correlation Heatmap with microViz Package
In this section, we utilize the cor_heatmap function available in the microViz package to generate a heatmap visualizing Pearson correlation coefficients between variables. This technique allows for the exploration of relationships and patterns within the data, providing insights into the strength and direction of linear associations between variables.



```r
ps_extra <- ps_dietswap %>%
tax_transform(
  trans = "clr",
  rank = "Genus",
  keep_counts = TRUE,
  zero_replace = 0,
  add = 0,
  transformation = NULL
)
```



```r
# Load required packages
library(phyloseq)
library(microViz)

# Set up the data with numerical variables and filter to top taxa
ps <- ps_dietswap %>%
  ps_mutate(
    bmi = recode(bmi_group, obese = 3, overweight = 2, lean = 1),
    female = if_else(sex == "female", true = 1, false = 0),
    african = if_else(nationality == "AFR", true = 1, false = 0)
  ) %>%
  tax_filter(
    tax_level = "Genus", min_prevalence = 1 / 10, min_sample_abundance = 1 / 10
  ) %>%
  tax_transform("identity", rank = "Genus")

# Randomly select 30 taxa from the 50 most abundant taxa (just for demo)
set.seed(123)
top_taxa <- tax_top(ps_dietswap, n = 10, by = sum, rank = "unique", use_counts = FALSE)

taxa <- sample(tax_top(ps_dietswap, n = 50), size = 30)

# Clean the data and draw the correlation heatmap 
cor_heatmap(
  data = ps, 
  taxa = taxa,
  taxon_renamer = function(x) stringr::str_remove(x, " [ae]t rel."),
  tax_anno = taxAnnotation(
    Prev. = anno_tax_prev(undetected = 50),
    Log2 = anno_tax_box(undetected = 50, trans = "log2", zero_replace = 1)
  )
)
```

<img src="./figures/unnamed-chunk-18-1.png" width="672" />

> Note: The parameter min_prevalence = 0.1 is set, equivalent to approximately 23 out of 222 samples. This value may vary depending on the dataset. Adjust as necessary.


## Hierarchical clustering

Hierarchical clustering is a method for grouping similar observations or variables together based on their similarity or dissimilarity. It creates a hierarchical tree-like structure (dendrogram) that illustrates the relationships between clusters, allowing for the identification of natural groupings or patterns within the data.


# (PART) PHYLOGENY ANALYSIS {-}

# Phylogenetic Analysis 

Phylogenetic analysis is essential in microbiome research for understanding the evolutionary relationships between microbial taxa. In this section, we will explore various tools and techniques for phylogenetic analysis. 

## View a tree by ape and ggtree R packages

```r
set.seed(1234)

library(ape)
library(ggtree)

tree <- rtree(20)
ggtree(tree)
```

<img src="./figures/unnamed-chunk-19-1.png" width="672" />

## Cladogram using ggtree() function


```r
library(ape)
library(ggtree)


ggtree(
  phy_tree(ps_caporaso),
  mapping = NULL,
  layout = "circular",
  open.angle = 0,
  mrsd = NULL,
  as.Date = FALSE,
  yscale = "none",
  yscale_mapping = NULL,
  ladderize = TRUE,
  right = FALSE,
  branch.length = "branch.length",
  root.position = 0,
  xlim = NULL
)
```

<img src="./figures/unnamed-chunk-20-1.png" width="576" />


