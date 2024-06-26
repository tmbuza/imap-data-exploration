# scripts/create_barplots.R
library(readr)
library(dplyr)
library(ggplot2)
library(svglite)

mo <- read_csv("data/mothur_composite.csv", show_col_types = FALSE) %>% 
  mutate(pipeline="MOTHUR", .before=2) %>% 
  select(-otu, -3)

q2 <- read_csv("data/qiime2_composite.csv", show_col_types = FALSE) %>% 
  mutate(pipeline="QIIME2", .before=2) %>% 
  select(-feature, -c(3:14))

rbind(mo, q2) %>% 
  group_by(sample_id) %>%
  mutate(total = sum(count)) %>%
  filter(total > 0) %>%
  group_by(Genus) %>%
  mutate(total = sum(count)) %>%
  # filter(total != 0) %>%
  filter(total >= 10) %>%
  ungroup() %>%
  select(-total) %>% 
  ggplot(aes(x=Genus, y=rel_abund, fill=pipeline)) +
  facet_grid(~pipeline) +
  geom_col() +
  coord_flip() +
  labs(fill="Pipeline")

ggsave(file="figures/taxon_barplot.png", width=10, height=10)
ggsave(file="figures/taxon_barplot.svg", width=10, height=10)