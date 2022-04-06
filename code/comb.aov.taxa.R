library(tidyverse)
library(readxl)
aov_res <- read.csv("data/metabo_anova_res.csv", row.names = 1) %>% 
  rownames_to_column(var = "ASV")
taxonomy <- read_excel("data/taxonomy.filtered.xlsx") %>% 
  select(ASV, Taxon)


comb.table <- left_join(aov_res, taxonomy, by = "ASV")

write.csv(comb.table, 
          file = "data/aov_res_taxa.csv")
