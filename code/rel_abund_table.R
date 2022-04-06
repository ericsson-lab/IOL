library(tidyverse)
library(readxl)

table <- read.csv("data/filtered_rarefied_table.csv")
metadata <- read.csv("data/metadata.csv") %>% 
  filter(Brand == "AnVision" |
           Brand == "CaraLens" |
           Brand == "Dioptrix") %>% 
  rename(sampleid = "Sample")
taxonomy <- read_excel("data/taxonomy.filtered.xlsx")



abund.list <- table %>% 
  select(-Taxon, featureid, contains("KD")) %>% 
  pivot_longer(-featureid, 
               names_to = "sampleid",
               values_to = "count") %>% 
  group_by(sampleid) %>% 
  mutate(rel_abund = count/sum(count) *100) %>%
  inner_join(., metadata) %>% 
  group_by(Brand, Sample.1, featureid) %>% 
  summarise(mean_rel_abund = mean(rel_abund),
            sd_rel_abund = sd(rel_abund)) %>% 
  filter(mean_rel_abund >= 5)

abund.list %>% 
  inner_join(., taxonomy) %>% 
  select(ASV, Taxon, Brand, Sample.1, mean_rel_abund, sd_rel_abund) %>% 
  write_csv(file = "data/rel_abund_table.csv")

