library(decontam); packageVersion("decontam")
library(tidyverse); packageVersion("tidyverse")
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(vegan); packageVersion("vegan")
library(ape); packageVersion("ape")
library(cowplot); packageVersion("cowplot")
library(EcolUtils); packageVersion("EcolUtils")

permutations = 9999
rarefaction_depth = 6067
set.seed(1851)

table <- read.delim("data/table.txt") %>%
  column_to_rownames(var = "Feature.ID")
OTU = otu_table(table, taxa_are_rows = TRUE)

metadata <- read.csv("data/metadata.csv", row.names = 1) %>% 
  as.data.frame() %>% 
  mutate(sample_type = case_when(Alt_SampleID == "NEG" ~ "NEG",
                                 TRUE ~ "POS"))

metadata <- metadata %>% mutate(sampleid = rownames(metadata))

metdata_ps = sample_data(metadata, errorIfNULL = T)

taxonomy.horiz <- read.delim("data/taxonomy.tsv") %>% 
  select(Feature.ID, Taxon) %>% 
  rename(featureid = "Feature.ID")

taxonomy = read.delim("data/taxonomy.tsv") %>% 
  select(Feature.ID, Taxon) %>% 
  rename(featureid = "Feature.ID") %>% 
  separate(col = "Taxon", 
           into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"),
           sep = "; ",
           fill = "right",
           extra = "drop") %>% 
  column_to_rownames(var = "featureid") %>% 
  as.matrix()
  
tax = tax_table(taxonomy)

ps <- phyloseq(OTU, tax, metdata_ps)

sample_data(ps)$is.neg <- sample_data(ps)$sample_type == "NEG"
contamdf.prev <- isContaminant(ps, method="prevalence", neg="is.neg")
table(contamdf.prev$contaminant)
  
true <- contamdf.prev %>% 
  filter(contaminant == "FALSE") %>% 
  row.names() %>% 
  as.array() 

filt.table <- table %>% filter(rownames(table) %in% true) %>% as.data.frame()

samples_pass_rarefaction <- filt.table %>% 
  rownames_to_column(var = "featureid") %>% 
  pivot_longer(-featureid,
               names_to = "sampleid",
               values_to = "count") %>% 
  group_by(sampleid) %>% 
  summarise(feature.count = sum(count)) %>% 
  arrange((feature.count)) %>% 
  filter(feature.count > rarefaction_depth) %>% 
  arrange(sampleid) %>% 
  pull(sampleid) %>% 
  as.array()

rarefied.table <- filt.table %>% 
  t() %>% 
  rrarefy(rarefaction_depth) %>% 
  t() %>% 
  as.data.frame()

table.to.save <- rarefied.table %>% 
  rownames_to_column(var = "featureid") %>% 
  select(featureid, samples_pass_rarefaction) %>% 
  left_join(., taxonomy.horiz)

# write_csv(table.to.save, file = "data/filtered_rarefied_table.csv")

###############################
### PCoA

## Filter feature table to only samples in metadata.filtered
## duplicate rownames and featureid in col 1, but select will solve that
metadata.filt <- metadata %>% filter(Brand == "AnVision" |
                                       Brand == "CaraLens" |
                                       Brand == "Dioptrix" |
                                       Brand == "Negative Control")

filt.table <- rarefied.table %>% 
  select(metadata.filt$sampleid)

table.t <- filt.table %>% 
  as.data.frame() %>%
  # select(metadata$sampleid) %>% 
  t()

#table.t <- table %>% t()
## PCoA plots for all sample types
#quarter root transformation of table
table.transformed <- table.t^1/4
  
plot.pcoa <- function(distance){
  dist <- vegdist(table.transformed, method= distance)
  permanova <- adonis(dist ~ Brand, 
                      permutations = permutations, 
                      data = metadata.filt)
  
  pcoa <- pcoa(dist, correction = "cailliez")
  
  percent_var = pcoa$values$Eigenvalues/pcoa$trace * 100
  
  pcoa_vectors <- pcoa$vectors %>% as_tibble(rownames = "sampleid") %>% 
    select(sampleid, Axis.1,Axis.2)
  
  colnames(pcoa_vectors) <- c("sampleid", "PCo1", "PCo2")
  
  hull <- as.data.frame(pcoa_vectors) %>% 
    left_join(metadata.filt) %>% 
    group_by(Brand) %>% 
    slice(chull(PCo1, PCo2))
  
  variance_rep <- round(percent_var[1:2],2)
  
  pcoa.vectors.metadata <- inner_join(metadata.filt, pcoa_vectors, by = "sampleid") 
  
  ggplot(pcoa.vectors.metadata, aes(x = PCo1, 
                                    y = PCo2, 
                                    color = Brand)) +
    
    geom_point(aes(shape = factor(Sample.1, levels = c("IOL", "IF",
                                                       "E-swab", "lysis buffer"))),
               size = 3) +
    stat_ellipse(show.legend = FALSE) +
    
    scale_color_brewer(palette = "Dark2") +
    
    scale_shape_manual(values = c(15,16,17,18),
                       labels = c("IOL", "IF",
                                  "E-Swab",
                                  "Lysis Buffer")) +
    
    labs(x = paste0("PCo1 - ", variance_rep[1], "%"),
         y = paste0("PCo1 - ", variance_rep[2], "%"),
         caption = paste0("F: ", round(permanova$aov.tab$F.Model,2),
                           "; p = ", round(permanova$aov.tab$`Pr(>F)`,7)))+
    
    theme_cowplot() +
    
    theme(
      axis.title = element_text(face = "bold"),
      axis.text = element_text(face = "bold"),
      
      legend.title = element_blank(),
      legend.text = element_text(face = "bold"),
      
      plot.caption = element_text(face = "bold",
                                  size = 14)
    ) +
    
    guides(color = guide_legend(override.aes = list(size = 3)))
  
}

bc.plot <- plot.pcoa("bray")
j.plot <- plot.pcoa("jaccard")

plot_grid(bc.plot + theme(plot.subtitle = element_blank(),),
          j.plot + theme(plot.subtitle = element_blank()),
          nrow = 2,
          ncol = 1,
          labels = "AUTO",
          scale = 0.95)

# ggsave("plots/pcoa.tiff",
#        height = 9,
#        width = 8,
#        units = c("in"),
#        dpi = 600,
#        bg = "white")

dist.bc <- vegdist(table.transformed, method= "bray")
dist.j <- vegdist(table.transformed, method= "jaccard")

permanova.bc <- adonis(dist.bc ~ Brand * Sample.1, permutations = permutations, data = metadata.filt)
permanova.j <- adonis(dist.j ~ Brand * Sample.1, permutations = permutations, data = metadata.filt)


pairwise.bc <- adonis.pair(dist.bc, as.factor(metadata.filt$Brand), nper = permutations)
pairwise.j <- adonis.pair(dist.j, as.factor(metadata.filt$Brand), nper = permutations)


# write.csv(pairwise.bc,
#           file = "data/bc.pairwise.csv")
# write.csv(pairwise.j,
#           file = "data/j.pairwise.csv")

