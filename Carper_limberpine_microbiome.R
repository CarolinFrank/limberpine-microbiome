library(tidyverse)
library(phyloseq)
library(DESeq2)
library(microbiome)
library(vegan)
library(picante)
library(glmnet)
library(ALDEx2)
library(HMP)
library(dendextend)
library(rms)
library(qiime2R)
library(colormap)
library(cowplot)
library(indicspecies)
library(metagenomeseq)

#Importing Qiime 2 artifacts to phyloseq objects
Rangewide<-qza_to_phyloseq(features="otu_table_500.qza",tree="tree.qza",taxonomy="taxonomy.qza",metadata="Mapping_500_q2R.txt")
Rangewide<-qza_to_phyloseq(features="otu_table_500.qza",tree="tree.qza",taxonomy="taxonomy_table.qza",metadata="Mapping_500_q2R.txt

#Subdividing phyloseq object into subsamples; species, sites, species/site combinations, etc (a few examples shown here)

all_surface <- subset_samples(Rangewide,SampleType%in%c("PHYLLOSPHERE"))
pifl <-subset_samples(Rangewide,HostSpecies%in%c("PIFL"))
rmbl <-subset_samples(Rangewide,SITE%in%c("RMBL"))

################################################

#Alpha diversity (example script used for all alpha diversity analyses)

pifl_rare <- phyloseq::rarefy_even_depth(pifl, rngseed = 123, replace = FALSE)
Observed = phyloseq::estimate_richness(pifl_rare, measures = "Observed")
Shannon = phyloseq::estimate_richness(pifl_rare, measures = "Shannon")
Site = phyloseq::sample_data(pifl_rare)$SITE
adiv <- data.frame(Observed,Shannon,Site)
kruskal.test(Observed ~ Site, data = adiv)
