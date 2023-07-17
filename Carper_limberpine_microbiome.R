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
nr_rmnp <-subset_samples(Rangewide,SITE%in%c("NR_FOR")|SITE%in%c("NR_TL")|SITE%in%c("RMNP"))


################################################

#Alpha diversity with plot and significance test (example script for pifl used for all alpha diversity analyses)

pifl_rare <- phyloseq::rarefy_even_depth(pifl, rngseed = 123, replace = FALSE)
Observed = phyloseq::estimate_richness(pifl_rare, measures = "Observed")
Shannon = phyloseq::estimate_richness(pifl_rare, measures = "Shannon")
Site = phyloseq::sample_data(pifl_rare)$SITE
adiv <- data.frame(Observed,Shannon,Site)
gather(key = metric, value = value, c(Observed, Shannon,)) %>%
mutate(metric = factor(metric, levels = c("Observed", "Shannon"))) %>%
ggplot(aes(x = Site, y = value)) +
 scale_color_manual(values=SiteColors) +
geom_boxplot(outlier.color = NA) +
geom_jitter(aes(color = Site), height = 0, width = .2) +
labs(x = "", y = "") +
facet_wrap(~ metric, scales = "free") +
theme(legend.position="none") +
theme(axis.text.x = element_text(angle = 90))
kruskal.test(Observed ~ Site, data = adiv)

################################################

#Beta diversity analysis, example is from set 1, table 2, but script was used for all analyses (sites, species, balanced sets, age)

#CSS normalization with metagenomeSeq
nr_rmnp_m <- phyloseq_to_metagenomeSeq(nr_rmnp)
nr_rmnp_mcss = cumNorm( nr_rmnp_m , p=cumNormStatFast(nr_rmnp_m) )
nr_rmnp_css <-nr_rmnp
otu_table(nr_rmnp_css) <- otu_table(MRcounts(nr_rmnp_mcss, norm = TRUE),taxa_are_rows = TRUE)

#Get distance matrices (Weighted unifrac and bray curtis)
nr_rmnp_css_uni <- phyloseq::distance(nr_rmnp_css, method = "wunifrac")
nr_rmnp_css_bray <- phyloseq::distance(nr_rmnp_css, method = "bray")

#PERMANOVA to test for dissimilarity among groups (sites and species example here)
vegan::adonis(nr_rmnp_css_bray ~ phyloseq::sample_data(nr_rmnp_css)$SITE*phyloseq::sample_data(nr_rmnp_css)$HostSpecies)
vegan::adonis(nr_rmnp_css_uni ~ phyloseq::sample_data(nr_rmnp_css)$SITE*phyloseq::sample_data(nr_rmnp_css)$HostSpecies)

#betadisper to test for homogeneity of dispersion among groups
beta <- betadisper(pifl_css_bray, sample_data(pifl_css)$SITE)
permutest(beta)

#PCoA: ordination and plots
nr_rmnp_css_ordu_bray <- ordinate(nr_rmnp_css,"PCoA","bray")
plot_ordination(nr_rmnp_css,nr_rmnp_css_ordu_bray,color="SITE") + scale_color_manual(values=SiteColors)+ geom_point(size = 2, alpha=0.5) + stat_ellipse(aes(group = SITE), linetype = 2)
plot_ordination(nr_rmnp_css,nr_rmnp_css_ordu_bray,color="HostSpecies") + scale_color_manual(values=SpeciesColors)+ geom_point(size = 2, alpha=0.5) + stat_ellipse(aes(group = HostSpecies), linetype = 2)

################################################

#indicator species analysis (example is for indicator species at the RMBL site)

#aggregate at family level
rmbl.fam.inv <- microbiome::aggregate_rare(rmbl, "Family",detection = 2/100, prevalence =30/100 )

#indicator species analysis at family level
OTU1 = as(otu_table(pifl_rmbl.fam.inv), "matrix")
if(taxa_are_rows(pifl_rmbl.fam.inv)){OTU1 <- t(OTU1)}
OTUdf = as.data.frame(OTU1)
pifl_rmbl_fam_inv = multipatt(OTUdf, phyloseq::sample_data(pifl_rmbl.fam.inv)$Age,func = "r.g", control = how(nperm=999))
print("pifl_rmbl_fam_inv")
summary(pifl_rmbl_fam_inv)

#variation partitioning analysis

#reading files
climate  <- read.delim("Mapping_climate_best_param.txt", row.names = 1)
landcover10 <-read.delim("Mapping_varpart_LC_cat_10km.txt", row.names = 1)
braycurtis <- read.delim("bc_distance.tsv" , row.names = 1)

varp10 <-varpart(braycurtis,climate,landcover10)
plot(varp10)
