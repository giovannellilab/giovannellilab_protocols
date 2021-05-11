library(dada2); packageVersion("dada2")

path <- "cordone_fastq" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

ls()

plotQualityProfile(fnFs[1:2])

plotQualityProfile(fnRs[1:2])

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(155,145),
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              trimLeft=17, trimRight=15,
              compress=TRUE, multithread=TRUE)
head(out)

errF <- learnErrors(filtFs, nbases = 5e8, multithread=TRUE, randomize=TRUE)

errR <- learnErrors(filtRs, nbases = 1e9, multithread=TRUE, randomize=TRUE)

plotErrors(errF, nominalQ=TRUE)

derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

dadaFs <- dada(derepFs, err=errF, pool="pseudo", multithread=TRUE)

dadaRs <- dada(derepRs, err=errR, pool="pseudo", multithread=TRUE)

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

seqtab <- makeSequenceTable(mergers)

dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

#seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(153,161)]

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

dim(seqtab)

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

taxa <- assignTaxonomy(seqtab.nochim, "~/silva_nr_v132_train_set.fa.gz", multithread=TRUE)

taxa <- addSpecies(taxa, "~/silva_species_assignment_v132.fa.gz")

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
write.csv(taxa.print, "taxa_print.csv")

### Load required libraries
library(microbiome) # data analysis and visualisation
library(phyloseq) # also the basis of data object. Data analysis and visualisation
library(RColorBrewer) # nice color options
library(dplyr) # data handling
library(network) # networks
library(intergraph)  # networks
library(ggnet)   # network plotting with ggplot
library(igraph)  # networks
library(ggplot2) # plotting library
library(gridExtra) # gridding plots
library(ape) # importing and handling phylogenetic trees
library(ggthemes) # additional themes fro ggplot2
library(magrittr) #
library(rioja) # plotting poackages for tabular bubbleplots
library(see) #viasualization and half-violin plots
library(ggpubr)
#library(ggtern)
library(plyr)
library(coda.base)
library(vegan)
library(propr)
library(msa)

theme_set(theme_bw())

ls()

rownames(seqtab.nochim)

# Make the phylogenetic tree
seqs <- getSequences(seqtab.nochim)

names(seqs) <- seqs # This propagates to the tip labels of the tree
mult <- msa(seqs, method="ClustalW", type="dna", order="input")

library(phangorn)

phang.align <- as.phyDat(mult, type="DNA", names=getSequence(seqtab))
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)

## negative edges length changed to 0!
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                       rearrangement = "stochastic", control = pml.control(trace = 0))

prok_sample <- read.csv("dataset_env.csv", header=T, sep=",", row.names=1)

prok_sample

prok_data <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = F), 
                     phy_tree(treeNJ), 
                      tax_table(taxa), 
                      sample_data(prok_sample))
prok_data

readcount(prok_data)

write_phyloseq(prok_data, type = "all")

# Removing the stationsfrom the deep stations 
prok_data <- subset_samples(prok_data, sample_names(prok_data) != "A2-Cordone" & sample_names(prok_data) !="D2-Cordone" &  sample_names(prok_data) != "8B-Cordone" & sample_names(prok_data) !="E2-Cordone")
prok_data = filter_taxa(prok_data, function(x) sum(x) > 0, TRUE) # After removing samples filter the taxa left with zero global abundance
prok_data # get stats on the dataset

readcount(prok_data)

# Clean up unwanted sequences from mitochrondria e chloroplast
#prok_data <- subset_taxa(prok_data, Family != "Mitochondria" & Order != "Chloroplast")
prok_data <- subset_taxa(prok_data, (Order!="Chloroplast") | is.na(Order))
prok_data <- subset_taxa(prok_data, (Family!="Mitochondria") | is.na(Family))
prok_data
readcount(prok_data)

# Removing the potential human pathogens and contaminants
prok_data2 <- subset_taxa(prok_data,  (Genus != "Lactococcus") | is.na(Genus))  
prok_data2 <- subset_taxa(prok_data2,  (Genus != "Lactobacillus") | is.na(Genus))  
prok_data2 <- subset_taxa(prok_data2,  (Genus != "Cutibacterium") | is.na(Genus))  
prok_data2 <- subset_taxa(prok_data2,  (Genus != "Enterococcus") | is.na(Genus)) 
prok_data2 <- subset_taxa(prok_data2,  (Genus != "Streptococcus") | is.na(Genus)) 
prok_data2 <- subset_taxa(prok_data2,  (Genus != "Acinetobacter") | is.na(Genus)) 
prok_data2 <- subset_taxa(prok_data2,  (Genus != "Citrobacter") | is.na(Genus)) 
prok_data2 <- subset_taxa(prok_data2,  (Genus != "Bifidobacterium") | is.na(Genus)) 
prok_data2 <- subset_taxa(prok_data2,  (Genus != "Proprionibacterium") | is.na(Genus)) 
prok_data2 <- subset_taxa(prok_data2,  (Order != "Enterobacteriales") | is.na(Order)) 
prok_data2 <- subset_taxa(prok_data2,  (Genus != "Corynebacterium_1") | is.na(Genus)) 
prok_data2 <- subset_taxa(prok_data2,  (Genus != "Corynebacterium") | is.na(Genus)) 
prok_data2 <- subset_taxa(prok_data2,  (Genus != "Escherichia/Shigella") | is.na(Genus))
prok_data2 <- subset_taxa(prok_data2, (Genus != "Abiotrophia") |  is.na(Genus))
prok_data2 <- subset_taxa(prok_data2, (Genus != "Achromobacter") | is.na(Genus))
prok_data2 <- subset_taxa(prok_data2, (Genus != "Actinobacillus") | is.na(Genus))
prok_data2 <- subset_taxa(prok_data2, (Genus != "Arcanobacterium") | is.na(Genus))
prok_data2 <- subset_taxa(prok_data2, (Genus != "Babesia") | is.na(Genus))
prok_data2 <- subset_taxa(prok_data2, (Genus != "Bacillus") | is.na(Genus))
prok_data2 <- subset_taxa(prok_data2, (Genus != "Bartonella") | is.na(Genus))
prok_data2 <- subset_taxa(prok_data2, (Genus != "Bordetella") |  is.na(Genus))
prok_data2 <- subset_taxa(prok_data2, (Genus != "Borrelia") | is.na(Genus))
prok_data2 <- subset_taxa(prok_data2, (Genus != "Brodetella") |  is.na(Genus))
prok_data2 <- subset_taxa(prok_data2, (Genus != "Brucella") | is.na(Genus))
prok_data2 <- subset_taxa(prok_data2, (Genus != "Capnocytophaga") | is.na(Genus))
prok_data2 <- subset_taxa(prok_data2, (Genus != "Chlamydia") | is.na(Genus))
prok_data2 <- subset_taxa(prok_data2, (Genus != "Comamonas") | is.na(Genus))
prok_data2 <- subset_taxa(prok_data2, (Genus != "Coxiella") | is.na(Genus))
prok_data2 <- subset_taxa(prok_data2, (Genus != "Cronobacter") | is.na(Genus))
prok_data2 <- subset_taxa(prok_data2, (Genus != "Deinococcus") | is.na(Genus))
prok_data2 <- subset_taxa(prok_data2, (Genus != "Dermatophilus") | is.na(Genus))
prok_data2 <- subset_taxa(prok_data2, (Genus != "Ehrlichia") | is.na(Genus))
prok_data2 <- subset_taxa(prok_data2, (Genus != "Enterococcus") | is.na(Genus))
prok_data2 <- subset_taxa(prok_data2, (Genus != "Erysipelothrix") |  is.na(Genus))
prok_data2 <- subset_taxa(prok_data2, (Genus != "Escherichia") | is.na(Genus))
prok_data2 <- subset_taxa(prok_data2, (Genus != "Francisella") | is.na(Genus))
prok_data2 <- subset_taxa(prok_data2, (Genus != "Gardnerella") | is.na(Genus))
prok_data2 <- subset_taxa(prok_data2, (Genus != "Granulicatella") | is.na(Genus))
prok_data2 <- subset_taxa(prok_data2, (Genus != "Haemophilus") | is.na(Genus))
prok_data2 <- subset_taxa(prok_data2, (Genus != "Hafnia") | is.na(Genus))
prok_data2 <- subset_taxa(prok_data2, (Genus != "Helicobacter") | is.na(Genus))
prok_data2 <- subset_taxa(prok_data2, (Genus != "Klebsiella") | is.na(Genus))
prok_data2 <- subset_taxa(prok_data2, (Genus != "Kocuria") | is.na(Genus))
prok_data2 <- subset_taxa(prok_data2, (Genus != "Lawsonia") | is.na(Genus))
prok_data2 <- subset_taxa(prok_data2, (Genus != "Legionella") | is.na(Genus))
prok_data2 <- subset_taxa(prok_data2, (Genus != "Leptospira") | is.na(Genus))
prok_data2 <- subset_taxa(prok_data2, (Genus != "Listeria") | is.na(Genus))
prok_data2 <- subset_taxa(prok_data2, (Genus != "Merkel_cell") | is.na(Genus))
prok_data2 <- subset_taxa(prok_data2, (Genus != "Micrococcus") | is.na(Genus))
prok_data2 <- subset_taxa(prok_data2, (Genus != "Morganella") | is.na(Genus))
prok_data2 <- subset_taxa(prok_data2, (Genus != "Mycoplasma") | is.na(Genus))
prok_data2 <- subset_taxa(prok_data2, (Genus != "Neisseria") | is.na(Genus))
prok_data2 <- subset_taxa(prok_data2, (Genus != "Nocardia") | is.na(Genus))
prok_data2 <- subset_taxa(prok_data2, (Genus != "Pasteurella") | is.na(Genus))
prok_data2 <- subset_taxa(prok_data2, (Genus != "Plesiomonas") | is.na(Genus))
prok_data2 <- subset_taxa(prok_data2, (Genus != "Propionibacterium") | is.na(Genus))
prok_data2 <- subset_taxa(prok_data2, (Genus != "Proteus") | is.na(Genus))
prok_data2 <- subset_taxa(prok_data2, (Genus != "Providencia") | is.na(Genus))
prok_data2 <- subset_taxa(prok_data2, (Genus != "Rothia") | is.na(Genus))
prok_data2 <- subset_taxa(prok_data2, (Genus != "Salmonella") | is.na(Genus))
prok_data2 <- subset_taxa(prok_data2, (Genus != "Serratia") | is.na(Genus))
prok_data2 <- subset_taxa(prok_data2, (Genus != "Shewanella") | is.na(Genus))
prok_data2 <- subset_taxa(prok_data2, (Genus != "Shigella") | is.na(Genus))
prok_data2 <- subset_taxa(prok_data2, (Genus != "Sphaerophorus") | is.na(Genus))
prok_data2 <- subset_taxa(prok_data2, (Genus != "Staphylococcus") | is.na(Genus))
prok_data2 <- subset_taxa(prok_data2, (Genus != "Stenotrophomonas") | is.na(Genus))
prok_data2 <- subset_taxa(prok_data2, (Genus != "Streptococcus") | is.na(Genus))
prok_data2 <- subset_taxa(prok_data2, (Genus != "Treponema") | is.na(Genus))
prok_data2 <- subset_taxa(prok_data2, (Genus != "Yersinia") | is.na(Genus))


prok_data2
readcount(prok_data2)

prok_data2 <- subset_taxa(prok_data2,  (Kingdom != "Eukaryota") | is.na(Kingdom)) 
prok_data2
readcount(prok_data2)

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

prok_data2 <- subset_taxa(prok_data2,  (Kingdom != "Eukaryota") | is.na(Kingdom))

prok_data2
sum(readcount(prok_data2))

# Filter ASV with less than 1 reads across all the samples
prok_data2 <- prune_taxa(taxa_sums(prok_data2) > 1, prok_data2)
prok_data2
sum(readcount(prok_data2))

# After these filtering step remove the ASV and samples that are left with all zeros
prok_data2 = filter_taxa(prok_data2, function(x) sum(x) > 0, TRUE)
prok_data2
sum(readcount(prok_data2))

sum(readcount(prok_data2))

# Normalize the counts across the different samples by converting the abundance to
# relative abundance and multiply by the median library size
bac_ra <- transform_sample_counts(prok_data2, function(x) ((x / sum(x))*median(readcount(prok_data2))))

median(readcount(prok_data2))

readcount(bac_ra)

# Transform Bacteria abundance to relative abundances for plotting and some stats
bac_ra = transform_sample_counts(bac_ra, function(x){x / sum(x)})

## Agglomerate at a specific taxonomic level at the Genus level
bac_ra_genus = tax_glom(bac_ra, "Genus", NArm = TRUE)
bac_ra_family = tax_glom(bac_ra, "Family", NArm = TRUE)
bac_ra_order = tax_glom(bac_ra, "Order", NArm = TRUE)
bac_ra_class = tax_glom(bac_ra, "Class", NArm = TRUE)
bac_ra_phyla = tax_glom(bac_ra, "Phylum", NArm = TRUE)

# Plot Relative Abundance by Class for each station
plot_bar(bac_ra_phyla, fill="Kingdom", x="station", title = "Diversity at Kingdom level") +
theme_hc()  

# Plot Relative Abundance by Class for each station
plot_bar(bac_ra_phyla, fill="Phylum", x="station", title = "Diversity at Phylum level") +
theme_hc()
#as_ggplot(get_legend(plot_bar(bac_ra_class, fill="Phylum", x="station", title = "Diversity at Phylum level")))
ggsave(file="phylum_diversity.svg", width=10, height=8)

# Plot Relative Abundance by Class for each station
plot_bar(bac_ra_class, fill="Class", x="station", title = "Diversity at Class level") +
theme_hc() # + theme(legend.position = "none")
#as_ggplot(get_legend(plot_bar(bac_ra_class, fill="Class", x="station", title = "Diversity at Class level")))
ggsave(file="class_diversity.svg", width=10, height=8)

# Plot Relative Abundance by Order for each station
plot_bar(bac_ra_order, fill="Order", x="station", title = "Diversity at Order level") +
theme_hc() + theme(legend.position = "none")
ggsave(file="order_diversity.svg", width=10, height=8)
as_ggplot(get_legend(plot_bar(bac_ra_class, fill="Class", x="station", title = "Diversity at Order level")))
ggsave(file="order_diversity_legend.svg", width=10, height=8)

# Plot Relative Abundance by Order for each station
plot_bar(bac_ra_family, fill="Family", x="station", title = "Diversity at Family level") +
theme_hc() + theme(legend.position = "none")
ggsave(file="family_diversity.svg", width=10, height=8)
as_ggplot(get_legend(plot_bar(bac_ra_family, fill="Family", x="station", title = "Diversity at Family level")))
ggsave(file="family_diversity_legend.svg", width=14, height=8)

# Top Families
names(sort(taxa_sums(bac_ra_family), TRUE)[1:10])

get_taxa_unique(bac_ra_family, "Family")

plot_bar(subset_taxa(bac_ra_order, Phylum == "Proteobacteria"), fill="Order", x="station", title = "Order level diversity within the Proteobacteria") +
    theme_hc() #   + theme(legend.position = "none")
#as_ggplot(get_legend(plot_bar(subset_taxa(bac_ra_order, Phylum == "Proteobacteria"), fill="Order", x="station", title = "Diversity at Class level")))
ggsave(file="proteobacteria_diversity.svg", width=10, height=8)

plot_bar(subset_taxa(bac_ra_family, Class == "Gammaproteobacteria"), fill="Family", x="station", title = "Family level diversity within the Gammaproteobacteria") +
    theme_hc()# + theme(legend.position = "none")
#as_ggplot(get_legend(plot_bar(subset_taxa(bac_ra_family, Class == "Gammaproteobacteria"), fill="Family", x="station", title = "Diversity at Class level")))
ggsave(file="Gammaproteobacteria_diversity.svg", width=10, height=8)

plot_bar(subset_taxa(bac_ra_genus, Class == "Gammaproteobacteria"), fill="Genus", x="station", title = "Family level diversity within the Gammaproteobacteria") +
    theme_hc()# + theme(legend.position = "none")
#as_ggplot(get_legend(plot_bar(subset_taxa(bac_ra_genus, Class == "Gammaproteobacteria"), fill="Genus", x="station", title = "Diversity at Class level")))
ggsave(file="Gammaproteobacteria_genus_diversity.svg", width=10, height=8)

plot_bar(subset_taxa(bac_ra_family, Class == "Alphaproteobacteria"), fill="Family", x="station", title = "Family level diversity within the Alphaproteobacteria") +
    theme_hc() # + theme(legend.position = "none")
#as_ggplot(get_legend(plot_bar(subset_taxa(bac_ra_family, Class == "Alphaproteobacteria"), fill="Family", x="station", title = "Diversity at Class level")))
ggsave(file="alphaproteobacteria_diversity.svg", width=10, height=8)

plot_bar(subset_taxa(bac_ra_family, Phylum == "Bacteroidetes"), fill="Family", x="station", title = "Family level diversity within the Bacteroidetes") +
    theme_hc()  #+ theme(legend.position = "none")
#as_ggplot(get_legend(plot_bar(subset_taxa(bac_ra_genus, Phylum == "Bacteroidetes"), fill="Genus", x="station", title = "Diversity at Class level")))
ggsave(file="bacteroidetes_diversity.svg", width=10, height=8)

plot_bar(subset_taxa(bac_ra_genus, Phylum == "Actinobacteria"), fill="Genus", x="station", title = "Genera level diversity within the Actinobacteria") +
    theme_hc()  #+ theme(legend.position = "none")
#as_ggplot(get_legend(plot_bar(subset_taxa(bac_ra_genus, Phylum == "Bacteroidetes"), fill="Genus", x="station", title = "Diversity at Class level")))
ggsave(file="actinobacteria_diversity.svg", width=10, height=8)

plot_bar(subset_taxa(bac_ra_genus, Phylum == "Firmicutes"), fill="Genus", x="station", title = "Genera level diversity within the Firmicutes") +
    theme_hc()  #+ theme(legend.position = "none")
#as_ggplot(get_legend(plot_bar(subset_taxa(bac_ra_genus, Phylum == "Bacteroidetes"), fill="Genus", x="station", title = "Diversity at Class level")))
ggsave(file="firmicutes_diversity.svg", width=10, height=8)

# Alpha diversity estimates
plot_richness(bac_ra, x="area", measures =c("Simpson"), color="area") + 
   geom_boxplot() +
   geom_jitter() +
   theme_hc()
ggsave(file="alpha_diversity.svg", width=10, height=8)

# Subset specific taxa involved in selected metabolisms for downstream analysis
prok_alk <- subset_taxa(bac_ra_genus, Genus == "Alcanivorax" | Genus == "Oleispira" 
                       | Genus == "Marinobacter"  | Genus == "Cycloclasticus"
                       | Genus == "Planomicrobium" | Genus == "Yeosuana "
                       | Genus == "Oleiphilus" | Genus == "Thalassolituus"
                    | Genus == "Oleiphilus" | Genus == "Thalassolituus"
                        | Genus == "Methylophaga" | Genus == "Methylobacillus"
                    ) # Subsetting only the olbbligate Alkane oxidizers
# Alkane Oxidizers
plot_bar(prok_alk, fill="Genus", x="station", title="Abundance of alkane oxidizers") +
  theme_hc() # Plot bar of alkane Oxidizing Bacteria
ggsave(file="hydrocarbon_degraders_diversity.svg", width=10, height=8)

### Inkdot plots
# Define the taxa level to be plotted 
bac_pha <- data.frame(otu_table(bac_ra_phyla)) # Phyla
bac_cla <- data.frame(otu_table(bac_ra_class)) # Class
bac_ord <- data.frame(otu_table(bac_ra_order)) # Order
bac_fam <- data.frame(otu_table(bac_ra_family)) # Family
bac_alk <- data.frame(otu_table(prok_alk)) # Alkane

names(prok_sample)

# Define environmental vectors to reorder plot
temp <- as.vector(sample_data(bac_ra_phyla))$temp
sal <- as.vector(sample_data(bac_ra_phyla))$sal
chla <-  as.vector(sample_data(bac_ra_phyla))$chla
fvfm <-  as.vector(sample_data(bac_ra_phyla))$fv_fm
biomass <-  as.vector(sample_data(bac_ra_phyla))$total
diato_ra <-  as.vector(sample_data(bac_ra_phyla))$diato_ra
hapto_ra <-  as.vector(sample_data(bac_ra_phyla))$hapto_ra

# Plot the bubbleplot for selected taxa and vector
#svg("inkspot_ord_temp.svg")
inkspot(bac_ord, temp, cex.max = 8,
          cex.axis = 0.3,
          legend.values =NULL,
          use.rank = T,
          x.axis.top = T,
          spec.names = get_taxa_unique(bac_ra_order, "Order"),
          site.names = data.frame(sample_data(bac_ra_order))$station,
          main = "Orders abundance w/ sites ordered by Temp on top axis"
        ) # Phyla with sites ordered by ph on top axis
#dev.off()

# Plot the bubbleplot for selected taxa and vector
#svg("inkspot_fam_sal.svg")
inkspot(bac_fam, sal, cex.max = 8,
          cex.axis = 0.3,
          legend.values =NULL,
          use.rank = T,
          x.axis.top = T,
          spec.names = get_taxa_unique(bac_ra_family, "Family"),
          site.names = data.frame(sample_data(bac_ra_family))$station,
          main = "Family abundance w/ sites ordered by salinity on top axis"
        )
#dev.off()

# Plot the bubbleplot for selected taxa and vector
#svg("inkspot_ord_chla.svg")
inkspot(bac_ord, chla, cex.max = 8,
          cex.axis = 0.3,
          legend.values =NULL,
          use.rank = T,
          x.axis.top = T,
          spec.names = get_taxa_unique(bac_ra_order, "Order"),
          site.names = data.frame(sample_data(bac_ra_order))$station,
          main = "Orders abundance w/ sites ordered by chl-a on top axis"
        ) 
#dev.off()

# Plot the bubbleplot for selected taxa and vector
#svg("inkspot_ord_fvfm.svg")
inkspot(bac_ord, fvfm, cex.max = 8,
          cex.axis = 0.3,
          legend.values =NULL,
          use.rank = T,
          x.axis.top = T,
          spec.names = get_taxa_unique(bac_ra_order, "Order"),
          site.names = data.frame(sample_data(bac_ra_order))$station,
          main = "Orders abundance w/ sites ordered by Fv/Fm on top axis"
        ) 
#dev.off()

plot_bar(subset_taxa(bac_ra_genus, Order == "Oceanospirillales"), fill="Genus", x="station") +
    theme_bw() 
ggsave(file="oceanospirillales_diversity.svg", width=10, height=8)

plot_bar(subset_taxa(bac_ra_genus, Order == "Alteromonadales"), fill="Genus", x="station") +
    theme_hc() 
ggsave(file="alteromonadales_diversity.svg", width=10, height=8)

plot_bar(subset_taxa(bac_ra_genus, Order == "Cellvibrionales"), fill="Genus", x="station") +
    theme_hc() 
ggsave(file="cellvibrionales_diversity.svg", width=10, height=8)

plot_bar(subset_taxa(bac_ra_genus, Order == "Flavobacteriales"), fill="Genus", x="station") +
    theme_hc() 
ggsave(file="flavobacteriales_diversity.svg", width=10, height=8)

# Trees can also be reduced usin tip_glom before plotting the tree
plot_tree(subset_taxa(bac_ra, Order == "Oceanospirillales"), color = "area", label.tips = "Genus",
          size = "abundance", plot.margin = 0.5, ladderize = TRUE)

set.seed(100)

## PCoA unifrac weighted with relative counts
bac_pcoa_w <- ordinate(bac_ra, method = "PCoA", distance = "unifrac", weighted=T)
evals_w <- bac_pcoa_w$values$Eigenvalues

#svg("pcoa_w.svg", height=8, width=8)
plot_ordination(bac_ra, bac_pcoa_w, type = "samples", color = "area", label="station", title="PCoA weighted Unifrac colored by area") +
  labs(col = "Sampling Area") +
  coord_fixed(sqrt(evals_w[2] / evals_w[1]))
#dev.off()

## PCoA unifrac unweighted with relative counts
bac_pcoa_un <- ordinate(prok_data, method = "PCoA", distance = "unifrac", weighted=F)
evals <- bac_pcoa_un$values$Eigenvalues

#svg("pcoa_unw.svg", height=8, width=8)
plot_ordination(prok_data, bac_pcoa_un, type = "samples", color = "area", label="station",title="PCoA unweighted Unifrac colored by area") +
  labs(col = "Sampling Area") +
  coord_fixed(sqrt(evals[2] / evals[1]))
#dev.off()

## nMDS with Jaccard and Bray-Curtis distance
bac_nmds_j <- ordinate(bac_ra, method = "NMDS", distance = "jaccard", weighted=T, trymax=100)
bac_nmds_bc <- ordinate(bac_ra, method = "NMDS", distance = "bray", weighted=T, trymax=100)

#svg("nmds_J.svg", height=8, width=8)
plot_ordination(bac_ra, bac_nmds_j, color="area", label="station", title="nMDS Jaccard diversity colored by Area") +
theme_bw()
#dev.off()

#svg("nmds_J.svg", height=8, width=8)
plot_ordination(bac_ra, bac_nmds_j, type = "species", color="Class", title="nMDS Jaccard diversity colored by Area") +
theme_bw()
#dev.off()

## Compare relative and absolute abundance wunifrac
grid.arrange(nrow = 2, ncol=2,
             plot_ordination(bac_ra, bac_pcoa_w, type = "samples", color = "area", label="station", title="PCoA weighted Unifrac") +
                 theme(legend.position = "none") +
               coord_fixed(sqrt(evals_w[2] / evals_w[1])),
             plot_ordination(bac_ra, bac_pcoa_un, type = "samples", color = "area",,title="PCoA unweighted Unifrac") +
                 theme(legend.position = "none")+
               coord_fixed(sqrt(evals[2] / evals[1])),
             plot_ordination(bac_ra, bac_nmds_j, color="area", label="station", title="nMDS Jaccard distances") +  theme(legend.position = "none"),
             plot_ordination(bac_ra, bac_nmds_bc, color="area", label="station", title="nMDS Bray-Curtis distances") +   theme(legend.position = "none")
             )

# convert the sample_data() within a phyloseq object to a vegan compatible data object
pssd2veg <- function(physeq) {
  sd <- sample_data(physeq)
  return(as(sd,"data.frame"))
}

# convert the otu_table() within a phyloseq object to a vegan compatible data object
psotu2veg <- function(physeq) {
  OTU <- otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}

# Vector fittin of the environmental variable. Starting with replotting the nMDS in vegan
bac_ra.v<-psotu2veg(bac_ra) # custom function to export phyloseq objects to vegan
bac.v_nmds_j<-metaMDS(bac_ra.v, methods="jaccard", trymax=50) #nMDS with jaccard distances
plot(bac.v_nmds_j, display = "sites", main = "nMDS Jaccard")

genus_names <- get_taxa_unique(bac_ra_genus, "Genus")
plot(bac.v_nmds_j, display = "species", type="n", main = "nMDS Jaccard")
ordilabel(bac.v_nmds_j, dis="sp", lab=genus_names)

save.image()

data.frame(sample_data(prok_data2))

# Vector fitting
environmental <- data.frame(sample_data(prok_data2))[,4:34]
env_fitting <- envfit(bac.v_nmds_j, environmental, perm = 9999, na.rm = T)
env_fitting

#svg("jaccard_vectors.svg", height=8, width=8)
plot(bac.v_nmds_j, display = "sites", main = "nMDS Jaccard")
plot(env_fitting, p.max = 0.05, col = "red")
#dev.off()

## Plotting of site score against selected environmental variables
# Salnity
#svg("dim2_sal.svg", height=8, width=3)
plot(data.frame(sample_data(prok_data2))$sal, scores(bac.v_nmds_j)[,2], main="Plot of the nMDS ordination scores against Salinity")
abline(lsfit(data.frame(sample_data(prok_data2))$sal, scores(bac.v_nmds_j)[,2]))
cor.test(data.frame(sample_data(prok_data2))$sal, scores(bac.v_nmds_j)[,2])
#dev.off()

## Plotting of site score against selected environmental variables
# N/P
#svg("dim2_np.svg", height=8, width=3)
plot(data.frame(sample_data(prok_data2))$n_p, scores(bac.v_nmds_j)[,2], main="Plot of the nMDS ordination scores against N/P")
abline(lsfit(data.frame(sample_data(prok_data2))$n_p, scores(bac.v_nmds_j)[,2]))
cor.test(data.frame(sample_data(prok_data2))$n_p, scores(bac.v_nmds_j)[,2])
#dev.off()

## Plotting of site score against selected environmental variables
# Chl_c2
plot(data.frame(sample_data(prok_data2))$chl_c2, scores(bac.v_nmds_j)[,2], main="Plot of the nMDS ordination scores against Chl-c2")
cor.test(data.frame(sample_data(prok_data2))$chl_c2, scores(bac.v_nmds_j)[,2])

## Plotting of site score against selected environmental variables
# fv/fm
plot(data.frame(sample_data(prok_data))$fv_fm, scores(bac.v_nmds_j)[,2], main="Plot of the nMDS ordination scores against Fv/Fm")
cor.test(data.frame(sample_data(prok_data))$fv_fm, scores(bac.v_nmds_j)[,2])

## Plotting of site score against selected environmental variables
# Latitude
#svg("dim1_long.svg", height=3, width=8)
plot(scores(bac.v_nmds_j)[,1], data.frame(sample_data(prok_data2))$long, main="Plot of the nMDS dim.1  scores against Longitude")
abline(lsfit(scores(bac.v_nmds_j)[,1], data.frame(sample_data(prok_data2))$long))
cor.test(scores(bac.v_nmds_j)[,1], data.frame(sample_data(prok_data2))$long)
#dev.off()

## Plotting of site score against selected environmental variables
# Hapto
#svg("dim1_hapto.svg", height=3, width=8)
plot(scores(bac.v_nmds_j)[,1], data.frame(sample_data(prok_data2))$hapto_ra, main="Plot of the nMDS dim.1  scores against Haptophyte")
abline(lsfit(scores(bac.v_nmds_j)[,1], data.frame(sample_data(prok_data2))$hapto_ra))
cor.test(scores(bac.v_nmds_j)[,1], data.frame(sample_data(prok_data2))$hapto_ra)
#dev.off()

## Plotting of site score against selected environmental variables
#svg("dim1_lat.svg", height=8, width=3)
plot(data.frame(sample_data(prok_data2))$lat, scores(bac.v_nmds_j)[,2], main="Plot of the nMDS dim 2 scores against Latitude")
abline(lsfit(data.frame(sample_data(prok_data2))$lat, scores(bac.v_nmds_j)[,2]))
cor.test(data.frame(sample_data(prok_data2))$lat, scores(bac.v_nmds_j)[,2])
#dev.off()

## Plotting of site score against selected environmental variables
# Chl_c2
plot(data.frame(sample_data(prok_data2))$pras_ra, scores(bac.v_nmds_j)[,1], main="Plot of the nMDS dim. 1 scores against Pras")
cor.test(data.frame(sample_data(prok_data2))$pras_ra, scores(bac.v_nmds_j)[,1])

## Plotting of site score against selected environmental variables
# Chl_c2
plot(data.frame(sample_data(prok_data2))$diato_ra, scores(bac.v_nmds_j)[,1], main="Plot of the nMDS dim. 1 scores against Diato")
cor.test(data.frame(sample_data(prok_data2))$diato_ra, c)

#svg("nmds_jaccard_salinity.svg", height=8, width=8)
plot(bac.v_nmds_j, display = "sites", main = "nMDS Jaccard with Salinity")
with(data.frame(sample_data(prok_data2)), ordisurf(bac.v_nmds_j, sal, add = TRUE, col = "green4"))
#dev.off()

#svg("nmds_jaccard_longitude.svg", height=8, width=8)
plot(bac.v_nmds_j, display = "sites", main = "nMDS Jaccard")
with(data.frame(sample_data(prok_data)), ordisurf(bac.v_nmds_j, long, add = TRUE, col = "red"))
#dev.off()

#svg("nmds_jaccard_haptophytes.svg", height=8, width=8)
plot(bac.v_nmds_j, display = "sites", main = "nMDS Jaccard")
with(data.frame(sample_data(prok_data2)), ordisurf(bac.v_nmds_j, hapto_ra, add = TRUE, col = "blue"))
#dev.off()

plot(bac.v_nmds_j, display = "sites", main = "nMDS Jaccard")
with(data.frame(sample_data(prok_data)), ordisurf(bac.v_nmds_j, pras_ra, add = TRUE, col = "green4"))

plot(bac.v_nmds_j, display = "sites", main = "nMDS Jaccard")
with(data.frame(sample_data(prok_data)), ordisurf(bac.v_nmds_j, fv_fm, add = TRUE, col = "green4"))

#svg("nmds_jaccard_NP.svg", height=8, width=8)
plot(bac.v_nmds_j, display = "sites", main = "nMDS Jaccard")
with(data.frame(sample_data(prok_data2)), ordisurf(bac.v_nmds_j, n_p, add = TRUE, col = "orange"))
#dev.off()

#svg("map_stations.svg", height=8, width=8)
plot(data.frame(sample_data(prok_data2))$lat~data.frame(sample_data(prok_data2))$long)
text(data.frame(sample_data(prok_data2))$lat~data.frame(sample_data(prok_data2))$long, labels = data.frame(sample_data(prok_data2))$station, data=data.frame(sample_data(prok_data2)),pos = 4)
#dev.off()

genus_names <- get_taxa_unique(bac_ra_genus, "Genus")
family_names <- get_taxa_unique(bac_ra_genus, "Family")
plot(bac.v_nmds_j, display = "species", type="n", main = "nMDS Jaccard")
ordilabel(bac.v_nmds_j, dis="sp", lab=genus_names)

# Functions for half-violin plots

"%||%" <- function(a, b) {
  if (!is.null(a)) a else b
}

geom_flat_violin <- function(mapping = NULL, data = NULL, stat = "ydensity",
                        position = "dodge", trim = TRUE, scale = "area",
                        show.legend = NA, inherit.aes = TRUE, fill=Color, ...) {
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomFlatViolin,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      trim = trim,
      scale = scale,
      ...
    )
  )
}

GeomFlatViolin <-
  ggproto("GeomFlatViolin", Geom,
          setup_data = function(data, params) {
            data$width <- data$width %||%
              params$width %||% (resolution(data$x, FALSE) * 0.9)

            # ymin, ymax, xmin, and xmax define the bounding rectangle for each group
            data %>%
              group_by(group) %>%
              mutate(ymin = min(y),
                     ymax = max(y),
                     xmin = x - width / 2,
                     xmax = x)
          },

          draw_group = function(data, panel_scales, coord) {
            # Find the points for the line to go all the way around
            data <- transform(data, 
                              xmaxv = x,
                              xminv = x + violinwidth * (xmin - x))

            # Make sure it's sorted properly to draw the outline
            newdata <- rbind(plyr::arrange(transform(data, x = xminv), y),
                             plyr::arrange(transform(data, x = xmaxv), -y))

            # Close the polygon: set first and last point the same
            # Needed for coord_polar and such
            newdata <- rbind(newdata, newdata[1,])

            ggplot2:::ggname("geom_flat_violin", GeomPolygon$draw_panel(newdata, panel_scales, coord))
          },

          draw_key = draw_key_polygon,

          default_aes = aes(weight = 1, colour = "grey20", fill = "white", size = 0.5,
                            alpha = NA, linetype = "solid"),

          required_aes = c("x", "y")
)

#Aurantivirga Polaribacter Polaribacter_1 SAR92_clade
#SUP05_cluster Pseudoalteromonas OM43_clade Clade_Ia 
#Psychrobacter Profundimonas Marinobacter Brumimicrobium 
# Halomonas Colwellia Oleispira OM60(NOR5)_clade
#Alcanivorax Bacteroides Fusobacterium

plot_abundance = function(prok_data2,title = "",
Facet = "Genus", Color = "area"){
p1f = subset_taxa(prok_data2, Genus %in% c("Fusobacterium") | Genus %in% c("Polaribacter")|
                 Genus %in% c("Polaribacter_1")|
                  Genus %in% c("SAR92_clade")|
                  Genus %in% c("SUP05_cluster")|
                  Genus %in% c("Pseudoalteromonas")|
                  Genus %in% c("OM43_clade")|
                  Genus %in% c("Clade_Ia")|
                  Genus %in% c("Profundimonas")|
                  Genus %in% c("Marinobacter")|
                  Genus %in% c("Brumimicrobium")|
                  Genus %in% c("Halomonas")|
                  Genus %in% c("Colwellia")|
                  Genus %in% c("Oleispira")|
                  Genus %in% c("OM60(NOR5)_clade")|
                  Genus %in% c("Alcanivorax")|
                  Genus %in% c("Bacteroides")|
                  Genus %in% c("Fusobacterium")               
                 )
mphyseq = psmelt(p1f)
mphyseq <- subset(mphyseq, Abundance > 0)
    
ggplot(data = mphyseq, mapping = aes_string(x = "area",y = "Abundance",
fill = Color)) +
geom_flat_violin(fill=Color, trim=F) +
stat_summary(fun.y=mean, geom="point", shape=3, size=8) +
coord_flip() +
geom_point(size = 1, alpha = 0.8) +
facet_wrap(facets = Facet) + scale_y_log10() +
theme_hc()
}

plot_abundance(bac_ra,"")
#ggsave("Fusobacterium_genus.svg", height=14, width=18)

#Aurantivirga Polaribacter Polaribacter_1 SAR92_clade
#SUP05_cluster Pseudoalteromonas OM43_clade Clade_Ia 
#Psychrobacter Profundimonas Marinobacter Brumimicrobium 
# Halomonas Colwellia Oleispira OM60(NOR5)_clade
#Alcanivorax Bacteroides Fusobacterium

plot_abundance = function(bac_ra,title = "",
Facet = "Genus", Color = "area"){
p1f = subset_taxa(bac_ra, Genus %in% c("Fusobacterium") | Genus %in% c("Polaribacter")|
                 Genus %in% c("Polaribacter_1")|
                  Genus %in% c("SAR92_clade")|
                  Genus %in% c("SUP05_cluster")|
                  Genus %in% c("Pseudoalteromonas")|
                  Genus %in% c("OM43_clade")|
                  Genus %in% c("Clade_Ia")|
                  Genus %in% c("Profundimonas")|
                  Genus %in% c("Marinobacter")|
                  Genus %in% c("Brumimicrobium")|
                  Genus %in% c("Halomonas")|
                  Genus %in% c("Colwellia")|
                  Genus %in% c("Oleispira")|
                  Genus %in% c("OM60(NOR5)_clade")|
                  Genus %in% c("Alcanivorax")|
                  Genus %in% c("Bacteroides")|
                  Genus %in% c("Fusobacterium")               
                 )
mphyseq = psmelt(p1f)
mphyseq <- subset(mphyseq, Abundance > 0)
    
ggplot(data = mphyseq, mapping = aes_string(x = "area",y = "Abundance",
fill = Color)) +
geom_flat_violin(fill=Color, trim=F) +
stat_summary(fun.y=mean, geom="point", shape=3, size=8) +
coord_flip() +
geom_point(size = 1, alpha = 0.8) +
facet_wrap(facets = Facet) + scale_y_log10() +
theme_hc()
}

plot_abundance(bac_ra,"")
#ggsave("Fusobacterium_genus.svg", height=14, width=18)

# Cryomorphaceae Nitrospinaceae

plot_abundance = function(prok_data2,title = "",
Facet = "Genus", Color = "area"){
p1f = subset_taxa(prok_data2, Genus %in% c("Aurantivirga")            
                 )
mphyseq = psmelt(p1f)
mphyseq <- subset(mphyseq, Abundance > 0)
    
ggplot(data = mphyseq, mapping = aes_string(x = "area",y = "Abundance",
fill = Color)) +
geom_flat_violin(fill=Color, trim=F) +
stat_summary(fun.y=mean, geom="point", shape=3, size=8) +
coord_flip() +
geom_point(size = 1, alpha = 0.8) +
facet_wrap(facets = Facet) + scale_y_log10() +
theme_hc()
}

plot_abundance(bac_ra,"")
ggsave("Aurantivirga_genus.svg", height=4, width=4)

# Cryomorphaceae Nitrospinaceae

plot_abundance = function(prok_data2,title = "",
Facet = "Family", Color = "area"){
p1f = subset_taxa(prok_data2, Family %in% c("Cryomorphaceae") | Family %in% c("Nitrospinaceae")             
                 )
mphyseq = psmelt(p1f)
mphyseq <- subset(mphyseq, Abundance > 0)
    
ggplot(data = mphyseq, mapping = aes_string(x = "area",y = "Abundance",
fill = Color)) +
geom_flat_violin(fill=Color, trim=F) +
stat_summary(fun.y=mean, geom="point", shape=3, size=8) +
coord_flip() +
geom_point(size = 1, alpha = 0.8) +
facet_wrap(facets = Facet) + scale_y_log10() +
theme_hc()
}

plot_abundance(bac_ra,"")
ggsave("Differential_families.svg", height=4, width=4)

plot_abundance = function(prok_data2,title = "",
Facet = "area", Color = "area"){
# Arbitrary subset, based on Phylum, for plotting
p1f = subset_taxa(prok_data2, Genus %in% c("Polaribacter_1"))
mphyseq = psmelt(p1f)
mphyseq <- subset(mphyseq, Abundance > 0)
    
ggplot(data = mphyseq, mapping = aes_string(x = "Abundance",
color = Color, fill = Color)) +
geom_density(kernel="gaussian") +
scale_x_log10()   
}

plot_abundance(bac_ra,"")

plot_density(bac_ra_genus, variable = "TGAGGAATATTGGACAATGGAGGAGACTCTGATCCAGCCATGCCGCGTGTAGGAAGAATGCCCTATGGGTTGTAAACTACTTTTATACAGGAAGAAACACTGGTATGTATACCAGCTTGACGGTACTGTAAGAATAAGGACCGGCTAACTCCGTG")

bac_ra_genus
tax_table(bac_ra_genus)
write.csv(as.matrix(tax_table(bac_ra_genus)), "ant17_genera_list.csv")

#Merge the phyloseq object based on the area parameters 
prok_data2_merged <- merge_samples(bac_ra,"area")

prok_data2_merged

# Subset all the taxa that appear in both the areas (core microbiome)
prok_data2_core = filter_taxa(prok_data2_merged, function(x) sum(x >= 1) == (2), TRUE)
prok_data2_core

# See how many ASV are in each sample
asv_df <- t(otu_table(prok_data2_merged))
colSums(asv_df != 0)

# Install and load the  Venn diagram package
#install.packages('VennDiagram')
library(VennDiagram)

# Plot Venn diagram
#svg("venn_asv_area.svg",height=4, width=4)
draw.pairwise.venn(area1 = (497), area2 = (571), cross.area = 353, category = c("Terranova Bay", 
    "Ross Sea"))
#dev.off()

# Kruskal-Wallis non-parametric test of one-way anova
p1f = subset_taxa(prok_data2, Family %in% c("Nitrospinaceae"))
mphyseq = psmelt(p1f)
mphyseq <- as.data.frame(subset(mphyseq, Abundance > 0))
test <- kruskal.test(mphyseq$Abundance~mphyseq$area)
#table_kruskal_test <- data.frame(test$statistic, test$p.value)
#names(table_kruskal_test) <- c("chi-sq","p-value")
table_kruskal_test[nrow(table_kruskal_test) + 1,] = list(test$statistic, test$p.value)
row.names(table_kruskal_test)<- c("Aurantivirga","Polaribacter","Polaribacter_1","SAR92_clade","SUP05_cluster","Pseudoalteromonas",
                                 "OM43_clade","Clade_Ia","Psychrobacter","Profundimonas","Marinobacter","Brumimicrobium","Halomonas",
                                 "Colwellia","Oleispira","OM60(NOR5)_clade","Alcanivorax","Bacteroides","Fusobacterium","Cryomorphaceae",
                                 "Nitrospinaceae")
table_kruskal_test



# Extract OTU table and Taxonomy tables from the fluid samples object, both absolute and relative
# Absolute counts
bac_form <- microbiomeutilities::format_to_besthit(bac_ra)
bac.otu <- t(otu_table(bac_form)@.Data)
bac.tax <- as.data.frame(tax_table(bac_form)@.Data)

#################################################################################################################
## Network from absolute count data
# Manual correlation construction test Spearman > 0.65
bac.cor_p6 <- cor(bac.otu, method="spearman")
bac.cor_p6[bac.cor_p6 < 0.65] = 0
bac.cor_p6.ig <- graph.adjacency(bac.cor_p6, mode='undirected', add.rownames = TRUE, weighted = TRUE)

# Set Phyla coloring
bac.cor.phyla <- map_levels(colnames(bac.otu), from = "best_hit", to = "Phylum", tax_table(bac_form))
bac.cor.family <- map_levels(colnames(bac.otu), from = "best_hit", to = "Family", tax_table(bac_form))

# Set up the colors for the phyla
c31 <- c("dodgerblue2","#E31A1C", "red", "green4", "#6A3D9A", "purple",
         "#FF7F00", "orange", "black","gold1", "skyblue2","#FB9A99", "pink",
         "palegreen2", "#CAB2D6", "purple","#FDBF6F", "orange","gray70", "khaki2",
         "maroon","orchid1","deeppink1","blue1","steelblue4","darkturquoise","green1","yellow4","yellow3",
         "darkorange4","brown")
my_color_p=c31[as.numeric(as.factor(bac.cor.phyla))]
my_color_g=c31[as.numeric(as.factor(bac.cor.family))]

# Manual correlation construction test Spearman > 0.55
bac.cor_p55 <- cor(bac.otu, method="spearman")
bac.cor_p55[bac.cor_p55 < 0.55] = 0
bac.cor_p55.ig <- graph.adjacency(bac.cor_p55, mode='undirected', add.rownames = TRUE, weighted = TRUE)
# Plot the computed netwrok using the FR force-directed layout
#svg("network_055_spearman.svg",height=8, width=8)
plot(igraph::simplify(bac.cor_p55.ig), layout=layout_with_fr, vertex.color=my_color, vertex.size = 4, vertex.label = NA)
legend("left", legend=levels(as.factor(bac.cor.phyla)) , col = c31 , bty = "n", pch=20 , pt.cex = 1, cex = 0.5, text.col="black", horiz = FALSE, inset = c(0.1, 0.1))
#dev.off()

# Manual correlation construction test Spearman > 0.6
bac.cor_p6 <- cor(bac.otu, method="spearman")
bac.cor_p6[bac.cor_p6 < 0.6] = 0
bac.cor_p6.ig <- graph.adjacency(bac.cor_p6, mode='undirected', add.rownames = TRUE, weighted = TRUE)
# Plot the computed netwrok using the FR force-directed layout
#svg("network_06_spearman.svg",height=8, width=8)
plot(igraph::simplify(bac.cor_p6.ig), layout=layout_with_fr, vertex.color=my_color, vertex.size = 4, vertex.label = NA)
legend("left", legend=levels(as.factor(bac.cor.phyla)) , col = c31 , bty = "n", pch=20 , pt.cex = 1, cex = 0.5, text.col="black", horiz = FALSE, inset = c(0.1, 0.1))
#dev.off()

# Manual correlation construction test Spearman > 0.65
bac.cor_p65 <- cor(bac.otu, method="spearman")
bac.cor_p65[bac.cor_p65 < 0.65] = 0
bac.cor_p65.ig <- graph.adjacency(bac.cor_p65, mode='undirected', add.rownames = TRUE, weighted = TRUE)
# Plot the computed netwrok using the FR force-directed layout
#svg("network_065_spearman.svg",height=8, width=8)
plot(igraph::simplify(bac.cor_p65.ig), layout=layout_with_fr, vertex.color=my_color, vertex.size = 4, vertex.label = NA)
legend("left", legend=levels(as.factor(bac.cor.phyla)) , col = c31 , bty = "n", pch=20 , pt.cex = 1, cex = 0.5, text.col="black", horiz = FALSE, inset = c(0.1, 0.1))
#dev.off()

# Manual correlation construction test Spearman > 0.7
bac.cor_p7 <- cor(bac.otu, method="spearman")
bac.cor_p7[bac.cor_p7 < 0.7] = 0
bac.cor_p7.ig <- graph.adjacency(bac.cor_p7, mode='undirected', add.rownames = TRUE, weighted = TRUE)
# Plot the computed netwrok using the FR force-directed layout
#svg("network_07_spearman.svg",height=8, width=8)
plot(igraph::simplify(bac.cor_p7.ig), layout=layout_with_fr, vertex.color=my_color, vertex.size = 4, vertex.label = NA)
legend("left", legend=levels(as.factor(bac.cor.phyla)) , col = c31 , bty = "n", pch=20 , pt.cex = 1, cex = 0.5, text.col="black", horiz = FALSE, inset = c(0.1, 0.1))
#dev.off()

# Manual correlation construction test Spearman > 0.75
bac.cor_p75 <- cor(bac.otu, method="spearman")
bac.cor_p75[bac.cor_p75 < 0.75] = 0
bac.cor_p75.ig <- graph.adjacency(bac.cor_p75, mode='undirected', add.rownames = TRUE, weighted = TRUE)
# Plot the computed netwrok using the FR force-directed layout
#svg("network_075_spearman.svg",height=8, width=8)
plot(igraph::simplify(bac.cor_p75.ig), layout=layout_with_fr, vertex.color=my_color, vertex.size = 4, vertex.label = NA)
legend("left", legend=levels(as.factor(bac.cor.phyla)) , col = c31 , bty = "n", pch=20 , pt.cex = 1, cex = 0.5, text.col="black", horiz = FALSE, inset = c(0.1, 0.1))
#dev.off()

# Manual correlation construction test Spearman > 0.8
bac.cor_p8 <- cor(bac.otu, method="spearman")
bac.cor_p8[bac.cor_p8 < 0.8] = 0
bac.cor_p8.ig <- graph.adjacency(bac.cor_p8, mode='undirected', add.rownames = TRUE, weighted = TRUE)
# Plot the computed netwrok using the FR force-directed layout
#svg("network_08_spearman.svg",height=8, width=8)
plot(igraph::simplify(bac.cor_p8.ig), layout=layout_with_fr, vertex.color=my_color, vertex.size = 4, vertex.label = NA)
legend("left", legend=levels(as.factor(bac.cor.phyla)) , col = c31 , bty = "n", pch=20 , pt.cex = 1, cex = 0.5, text.col="black", horiz = FALSE, inset = c(0.1, 0.1))
#dev.off()

# Manual correlation construction test Spearman > 0.85
bac.cor_p85 <- cor(bac.otu, method="spearman")
bac.cor_p85[bac.cor_p85 < 0.85] = 0
bac.cor_p85.ig <- graph.adjacency(bac.cor_p85, mode='undirected', add.rownames = TRUE, weighted = TRUE)
# Plot the computed netwrok using the FR force-directed layout
svg("network_085_spearman.svg",height=8, width=8)
plot(igraph::simplify(bac.cor_p85.ig), layout=layout_with_fr, vertex.color=my_color, vertex.size = 4, vertex.label = NA)
legend("left", legend=levels(as.factor(bac.cor.phyla)) , col = c31 , bty = "n", pch=20 , pt.cex = 1, cex = 0.5, text.col="black", horiz = FALSE, inset = c(0.1, 0.1))
dev.off()

save.image()

colnames(sample_data(prok_data2))

env_data_spearman <- cor(as.data.frame(sample_data(prok_data2))[,c(4:42)], method = "spearman", use = "complete.obs") 
env_data_pearson <- cor(as.data.frame(sample_data(prok_data2))[,c(4:42)], use = "complete.obs")

#svg("env_heatmap_pearson.svg", height=8, width=8)
heatmap(env_data_spearman, symm=T, revC=F)
#dev.off()

#svg("env_heatmap_pearson.svg", height=8, width=8)
heatmap(env_data_pearson, symm=T, revC=F)
#dev.off()

env_data_pearson

env_data<-as.data.frame(sample_data(prok_data2))

for (i in 4:42) {
    test.i <- cor(scores(bac.v_nmds_j)[,1], env_data[,i], use="complete.obs")
    print(test.i)
     }

for (i in 4:42) {
    test.i <- cor(scores(bac.v_nmds_j)[,2], env_data[,i], use="complete.obs")
    print(test.i)
     }

bac_ra_otu <- t(as.matrix(otu_table(bac_ra)))

bac_ra_tax <- as.matrix(tax_table(bac_ra))

head(bac_ra_otu)

head(bac_ra_tax)

bac_ra_otutax <- cbind(bac_ra_otu,bac_ra_tax)

write.csv(bac_ra_otutax, "otu_table_w_taxa.csv")

save.image()


