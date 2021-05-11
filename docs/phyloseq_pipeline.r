## Basci phyloseq pyypeline to get the 16S rRNA analysis started
## Giovannelli Lab - Jun 2020

## This pipeline assumes you ave already preprocessed you 16S rRNA reads
## using dada2 or any other tool of your choice. Keep in mind that the
## import procedure to create the basic phyloseq object is designed assuming
## you have used our basic dada2_pipeline. Check out our github to see it.
## It has been checked using phyloseq v.1.30.0

## Not all the plots or analyses reported here will make sense for your data,
## especially those requiring informations about factors or specific environmental
## variables. We give here some basic examples, modify them according to your
## specific research needs

## Get in touch for any question!!! dgiovannelli.github.io

### Load required libraries
library(microbiome) # data analysis and visualisation
library(phyloseq) # also the basis of data object. Data analysis and visualisation
library(ggthemes) # additional themes fro ggplot2
#library(network) # networks
#library(intergraph)  # networks
#library(ggnet)   # network plotting with ggplot
library(igraph)  # networks
library(vegan)
library(ggpubr)

## Verify that your accompaning environmental file has AT LEAST the following:
# - the first column contains the same names of the columns of your otu_table.
# the easiest way is to read the column name of the otu_table and use it to
# create you env_dataset file in excel
# - you have a column containing the name of your samples/station/site
# - you have a column clearly marking the blanks and the samples in your file.
# for this pipeline we assume that you have flagged your blanks and sample with
# these exact names
# - one or more factors meaningful for your experimentaal design (something
# like area, treatment, group, etc...). We use here for the pipeline factor1
# and factor2, of two levels each and with factor2 nested into factor1. If you
# do not know what i'm talking about, go read https://en.wikipedia.org/wiki/Design_of_experiments
# and https://en.wikipedia.org/wiki/Factorial_experiment
# - you have at least one other measured variables. Here we use temperature (temp)
# and pH as examples for the pipeline

## Remember to save periodically your work using the save.image() function. And
# remember that each time you resume the analysis, beside calling back your session
# with the load(".RData") command you need to reload the libraries.

## Also, remember that you can save any plot by using ggsave() or the
# svg("name.svg", height=8, width=8) and dev.off() combination. Google it for more help.

# import your environmental data file
env_dataset <- read.csv("/home/minion/Desktop/pipeline_test/phyloseq_test/env_data.csv", header=T, sep=",", row.names=1)

# check the imported file to see if factors have been impoted correctly
summary(env_dataset)

# Create the phyloseq object
prok_data_raw <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = F),
                      sample_data(env_dataset),
                      tax_table(taxa))

# Pre-clean up and normalization statistics
prok_data_raw
sum(readcount(prok_data_raw))
readcount(prok_data_raw)
sample_names(prok_data_raw)

# Checking diversity in  the blank
plot_bar(subset_samples(prok_data_raw, sample_names(prok_data_raw) == "blk"), fill="Phylum", title = "Diversity in the blank") +
theme_hc() +
theme(legend.position = "bottom")

# Define a function to remove negative controls
prune_negatives = function(physeq, negs, samps) {
  negs.n1 = prune_taxa(taxa_sums(negs)>=1, negs)
  samps.n1 = prune_taxa(taxa_sums(samps)>=1, samps)
  allTaxa <- names(sort(taxa_sums(physeq),TRUE))
  negtaxa <- names(sort(taxa_sums(negs.n1),TRUE))
  taxa.noneg <- allTaxa[!(allTaxa %in% negtaxa)]
  return(prune_taxa(taxa.noneg,samps.n1))
}

#Removing the ASV found in the negative controls from all the samples
blanks = subset_samples(prok_data_raw, sample_names(prok_data_raw) == "blk")
samples = subset_samples(prok_data_raw, sample_names(prok_data_raw) != "blk")
prok_data = prune_negatives(prok_data_raw,blanks,samples)

# Check the stats after removing the blanks
prok_data
sum(readcount(prok_data))
(sum(readcount(prok_data))/sum(readcount(prok_data_raw)))*100 # check percentage of reads after blank removal
readcount(prok_data)

# Clean up unwanted sequences from Eukarya, mitochrondria and chloroplast
prok_data <- subset_taxa(prok_data,  (Kingdom != "Eukaryota") | is.na(Kingdom))
prok_data <- subset_taxa(prok_data, (Order!="Chloroplast") | is.na(Order))
prok_data <- subset_taxa(prok_data, (Family!="Mitochondria") | is.na(Family))
prok_data
sum(readcount(prok_data))

# Removing the potential human pathogens and contaminants
# This step needs to be evaluated with attention since many of these genera
# might be relevant in many environmental settings. By default this step is
# commented out. Feel free to experiment and evaluate the results.
#prok_data <- subset_taxa(prok_data,  (Genus != "Lactococcus") | is.na(Genus))
#prok_data <- subset_taxa(prok_data,  (Genus != "Lactobacillus") | is.na(Genus))
#prok_data <- subset_taxa(prok_data,  (Genus != "Cutibacterium") | is.na(Genus))
#prok_data <- subset_taxa(prok_data,  (Genus != "Enterococcus") | is.na(Genus))
#prok_data <- subset_taxa(prok_data,  (Genus != "Streptococcus") | is.na(Genus))
#prok_data <- subset_taxa(prok_data,  (Genus != "Acinetobacter") | is.na(Genus))
#prok_data <- subset_taxa(prok_data,  (Genus != "Citrobacter") | is.na(Genus))
#prok_data <- subset_taxa(prok_data,  (Genus != "Bifidobacterium") | is.na(Genus))
#prok_data <- subset_taxa(prok_data,  (Genus != "Proprionibacterium") | is.na(Genus))
#prok_data <- subset_taxa(prok_data,  (Order != "Enterobacteriales") | is.na(Order))
#prok_data <- subset_taxa(prok_data,  (Genus != "Corynebacterium_1") | is.na(Genus))
#prok_data <- subset_taxa(prok_data,  (Genus != "Corynebacterium") | is.na(Genus))
#prok_data <- subset_taxa(prok_data,  (Genus != "Escherichia/Shigella") | is.na(Genus))
#prok_data <- subset_taxa(prok_data, (Genus != "Abiotrophia") |  is.na(Genus))
#prok_data <- subset_taxa(prok_data, (Genus != "Achromobacter") | is.na(Genus))
#prok_data <- subset_taxa(prok_data, (Genus != "Actinobacillus") | is.na(Genus))
#prok_data <- subset_taxa(prok_data, (Genus != "Arcanobacterium") | is.na(Genus))
#prok_data <- subset_taxa(prok_data, (Genus != "Babesia") | is.na(Genus))
#prok_data <- subset_taxa(prok_data, (Genus != "Bacillus") | is.na(Genus))
#prok_data <- subset_taxa(prok_data, (Genus != "Bartonella") | is.na(Genus))
#prok_data <- subset_taxa(prok_data, (Genus != "Bordetella") |  is.na(Genus))
#prok_data <- subset_taxa(prok_data, (Genus != "Borrelia") | is.na(Genus))
#prok_data <- subset_taxa(prok_data, (Genus != "Brodetella") |  is.na(Genus))
#prok_data <- subset_taxa(prok_data, (Genus != "Brucella") | is.na(Genus))
#prok_data <- subset_taxa(prok_data, (Genus != "Capnocytophaga") | is.na(Genus))
#prok_data <- subset_taxa(prok_data, (Genus != "Chlamydia") | is.na(Genus))
#prok_data <- subset_taxa(prok_data, (Genus != "Comamonas") | is.na(Genus))
#prok_data <- subset_taxa(prok_data, (Genus != "Coxiella") | is.na(Genus))
#prok_data <- subset_taxa(prok_data, (Genus != "Cronobacter") | is.na(Genus))
#prok_data <- subset_taxa(prok_data, (Genus != "Deinococcus") | is.na(Genus))
#prok_data <- subset_taxa(prok_data, (Genus != "Dermatophilus") | is.na(Genus))
#prok_data <- subset_taxa(prok_data, (Genus != "Ehrlichia") | is.na(Genus))
#prok_data <- subset_taxa(prok_data, (Genus != "Enterococcus") | is.na(Genus))
#prok_data <- subset_taxa(prok_data, (Genus != "Erysipelothrix") |  is.na(Genus))
#prok_data <- subset_taxa(prok_data, (Genus != "Escherichia") | is.na(Genus))
#prok_data <- subset_taxa(prok_data, (Genus != "Francisella") | is.na(Genus))
#prok_data <- subset_taxa(prok_data, (Genus != "Gardnerella") | is.na(Genus))
#prok_data <- subset_taxa(prok_data, (Genus != "Granulicatella") | is.na(Genus))
#prok_data <- subset_taxa(prok_data, (Genus != "Haemophilus") | is.na(Genus))
#prok_data <- subset_taxa(prok_data, (Genus != "Hafnia") | is.na(Genus))
#prok_data <- subset_taxa(prok_data, (Genus != "Helicobacter") | is.na(Genus))
#prok_data <- subset_taxa(prok_data, (Genus != "Klebsiella") | is.na(Genus))
#prok_data <- subset_taxa(prok_data, (Genus != "Kocuria") | is.na(Genus))
#prok_data <- subset_taxa(prok_data, (Genus != "Lawsonia") | is.na(Genus))
#prok_data <- subset_taxa(prok_data, (Genus != "Legionella") | is.na(Genus))
#prok_data <- subset_taxa(prok_data, (Genus != "Leptospira") | is.na(Genus))
#prok_data <- subset_taxa(prok_data, (Genus != "Listeria") | is.na(Genus))
#prok_data <- subset_taxa(prok_data, (Genus != "Merkel_cell") | is.na(Genus))
#prok_data <- subset_taxa(prok_data, (Genus != "Micrococcus") | is.na(Genus))
#prok_data <- subset_taxa(prok_data, (Genus != "Morganella") | is.na(Genus))
#prok_data <- subset_taxa(prok_data, (Genus != "Mycoplasma") | is.na(Genus))
#prok_data <- subset_taxa(prok_data, (Genus != "Neisseria") | is.na(Genus))
#prok_data <- subset_taxa(prok_data, (Genus != "Nocardia") | is.na(Genus))
#prok_data <- subset_taxa(prok_data, (Genus != "Pasteurella") | is.na(Genus))
#prok_data <- subset_taxa(prok_data, (Genus != "Plesiomonas") | is.na(Genus))
#prok_data <- subset_taxa(prok_data, (Genus != "Propionibacterium") | is.na(Genus))
#prok_data <- subset_taxa(prok_data, (Genus != "Proteus") | is.na(Genus))
#prok_data <- subset_taxa(prok_data, (Genus != "Providencia") | is.na(Genus))
#prok_data <- subset_taxa(prok_data, (Genus != "Rothia") | is.na(Genus))
#prok_data <- subset_taxa(prok_data, (Genus != "Salmonella") | is.na(Genus))
#prok_data <- subset_taxa(prok_data, (Genus != "Serratia") | is.na(Genus))
#prok_data <- subset_taxa(prok_data, (Genus != "Shewanella") | is.na(Genus))
#prok_data <- subset_taxa(prok_data, (Genus != "Shigella") | is.na(Genus))
#prok_data <- subset_taxa(prok_data, (Genus != "Sphaerophorus") | is.na(Genus))
#prok_data <- subset_taxa(prok_data, (Genus != "Staphylococcus") | is.na(Genus))
#prok_data <- subset_taxa(prok_data, (Genus != "Stenotrophomonas") | is.na(Genus))
#prok_data <- subset_taxa(prok_data, (Genus != "Streptococcus") | is.na(Genus))
#prok_data <- subset_taxa(prok_data, (Genus != "Treponema") | is.na(Genus))
#prok_data <- subset_taxa(prok_data, (Genus != "Yersinia") | is.na(Genus))

## Check the results of the filtering
# prok_data
# sum(readcount(prok_data))

# Filter ASV with less than 5 reads across all the samples
  prok_data <- prune_taxa(taxa_sums(bac_data) > 5, bac_data)
  prok_data

  # After these filtering step remove the ASV and samples that are left with all zeros
  prok_data = filter_taxa(prok_data, function(x) sum(x) > 0, TRUE)
  prok_data
  sum(readcount(prok_data))

# Normalize the counts across the different samples by converting the abundance to
# relative abundance and multiply by the median library size
prok_ndata <- transform_sample_counts(prok_data, function(x) ((x / sum(x))*median(readcount(prok_data))))

# Transform normalized abundance to relative abundances for plotting and some stats
prok_ra = transform_sample_counts(prok_ndata, function(x){x / sum(x)})

save.image() # save the work done so far!

# from here on prok_ndata and prok_ra are the two main objects to be used for
# downstream analysis

## Agglomerate at a specific taxonomic level at the Genus level
prok_ra_genus = tax_glom(prok_ra, "Genus", NArm = TRUE)
prok_ra_family = tax_glom(prok_ra, "Family", NArm = TRUE)
prok_ra_order = tax_glom(prok_ra, "Order", NArm = TRUE)
prok_ra_class = tax_glom(prok_ra, "Class", NArm = TRUE)
prok_ra_phyla = tax_glom(prok_ra, "Phylum", NArm = TRUE)
prok_ra_kingdom = tax_glom(prok_ra, "Kingdom", NArm = TRUE)

# Plot Relative Abundance by Class for each station
plot_bar(prok_ra_phyla, fill="Phylum", title = "Diversity at Phylum level") +
theme_hc()

# Plot Relative Abundance by Class for each station divide by factor1
plot_bar(prok_ra_class, fill="Class", x="station", title = "Diversity at Phylum level by factor1") +
facet_wrap(~depth) +
theme_hc()

# Alpha diversity estimates by factor 1
plot_richness(prok_ndata, x="station", measures =c("Simpson"), color="depth") +
   geom_boxplot() +
   geom_jitter() +
   theme_hc()

### VERIFY FUNCTIONALITY
# Extract the Top 10 Phyla by relative abundance
top10_phyla <- sort(taxa_sums(prok_ra_phyla), TRUE)[1:10]
phyla_top10 <- prune_taxa(names(top10_phyla), prok_ra_phyla)

# Plot the top 10 families
plot_bar(phyla_top10, fill="Phylum", x="station", title = "Diversity of the top 10 phyla") +
    theme_hc()

# List the top 10 phyla. The operation can be repeated for any N at any taxonomic levels
get_taxa_unique(phyla_top10, "Phylum")

save.image(file="/home/minion/Desktop/pipeline_test/phyloseq_test/RData/phyla_top10.RData")

# Plotting the diversity of one of the top phyla. The Proteobacteria are used as example

##### prok_ra instead of prok_ra_phyla
##### to use as_ggplot → library(ggpubr) is necessary
plot_bar(subset_taxa(prok_ra, Phylum == "Proteobacteria"), fill="Order", x="station", title = "Order level diversity within the Proteobacteria") +
    theme_hc() + theme(legend.position = "none")
as_ggplot(get_legend(plot_bar(subset_taxa(prok_ra, Phylum == "Proteobacteria"), fill="Order", x="station", title = "Diversity at Class level")))


## Starting the multivariate analysis on the dataset
set.seed(100) # to make the analysis reproducible

## nMDS with Bray curtis weighted and unweighted
prok_nmds_w <- ordinate(prok_ndata, method = "NMDS", distance = "bray", weighted=T, trymax=100)
prok_nmds_uw <- ordinate(prok_ndata, method = "NMDS", distance = "bray", weighted=F, trymax=100)

# Plot the two ordinations
plot_ordination(prok_ndata, prok_nmds_w, color="factor1", label="factor2", title="nMDS weighted Bray similarity colored by factor1") +
theme_bw()

plot_ordination(prok_ndata, prok_nmds_uw, color="factor1", label="factor2", title="nMDS unweighted Bray similarity colored by factor1") +
theme_bw()

## Create two functions for vector fitting in Vegan
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


## Vector fitting of the environmental variable. Starting with replotting the nMDS in vegan
prok_ra.v<-psotu2veg(prok_ndata) # custom function to export phyloseq objects to vegan
prok.v_nmds_w<-metaMDS(prok_ra.v, methods="bray", trymax=50) #nMDS with Bray-Curtis distances
plot(prok.v_nmds_w, display = "sites", main = "nMDS Bray")

# Vector fitting. Assuming temp and pH are variable 4 and 5
environmental <- pssd2veg(prok_ndata)
env_fitting <- envfit(prok.v_nmds_w, environmental[,4:5], perm = 9999, na.rm = T)
env_fitting

plot(prok.v_nmds_w, display = "sites", main = "nMDS Jaccard")
plot(env_fitting, p.max = 0.0499, col = "red")

## Plotting of site score against selected environmental variables
# Temperature and dimension 1
plot(scores(prok.v_nmds_w)[,1], data.frame(sample_data(prok_ndata))$temperature, main="Plot of the nMDS ordination scores against Temperature")
abline(lsfit(scores(prok.v_nmds_w)[,1], data.frame(sample_data(prok_data))$temperature))
cor.test(scores(prok.v_nmds_w)[,1], data.frame(sample_data(prok_data))$temperature)

# Temperature and dimension 2
plot(scores(prok.v_nmds_w)[,2], data.frame(sample_data(prok_ndata))$temperature, main="Plot of the nMDS ordination scores against Temperature")
abline(lsfit(scores(prok.v_nmds_w)[,2], data.frame(sample_data(prok_data))$temperature))
cor.test(scores(prok.v_nmds_w)[,2], data.frame(sample_data(prok_data))$temperature)

# Isoline fitting to check for linearity of the environmental variable with the ordination gradient
plot(prok.v_nmds_w, display = "sites", main = "nMDS Jaccard with Temperature")
with(data.frame(sample_data(prok_ndata)), ordisurf(prok.v_nmds_w, temperature, add = TRUE, col = "red"))

## Testing for the influence of factors in explaining the diversity using ADONIS
adn.test <- adonis(distance(prok_ndata, method="bray") ~ factor1*factor2, data = data.frame(sample_data(prok_ndata)), perm = 999)
adn.test

# Using the Kruskal-Wallis non-parametric test of one-way anova to test for differential abubndance of selected taxa between factor1
p1f = subset_taxa(prok_ndata, Phylum %in% c("Proteobacteria"))
mphyseq = psmelt(p1f)
mphyseq <- as.data.frame(subset(mphyseq, Abundance > 0))
kruskal.test(mphyseq$Abundance~mphyseq$factor1)


#################################################################################################################
## Network from absolute count data

# Extract OTU table and Taxonomy tables from the dataset for network construction
prok_form <- microbiomeutilities::format_to_besthit(prok_ndata)
prok.otu <- t(otu_table(prok_form)@.Data)
prok.tax <- as.data.frame(tax_table(prok_form)@.Data)

# You can test other correlation level cutoff beside 0.65. Beware of caveats in interpreting the resulting networks
prok.cor_p6 <- cor(prok.otu, method="spearman")
prok.cor_p6[prok.cor_p6 < 0.65] = 0
prok.cor_p6.ig <- graph.adjacency(prok.cor_p6, mode='undirected', add.rownames = TRUE, weighted = TRUE)

# Set Phyla coloring
prok.cor.phyla <- map_levels(colnames(prok.otu), from = "best_hit", to = "Phylum", tax_table(prok_form))

# Set up the colors for the phyla
c31 <- c("dodgerblue2","#E31A1C", "red", "green4", "#6A3D9A", "purple",
         "#FF7F00", "orange", "black","gold1", "skyblue2","#FB9A99", "pink",
         "palegreen2", "#CAB2D6", "purple","#FDBF6F", "orange","gray70", "khaki2",
         "maroon","orchid1","deeppink1","blue1","steelblue4","darkturquoise","green1","yellow4","yellow3",
         "darkorange4","brown")
my_color_p=c31[as.numeric(as.factor(prok.cor.phyla))]

# Plot the computed netwrok using the FR force-directed layout
plot(igraph::simplify(prok.cor_p6.ig), layout=layout_with_fr, vertex.color=my_color_p, vertex.size = 4, vertex.label = NA)
legend("left", legend=levels(as.factor(prok.cor.phyla)) , col = c31 , bty = "n", pch=20 , pt.cex = 1, cex = 0.5, text.col="black", horiz = FALSE, inset = c(0.1, 0.1))

# Construct a matrix of environmental covariates. For only two environmental variables, like in htis example, the results
# are not really meaningful. For larger dataset this can be important in interpreting results
env_data_spearman <- cor(as.data.frame(sample_data(prok_ndata))[,c(10:11)], method = "spearman", use = "complete.obs")
heatmap(env_data_spearman, symm=T, revC=F)

# After these first analysis, the specific of the next steps will depend upon the scientific questions.
# Take a look at published papers, read through the aanalysis we have published on our GitHub for
# other projects, try to replicate some of the analysis adapting them to your specific questions and
# dataset. When we meet tyo discuss your data, I expect for you to bring some of the basic plots
# described here and a list of questions to discuss together.
## Happy phyloseqing!
