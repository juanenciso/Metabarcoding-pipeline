**Pipeline for 16S metabarcoding assembly**

**QIIME2 activation in Linux through conda-environment**
conda activate qiime2-amplicon-2023.9
**Importing data**                                 
qiime tools import \ 
  --type 'SampleData[PairedEndSequencesWithQuality]' \ 
  --input-path pe-64-manifest.tsv \ 
  --output-path paired-end-demux.qza \ 
  --input-format PairedEndFastqManifestPhred33 

**Visualization of the file paired-end-demux.qza**
qiime demux summarize \ 
  --i-data paired-end-demux.qza \ 
  --o-visualization demux.qzv 

  **DADA2 denoising of high-quality reads**
  qiime dada2 denoise-paired \ 
  --i-demultiplexed-seqs paired-end-demux.qza \ 
  --p-trim-left-f 21 \ 
  --p-trim-left-r 21 \ 
  --p-trunc-len-f 235 \ 
  --p-trunc-len-r 180 \ 
  --o-table table.qza \ 
  --o-representative-sequences rep-seqs.qza \ 
  --o-denoising-stats denoising-stats.qza \ 
  --p-n-reads-learn 50000 \ 
  --p-max-ee-f 2 \ 
  --p-max-ee-r 2 \ 
  --verbose 

**Visualization of the denoised data**
qiime metadata tabulate \ 
  --m-input-file denoising-stats.qza \ 
  --o-visualization denoising-stats.qzv 

**Visualization of the frequency per sample and feature**
qiime feature-table summarize \ 
  --i-table table.qza \ 
  --o-visualization table.qzv 

**Visualization of representative sequences**
qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \ 
  --o-visualization rep-seqs.qzv 

**Clustering features using de novo VSEARCH algorithm at 97% identity**
qiime vsearch cluster-features-de-novo \ 
  --i-table table.qza \ 
  --i-sequences rep-seqs.qza \ 
  --p-perc-identity 0.97 \ 
  --o-clustered-table table-dn-97.qza \
  --o-clustered-sequences rep-seqs-dn-97.qza 

**Download silva-138-99-nb-classifier.qza**
Download Silva 138 99% OTUs from 515F/806R region of sequences from the QIIME2 data resources. https://docs.qiime2.org/2023.9/data-resources/. The name of the corresponding files is silva-138-99-seqs.qza and silva-138-99-tax.qza.

**Taxonomical classification using VSEARCH algorithm**
qiime feature-classifier classify-consensus-vsearch \
  --i-query rep-seqs-dn-97.qza \ 
  --i-reference-reads silva-138-99-seqs.qza \ 
  --i-reference-taxonomy silva-138-99-tax.qza \ 
  --o-classification taxonomy-dn-97.qza \ 
  --o-search-results search-results.qza 

**Visualization of the classified OTUs**
qiime taxa barplot \ 
  --i-table table-dn-97.qza \
  --i-taxonomy taxonomy-dn-97.qza \ 
  --m-metadata-file metadata.tsv \ 
  --o-visualization tax-barplot.qzv

**Export of the BIOM file from table-dn-97.qza**
qiime tools export \ 
  --input-path table-dn-97.qza \ 
  --output-path exported-table 

**Export of the taxonomy.tsv from taxonomy-dn-97.qza into exported-taxonomy folder**
qiime tools export \
  --input-path taxonomy-dn-97.qza \ 
  --output-path exported-taxonomy 

**Export the filtered sequences in. tsv format from denoising-stats.qza**
qiime tools export \ 
  --input-path denoising-stats.qza \ 
  --output-path statistics 

**Export the DNA sequences in a fasta format into sequences-rep-seq folder**
qiime tools export \ .
  --input-path rep-seqs-dn-97.qza \ 
  --output-path sequences-rep-seqs 

**Conversion of the BIOM file into tab-separated**
biom convert -i feature-table.biom -o feature-table.txt --to-tsv


**Data processing in RStudio**

**Install the following packages in R**
install.packages("vegan") 
install.packages("phyloseq") 
install.packages("tidyverse") 
install.packages("patchwork")
install.packages("agricolae") 
install.packages("FSA") 
install.packages("rcompanion")
install.packages("ggplot2")
install.packages("viridis")
**Load the required packages**
library(vegan)
library(phyloseq)
library(tidyverse)
library(patchwork)
library(agricolae)
library(FSA)
library(rcompanion)
library(ggplot2)
library(viridis)

**Import the OTU, metadata and taxonomy tables** 
data_otu <- read.table("feature-table.txt", header = TRUE) 
data_grp <- read.table("metadata.txt", header=TRUE, stringsAsFactors = TRUE) 
data_taxo <- read.table("taxonomy.txt", header = TRUE, fill = TRUE) 

**Create a phyloseq object**
OTU <- otu_table(as.matrix(data_otu), taxa_are_rows = FALSE) 
SAM <- sample_data(data_grp, errorIfNULL = TRUE 
TAX <- tax_table(as.matrix(data_taxo)) 
data_phylo <- phyloseq(OTU, TAX, SAM) 

**Filter OTUs to exclude Archaea, Eukarya and Mitochondria**
data_phylo_filtered <- subset_taxa(
  data_phylo,
  !LCA_simplified %in% c("Archaea", "Eukaryota") & !Family %in% c("Mitochondria")
)
data_otu_filtered <- as.data.frame(otu_table(data_phylo_filtered))
data_taxo_filtered <- as.data.frame(tax_table(data_phylo_filtered))

**Calculation of richness**
data_richness <- estimateR(data_otu_filtered) 

**Calculation of Pielou´s evenness** 
data_evenness <- diversity(data_otu_filtered) / log(specnumber(data_otu_filtered))  

**Calculation of Shannon index**
data_shannon <- diversity(data_otu_filtered, index = "shannon") 

**Combination of indices** 
data_alphadiv <- cbind(data_grp, t(data_richness), data_shannon, data_evenness) 

**Removal of single data frames**
rm(data_richness, data_evenness, data_shannon)

**Transformation of the data in tidy format**
data_alphadiv_tidy <- 
  data_alphadiv %>%
  mutate(sample_id = rownames(data_alphadiv)) %>%
  gather(key = alphadiv_index,
         value = obs_values,
         -locality) 
**Plotting of indices **
P1 <- ggplot(data_alphadiv, aes(x = locality, y = S.obs)) +
  geom_boxplot(fill=c("blue")) +
  labs(title= 'Richness', x = ' ', y = '', tag = "A") +
  geom_point() 
P2 <- ggplot(data_alphadiv, aes(x = locality, y = S.chao1)) +
  geom_boxplot(fill=c("blue")) +
  labs(title= 'Chao1', x = ' ', y = '', tag = "B") +
  geom_point() 
P3 <- ggplot(data_alphadiv, aes(x = locality, y = data_evenness)) +
  geom_boxplot(fill=c("blue")) +
  labs(title= 'Evenness', x = ' ', y = '', tag = "C") +
  geom_point() 
P4 <- ggplot(data_alphadiv, aes(x = locality, y = data_shannon)) +
  geom_boxplot(fill=c("blue")) +
  labs(title= 'Shannon', x = ' ', y = '', tag = "D") +
  geom_point() 
**Plot all the alpha indexes in a single graph**
P1 | P2 | P3 | P4 

**Run a non-parametric test (e.g Kruskal-Wallis test and post-hoc Dunn’s test)**
kruskal.test(data_shannon ~ locality, data = data_alphadiv) 
PT <- dunnTest(data_shannon ~ locality, data = data_alphadiv, method="bh") 
PT2 <- PT$res  
cldList(comparison = PT2$Comparison, p.value = PT2$P.adj, threshold  = 0.05) 

**Run a parametric test for every index (e.g ANOVA)** 
Shapiro.test(data_alphadiv$S.obs) # replace “S.obs” with the variable of interest.
aov_test_shannon <- aov(data_shannon ~ locality, data = data_alphadiv) 
aov_test_evenness <- aov(data_evenness~ locality, data = data_alphadiv) 
aov_test_richness <- aov(S.obs~ locality, data = data_alphadiv) 

**Running the Tukey-Kramer test.** 
hsd_test <- TukeyHSD(aov_test_richness) 
hsd_res <- HSD.test(aov_test_richness, "locality", group=T)$groups 

**Beta diversity analyses**
Import the OTU table, the metadata and taxonomy files, if they are not imported from the previous steps
Load data (Ensure OTU table, metadata and taxonomy data are loaded).

**Data filtering of low number of reads**
data_phylo_filt = filter_taxa(data_phylo, function(x) sum(x > 2) > (0.15 * length(x)), TRUE) 
data_otu_filt = data.frame(otu_table(data_phylo_filt))
dim(data_otu_filt) 

**Set for reproducibility**
set.seed(1892) 

**Rarefaction**
OTU_filt_rar = rarefy_even_depth(otu_table(data_phylo_filt), rngseed = TRUE, replace = FALSE) 

**Construction of a separate file**
data_otu_filt_rar = data.frame(otu_table(OTU_filt_rar)) 

**Creation of a phyloseq object**
data_phylo_filt_rar <- phyloseq(OTU_filt_rar, TAX, SAM) 

**Calculation of Bray-Curtis distances**
dist_bc <- as.matrix(vegdist(data_otu_filt_rar, method = "bray")) 

**Principal Component Analyses (PCoA) Calculation**
pcoa_bc = ordinate(data_phylo_filt_rar, "PCoA", "bray") 

**PCoA visualization**
sample_metadata <- as(sample_data(data_phylo_filt_rar), "data.frame") 
ordination_scores <- as.data.frame(pcoa_bc$vectors) 
ordination_data <- cbind(ordination_scores, sample_metadata) 
ggplot(ordination_data, aes(x = Axis.1, y = Axis.2, color = locality)) +
  geom_point(size = 3) +
  stat_ellipse(aes(group = locality), type = "norm", level = 0.95) +
  labs(
    title = "PCoA with Bray-Curtis Distances",
    color = "Locality",
    x = "PCoA Axis 1",
    y = "PCoA Axis 2"
  ) +
  theme_minimal() 
  
**Principal Component Analyses (PCoA) Calculation**
adonis2(data_otu_filt_rar~locality,data=data_grp, permutations=9999 , method="bray") 

**Gamma diversity**
Load the OTU table (see Note 15)
Tab <- read_tsv("OTU_table.txt",
                col_types = cols(otu_id = col_character(),
                                 .default = col_number())) 

**Matrix organization with pivot longer function**
dat <- Tab %>%
  pivot_longer(-otu_id, names_to = "sample_id", values_to = "count") 
Load the taxonomy table
Tax <- read_tsv("taxonomy.txt", col_types = cols(.default = "character")) 

**Add taxonomy table to the conjoin matrix**
dat <- dat %>%
  left_join(Tax, by = "otu_id") 
Load the metadata table
Meta <- read_tsv("metadata.txt",
                 col_types = cols(.default = col_character())) 

**Add metadata table to the general matrix**
dat <- dat %>%
  left_join(Meta, by = "sample_id") 
  
**Visualization of the gamma diversity in ggplot2**
dat %>%
  ggplot(aes(x = sample_id, y = count)) +
  facet_grid(~ locality, scales = "free_x", space = "free_x") +
  geom_bar(aes(fill =Genus), stat = "identity", position = "fill", width = 1) +
  scale_fill_brewer(palette = "Paired") 

