# 🧬 Pipeline for 16S Metabarcoding Assembly

This pipeline performs 16S rRNA metabarcoding data analysis using QIIME2 and RStudio. As input file, I used the FASTQ of the amplicon sequencing products of the V4 hypervariable region of this gene.

---

##  QIIME2

##  Step 1: Activate QIIME2 in conda environment

```bash
conda activate qiime2-amplicon-2023.9
```

##  Step 2: Import Data

```bash
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path pe-64-manifest.tsv \
  --output-path paired-end-demux.qza \
  --input-format PairedEndFastqManifestPhred33
```

##  Step 3: Visualize Imported Data

```bash
qiime demux summarize \
  --i-data paired-end-demux.qza \
  --o-visualization demux.qzv
```

##  Step 4: DADA2 Denoising

```bash
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
```

##  Step 5: Visualizations

```bash
qiime metadata tabulate \
  --m-input-file denoising-stats.qza \
  --o-visualization denoising-stats.qzv

qiime feature-table summarize \
  --i-table table.qza \
  --o-visualization table.qzv

qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --o-visualization rep-seqs.qzv
```

##  Step 6: Cluster Features

```bash
qiime vsearch cluster-features-de-novo \
  --i-table table.qza \
  --i-sequences rep-seqs.qza \
  --p-perc-identity 0.97 \
  --o-clustered-table table-dn-97.qza \
  --o-clustered-sequences rep-seqs-dn-97.qza
```

##  Step 7: Taxonomic Classification

```bash
qiime feature-classifier classify-consensus-vsearch \
  --i-query rep-seqs-dn-97.qza \
  --i-reference-reads silva-138-99-seqs.qza \
  --i-reference-taxonomy silva-138-99-tax.qza \
  --o-classification taxonomy-dn-97.qza \
  --o-search-results search-results.qza
```

##  Step 8: Taxa Barplot

```bash
qiime taxa barplot \
  --i-table table-dn-97.qza \
  --i-taxonomy taxonomy-dn-97.qza \
  --m-metadata-file metadata.tsv \
  --o-visualization tax-barplot.qzv
```

##  Step 9: Export results

```
qiime tools export --input-path table-dn-97.qza --output-path exported-table
qiime tools export --input-path taxonomy-dn-97.qza --output-path exported-taxonomy
biom convert -i feature-table.biom -o feature-table.txt --to-tsv
```

📊 RStudio Analysis

##  Step 10: Install R Packages

```
install.packages(c("vegan", "phyloseq", "tidyverse", "patchwork",
                   "agricolae", "FSA", "rcompanion", "ggplot2", "viridis"))
```

##  Step 11: Load and Prepare Data

```
library(phyloseq)
data_otu <- read.table("feature-table.txt", header = TRUE)
data_grp <- read.table("metadata.txt", header = TRUE, stringsAsFactors = TRUE)
data_taxo <- read.table("taxonomy.txt", header = TRUE, fill = TRUE)

OTU <- otu_table(as.matrix(data_otu), taxa_are_rows = FALSE)
SAM <- sample_data(data_grp)
TAX <- tax_table(as.matrix(data_taxo))
data_phylo <- phyloseq(OTU, TAX, SAM)
```

##  Step 12: Filter Taxa

```
data_phylo_filtered <- subset_taxa(
  data_phylo,
  !LCA_simplified %in% c("Archaea", "Eukaryota") & !Family %in% c("Mitochondria")
)

P1 <- ggplot(data_alphadiv, aes(x = locality, y = S.obs)) +
18
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

P1 | P2 | P3 | P4
```

##  Step 13: Alpha Diversity

```
data_otu_filtered <- as.data.frame(otu_table(data_phylo_filtered))
data_richness <- estimateR(data_otu_filtered)
data_evenness <- diversity(data_otu_filtered) / log(specnumber(data_otu_filtered))
data_shannon <- diversity(data_otu_filtered, index = "shannon")
data_alphadiv <- cbind(data_grp, t(data_richness), data_shannon, data_evenness)

kruskal.test(data_shannon ~ locality, data = data_alphadiv) 
PT <- dunnTest(data_shannon ~ locality, data = data_alphadiv,
method="bh") 
PT2 <- PT$res 
cldList(comparison = PT2$Comparison, p.value = PT2$P.adj, threshold =
0.05)

aov_test_shannon <- aov(data_shannon ~ locality, data = data_alphadiv)
aov_test_evenness <- aov(data_evenness~ locality, data = data_alphadiv) 
aov_test_richness <- aov(S.obs~ locality, data = data_alphadiv) 

hsd_test <- TukeyHSD(aov_test_richness) 
hsd_res <- HSD.test(aov_test_richness, "locality", group=T)$groups
```

##  Step 14: Beta Diversity

```
data_phylo_filt <- filter_taxa(data_phylo, function(x) sum(x > 2) > (0.15 * length(x)), TRUE)
set.seed(1892)
OTU_filt_rar <- rarefy_even_depth(otu_table(data_phylo_filt), rngseed = TRUE, replace = FALSE)
data_phylo_filt_rar <- phyloseq(OTU_filt_rar, TAX, SAM)

dist_bc <- vegdist(as.data.frame(otu_table(data_phylo_filt_rar)), method = "bray")
pcoa_bc <- ordinate(data_phylo_filt_rar, "PCoA", "bray")
ordination_scores <- as.data.frame(pcoa_bc$vectors)
ordination_data <- cbind(ordination_scores, as(sample_data(data_phylo_filt_rar), "data.frame"))

ggplot(ordination_data, aes(x = Axis.1, y = Axis.2, color = locality)) +
  geom_point(size = 3) +
  stat_ellipse(aes(group = locality), type = "norm", level = 0.95) +
  labs(title = "PCoA with Bray-Curtis Distances", x = "PCoA Axis 1", y = "PCoA Axis 2") +
  theme_minimal()
```

##  Step 15: PERMANOVA

```
adonis2(as.data.frame(otu_table(data_phylo_filt_rar)) ~ locality, data = data_grp, permutations = 9999, method = "bray")
```

##  Step 16: Gamma Diversity

```
Tab <- read_tsv("OTU_table.txt", col_types = cols(otu_id = col_character(), .default = col_number()))
dat <- Tab %>% pivot_longer(-otu_id, names_to = "sample_id", values_to = "count")
Tax <- read_tsv("taxonomy.txt", col_types = cols(.default = "character"))
dat <- dat %>% left_join(Tax, by = "otu_id")
Meta <- read_tsv("metadata.txt", col_types = cols(.default = col_character()))
dat <- dat %>% left_join(Meta, by = "sample_id")

ggplot(dat, aes(x = sample_id, y = count)) +
  facet_grid(~ locality, scales = "free_x", space = "free_x") +
  geom_bar(aes(fill = Genus), stat = "identity", position = "fill", width = 1) +
  scale_fill_brewer(palette = "Paired")
```



