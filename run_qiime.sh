#!/usr/bin/env bash
set -euo pipefail

# Activate QIIME 2
source ~/miniconda3/etc/profile.d/conda.sh
conda activate qiime2-2024.10

# Project path (opcional)
BASE_DIR="$PWD"

#Activate conda#

conda activate qiime2-2022.8

#Import data using the manifest file#

qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path se-33-manifest.tsv \
  --output-path single-end-demux.qza \
  --input-format SingleEndFastqManifestPhred33V2

#Visualization of the file sigle-end-demux.gza#

qiime demux summarize \
--i-data single-end-demux.qza \
--o-visualization demux.qzv

qiime dada2 denoise-single \
  --p-trim-left 0 \
  --p-trunc-len 230 \
  --i-demultiplexed-seqs single-end-demux.qza \
  --o-representative-sequences rep-seq.qza \
  --o-table table.qza \
  --o-denoising-stats stats.qza

qiime metadata tabulate \
  --m-input-file stats.qza \
  --o-visualization stats.qzv

qiime feature-table summarize \
  --i-table table.qza \
  --o-visualization table.qzv \
  --m-sample-metadata-file Juan-meta-data-file.tsv


qiime feature-table tabulate-seqs \
  --i-data rep-seq.qza \
  --o-visualization  rep-seq.qzv

qiime feature-classifier classify-sklearn \
  --i-classifier silva-138-99-nb-classifier.qza \
  --i-reads rep-seq.qza \
  --o-classification taxonomy.qza

qiime taxa barplot \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file Juan-meta-data-file.tsv \
  --o-visualization tax-barplot.qzv

#Phylogenetic tree for the diversity analyses#
qiime phylogeny align-to-tree-mafft-fasttree --i-sequences rep-seq.qza --o-alignment aligned-rep-seq.qza --o-masked-alignment masked-aligned-rep-seq.qza --o-tree unrooted-tree.qza --o-rooted-tree rooted-tree.qza

#Alpha rarefaction curve#
qiime diversity alpha-rarefaction --i-table table.qza --i-phylogeny rooted-tree.qza --p-max-depth 1172 --m-metadata-file Juan-meta-data-file.tsv --o-visualization alpha-rarefaction.qzv
#Alpha and beta diversity analysis#
qiime diversity core-metrics-phylogenetic --i-phylogeny rooted-tree.qza --i-table table.qza --p-sampling-depth 1103 --m-metadata-file Juan-meta-data-file.tsv --output-dir metrics-results
#Alpha group significance analysis#
qiime diversity alpha-group-significance --i-alpha-diversity core-metrics-results/faith_pd_vector.qza --m-metadata-file Juan-meta-data-file.tsv --o-visualization core-metrics-results/observed_features_significance.qzv

qiime diversity alpha-group-significance --i-alpha-diversity core-metrics-results/observed_features_vector.qza --m-metadata-file Juan-meta-data-file.tsv --o-visualization core-metrics-results/observed_features_significance.qza
#Beta group significance#
qiime diversity beta-group-significance --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza --m-metadata-file Juan-meta-data-file.tsv --m-metadata-column species --o-visualization core-metrics-results/unweighted-unifrac-body-site-significance.qzv --p-pairwise 

qiime diversity beta-group-significance --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza --m-metadata-file Juan-meta-data-file.tsv --m-metadata-column species --o-visualization core-metrics-results/unweighted-unifrac-species-significance.qzv --p-pairwise

qiime diversity beta-group-significance --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza --m-metadata-file Juan-meta-data-file.tsv --m-metadata-column locality --o-visualization core-metrics-results/unweighted-unifrac-locality-significance.qzv --p-pairwise
qiime diversity beta-group-significance --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza --m-metadata-file Juan-meta-data-file.tsv --m-metadata-column locality --o-visualization core-metrics-results/unweighted-unifrac-locality-significance.qzv --p-pairwise




