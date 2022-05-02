#!/bin/sh -l

#SBATCH --nodes=1                 ##number of requested nodes (a node is a group of CPUs)
#SBATCH -n 12		           ##number of requested CPUs
#SBATCH --time=04:00:00             ##time requested
#SBATCH --job-name qiime_pipeline    ##name of job
#SBATCH -A microbiome                 ##name of the queue you are using. Could be scholar or microbiome if you are on snyder

#########################################################################

echo "Start time"
date +"%d %B %Y %H:%M:%S"

#step 1, cd into the proper directory. This directory must already exist

cd $RCAC_SCRATCH
cd qiime/manifest_import

pwd

#step 2, Load Modules

echo "loading bioinfo"
module load bioinfo

echo "loading qiime2-2021.11"
module load Qiime/2-2021.11

#step 3, Run the code you want

#echo "initiating qiime2 moving pictures tutorial pipeline"

echo "Data Import"
#Make the sequences into a qiime recognizable format

#qiime tools import \
 # --type EMPSingleEndSequences \
  #--input-path emp-single-end-sequences \
  #--output-path emp-single-end-sequences.qza

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path manifest_file.tsv \
  --output-path demux.qza \
  --input-format PairedEndFastqManifestPhred33V2

qiime demux summarize \
  --i-data demux.qza \
  --o-visualization demux.qzv

##De-noise. Error corrections.

echo "dada2"

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs demux.qza \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-trunc-len-f 220 \
  --p-trunc-len-r 220 \
  --o-representative-sequences rep-seqs-full.qza \
  --o-table table-full.qza \
  --o-denoising-stats stats-dada2-full.qza

#Convert .qza to .qzv format

qiime metadata tabulate \
  --m-input-file stats-dada2-full.qza \
  --o-visualization stats-dada2-full.qzv


#Create a .qzv from the output of denoise

qiime feature-table summarize \
  --i-table table-full.qza \
  --o-visualization table-full.qzv \
  --m-sample-metadata-file MetadataEM.tsv

qiime feature-table tabulate-seqs \
  --i-data rep-seqs-full.qza \
  --o-visualization rep-seqs-full.qzv

#Phylogenetic tree

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs-full.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza

#Alpha and Beta diversity analysis

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree.qza \
  --i-table table-full.qza \
  --p-sampling-depth 31560 \
  --m-metadata-file MetadataEM.tsv \
  --output-dir core-metrics-results

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/faith_pd_vector.qza \
  --m-metadata-file MetadataEM.tsv \
  --o-visualization core-metrics-results/faith-pd-group-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/evenness_vector.qza \
  --m-metadata-file MetadataEM.tsv \
  --o-visualization core-metrics-results/evenness-group-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/shannon_vector.qza \
  --m-metadata-file MetadataEM.tsv \
  --o-visualization core-metrics-results/shannon-group-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/observed_features_vector.qza \
  --m-metadata-file MetadataEM.tsv \
  --o-visualization core-metrics-results/observed_features-group-significance.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file MetadataEM.tsv \
  --m-metadata-column sample_type \
  --o-visualization core-metrics-results/unweighted-unifrac-transect-significance.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file MetadataEM.tsv \
  --m-metadata-column sample_type \
  --o-visualization core-metrics-results/weighted-unifrac-transect-group-significance.qzv \
  --p-pairwise

#Alpha Rarefation Plotting

qiime diversity alpha-rarefaction \
  --i-table table-full.qza \
  --i-phylogeny rooted-tree.qza \
  --p-max-depth 4000 \
  --m-metadata-file MetadataEM.tsv \
  --o-visualization alpha-rarefaction.qzv

#Downloading Classifier
#wget \
 # -O "gg-13-8-99-515-806-nb-classifier.qza" \
 # "https://data.qiime2.org/2021.11/common/gg-13-8-99-515-806-nb-classifier.qza"

#Taxonomy
qiime feature-classifier classify-sklearn \
  --i-classifier gg-13-8-99-515-806-nb-classifier.qza \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy.qza

qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv

#Taxa Bar plots
qiime taxa barplot \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file sample-metadata.tsv \
  --o-visualization taxa-bar-plots.qzv

#ANCOM
qiime feature-table filter-samples \
  --i-table table.qza \
  --m-metadata-file sample-metadata.tsv \
  --p-where "[body-site]='gut'" \
  --o-filtered-table gut-table.qza

qiime composition add-pseudocount \
  --i-table gut-table.qza \
  --o-composition-table comp-gut-table.qza

qiime composition ancom \
  --i-table comp-gut-table.qza \
  --m-metadata-file sample-metadata.tsv \
  --m-metadata-column subject \
  --o-visualization ancom-subject.qzv

qiime taxa collapse \
  --i-table gut-table.qza \
  --i-taxonomy taxonomy.qza \
  --p-level 6 \
  --o-collapsed-table gut-table-l6.qza

qiime composition add-pseudocount \
  --i-table gut-table-l6.qza \
  --o-composition-table comp-gut-table-l6.qza

qiime composition ancom \
  --i-table comp-gut-table-l6.qza \
  --m-metadata-file sample-metadata.tsv \
  --m-metadata-column subject \
  --o-visualization l6-ancom-subject.qzv

echo "End time"
date +"%d %B %Y %H:%M:%S"

#to run: sbatch qiime_pipeline.slurm

