#!/bin/sh -l

#SBATCH --nodes=1                 ##number of requested nodes (a node is a group of CPUs)
#SBATCH -n 12		           ##number of requested CPUs
#SBATCH --time=24:00:00             ##time requested
#SBATCH --job-name soil_pipeline2    ##name of job
#SBATCH -A microbiome                 ##name of the queue you are using. Could be scholar or microbiome if you are on snyder

#########################################################################

echo "Start time"
date +"%d %B %Y %H:%M:%S"

#step 1, cd into the proper directory. This directory must already exist

cd $RCAC_SCRATCH
cd qiime/qiime2-a-tutorial

pwd

#step 2, Load Modules

echo "loading bioinfo"
module load bioinfo

echo "loading qiime2-2022.8"

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/faith_pd_vector.qza \
  --m-metadata-file /scratch/bell/pandinoc/invitro/MetadataS2.txt \
  --o-visualization core-metrics-results/faith-pd-group-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/shannon_vector.qza \
  --m-metadata-file /scratch/bell/pandinoc/invitro/MetadataS2.txt \
  --o-visualization core-metrics-results/shannon-group-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/observed_features_vector.qza \
  --m-metadata-file /scratch/bell/pandinoc/invitro/MetadataS2.txt \
  --o-visualization core-metrics-results/observed-features-group-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/evenness_vector.qza \
  --m-metadata-file /scratch/bell/pandinoc/invitro/MetadataS2.txt \
  --o-visualization core-metrics-results/evenness-group-significance.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/bray_curtis_distance_matrix.qza \
  --m-metadata-file  /scratch/bell/pandinoc/invitro/MetadataS2.txt \
  --m-metadata-column treatment \
  --o-visualization core-metrics-results/bray-curtis-treatment-significance.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file  /scratch/bell/pandinoc/invitro/MetadataS2.txt \
  --m-metadata-column treatment \
  --o-visualization core-metrics-results/weighted-unifrac-subject-treatment-significance.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/jaccard_distance_matrix.qza \
  --m-metadata-file  /scratch/bell/pandinoc/invitro/MetadataS2.txt \
  --m-metadata-column treatment \
  --o-visualization core-metrics-results/jaccard-treatment-significance.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file  /scratch/bell/pandinoc/invitro/MetadataS2.txt \
  --m-metadata-column treatment \
  --o-visualization core-metrics-results/unweighted-unifrac-subject-treatment-significance.qzv \
  --p-pairwise

echo "End time"
date +"%d %B %Y %H:%M:%S"

#to run: sbatch qiime_pipeline_step2.slurm
