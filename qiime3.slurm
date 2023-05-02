#!/bin/sh -l

#SBATCH --nodes=1                 ##number of requested nodes (a node is a group of CPUs)
#SBATCH -n 12		           ##number of requested CPUs
#SBATCH --time=24:00:00             ##time requested
#SBATCH --job-name soil_qiime_pipeline    ##name of job
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
module load Qiime/2-2022.8

#step 3, Run the code you want

echo "initiating qiime2 pipeline"


echo "Classify and get .qzv"
qiime feature-classifier classify-sklearn \
  --i-classifier classifier.qza \
  --i-reads rep-seqsS.qza \
  --o-classification taxonomyS.qza

qiime metadata tabulate \
  --m-input-file taxonomyS.qza \
  --o-visualization taxonomyS.qzv


echo "Make a stacked barplot"
qiime taxa barplot \
  --i-table tableS.qza \
  --i-taxonomy taxonomyS.qza \
  --m-metadata-file /scratch/bell/pandinoc/invitro/MetadataS2.txt \
  --o-visualization taxa-bar-plotsS.qzv



echo "End time"
date +"%d %B %Y %H:%M:%S"

#to run: sbatch qiime_pipeline.slurm
