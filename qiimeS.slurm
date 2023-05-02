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

qiime tools import \
--type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path /scratch/bell/pandinoc/invitro/Manifest_file2.txt \
--input-format PairedEndFastqManifestPhred33V2 \
--output-path demuxS.qza

qiime demux summarize \
  --i-data demuxS.qza \
  --o-visualization demuxS.qzv

qiime dada2 denoise-paired \
--i-demultiplexed-seqs demux.qza \
--p-trim-left-f 13 \
--p-trim-left-r 13 \
--p-trunc-len-f 278 \
--p-trunc-len-r 220 \
--o-table tableS.qza \
--o-representative-sequences rep-seqsS.qza \
--o-denoising-stats stats-dada2S.qza

qiime metadata tabulate \
  --m-input-file stats-dada2S.qza \
  --o-visualization stats-dada2S.qzv


qiime feature-table summarize \
  --i-table tableS.qza \
  --o-visualization tableS.qzv \
  --m-sample-metadata-file /scratch/bell/pandinoc/invitro/MetadataS2.txt

qiime feature-table tabulate-seqs \
  --i-data rep-seqsS.qza \
  --o-visualization rep-seqsS.qzv


echo "End time"
date +"%d %B %Y %H:%M:%S"

#to run: sbatch qiime_pipeline_step2.slurm
