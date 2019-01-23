#!/usr/bin/env bash

#SBATCH --time=3-12:00:00
#SBATCH --nodes=1
#SBATCH -o annotateVEP-%j.out
#SBATCH -e annotateVEP-%j.err
#SBATCH --mail-user=scott.frodsham@hsc.utah.edu
#SBATCH --mail-type=END
#SBATCH --account=pezzolesi-np
#SBATCH --partition=pezzolesi-np

module load vep

export RANGE=$(cat $bedfile)
export vcf="$pez/WES_data/RFD_families"
export vep_cache="/scratch/ucgd/serial/common/apps/notchpeak.peaks/vep/release-92/vep_cache"

tabix -h $vcf $RANGE | \
vep --cache --dir_cache $vep_cache --format $vcf --fork 32 --output_file "UCCDS_RFD_vep.vcf.gz" --vcf --offline --assembly "GRCh37" --everything

