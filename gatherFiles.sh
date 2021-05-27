#!/usr/bin/env bash

#SBATCH --time=0-05:00:00
#SBATCH -o findFiles-%j.out
#SBATCH -e findFiles-%j.err
#SBATCH --mail-user=scott.frodsham@hsc.utah.edu
#SBATCH --mail-type=END
#SBATCH --account=pezzolesi-np
#SBATCH --partition=pezzolesi-np

find $pez -type f -name "*fastq.gz" > ~/fastqDotGz.txt
find $pez -type f -name "*fasta" > ~/fasta.txt
find $pez -type f -name "*fasta.gz" > ~/fastaDotGz.txt
find $pez -type f -name "*fq" > ~/fq.txt
find $pez -type f -name "*fq.gz" > ~/fqDotGg.txt
find $pez -type f -name "*fa" > ~/fa.txt
find $pez -type f -name "*fa.gz" > ~/faDotGz.txt
find $pez -type f -name "*bam" > ~/bam.txt
find $pez -type f -name "*sam" > ~/sam.txt
find $pez -type f -name "*g\.vcf" > ~/gDotVcf.txt
find $pez -type f -name "*gvcf" > ~/gvcf.txt
find $pez -type f -name "*cram" > ~/cram.txt
find $pez -type f -name "*g\.vcf\.gz" > ~/gDotVcfDotGz.txt
find $pez -type f -name "*gvcf\.gz" > ~/gvcfDotGz.txt
