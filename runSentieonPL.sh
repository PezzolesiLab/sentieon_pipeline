#!/usr/bin/env bash

#SBATCH --time=3-00:00:00
#SBATCH --nodes=1
#SBATCH -o pipelineKickOff-%j.out
#SBATCH -e pipelineKickOff-%j.err
#SBATCH --mail-user=scott.frodsham@hsc.utah.edu
#SBATCH --mail-type=END
#SBATCH --account=pezzolesi-np
#SBATCH --partition=pezzolesi-np
##SBATCH --partition=ember

resume=$1
scratchDir=$scr/missingBams
#SLURM_CLUSTERS="kingspeak"
#export SLURM_CLUSTERS

if [[ $resume == "resume" ]]; then
    if [ -d $scratchDir ]; then
        cd $scratchDir
        echo "RESUMING"
        nextflow run -with-report -with-trace -with-timeline -with-dag dag.html sentieon.nf -resume
    else
        echo "There's nothing to resume (scratch directory doesn't exist)"
    fi
elif [[ $resume == "new" ]]; then
    if [ ! -d $scratchDir ]; then
        echo "STARTING FRESH"
        mkdir -p $scratchDir/{results/{fastp,fastqc,bqsr,bam/{stats,coverage},gvcf,vcf/stats},bin}
        cp ./sentieon.nf ./nextflow.config $scratchDir
        cp ./bin/pezzAlign $scratchDir/bin
        cd $scratchDir
        nextflow run -with-report -with-trace -with-timeline -with-dag dag.html sentieon.nf
    else
        echo "A project already exists. Delete it and start over or resume it"
    fi
else
    echo "takes string argument 'resume' or 'new'"
fi
