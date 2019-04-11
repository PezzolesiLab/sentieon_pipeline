#!/usr/bin/env bash

#SBATCH --time=5-00:00:00
#SBATCH --nodes=1
#SBATCH -o pipelineKickOff-%j.out
#SBATCH -e pipelineKickOff-%j.err
#SBATCH --mail-user=julio.fierro@hsc.utah.edu,scott.frodsham@hsc.utah.edu
#SBATCH --mail-type=END
#SBATCH --account=pezzolesi-np
#SBATCH --partition=pezzolesi-shared-np
##SBATCH --partition=ember
#SBATCH --ntasks=6
#SBATCH --mem=32G

# this vvvv should be either 'resume' or 'new'
resume=$1
# rename this vvvv to create a new directory for your project (i.e. if you want to start over without deleting what you've already done)
scratchDir=$scr/run-node-nf
# rename this vvvv to change clusters (e.g. kingspeak, notchpeak, ember, lonepeak, etc.)
#SLURM_CLUSTERS="notchpeak"
#export SLURM_CLUSTERS

if [[ $resume == "resume" ]]; then
    if [ -d $scratchDir ]; then
        echo -e "\nResuming your previous run\n"
        #cp ./sentieon.nf ./nextflow.config $scratchDir
        #cp -r ./bin $scratchDir
        cd $scratchDir
        nextflow run -with-report -with-trace -with-timeline -with-dag dag.html sentieon.nf -resume
    else
        echo "There's nothing to resume (scratch directory doesn't exist)"
    fi
elif [[ $resume == "new" ]]; then
    if [ ! -d $scratchDir ]; then
        echo -e "\nStarting a new project\n"
        mkdir -p $scratchDir/{results/{fastp,fastqc,bqsr,bam/{stats,coverage},gvcf,vcf/stats},bin}
        cp ./sentieon.nf ./nextflow.config $scratchDir
        cp -r ./bin $scratchDir
        cd $scratchDir
        nextflow run -with-report -with-trace -with-timeline -with-dag dag.html sentieon.nf
    else
        echo "A project already exists. Delete it and start over or resume it"
    fi
else
    echo "takes string argument 'resume' or 'new'"
fi
