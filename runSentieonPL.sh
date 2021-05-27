#!/usr/bin/env bash

#SBATCH --time=3-00:00:00
#SBATCH --nodes=1
#SBATCH -o pipelineKickOff-%j.out
#SBATCH -e pipelineKickOff-%j.err
#SBATCH --mail-user=scott.frodsham@hsc.utah.edu
#SBATCH --mail-type=BEGIN,END
#SBATCH --account=pezzolesi-np
#SBATCH --partition=pezzolesi-np
#SBATCH --ntasks=6
#SBATCH --mem=32G

# this vvvv should be either 'resume' or 'new'
resume=$1
configFile=$2

## finds the scratch file around line 80 of the config file
#configScratch=$(awk '/^    scratch/ {print}' nextflow.config)
## split it on the '=' and stick it in the scratchDir variable
#scratchVar3=${configScratch##*= }
## remove the quotes from the path
#scratchVar2=$(sed -e 's/"//' -e 's/"$//' <<< "$scratchVar3")
## replace $USER with the actual username
#scratchDir=$(sed "s/\$USER/$USER/" <<< "$scratchVar2")

# OR you can name it manually (just be sure that this is the same scratch dir given in the config file)
scratchDir="$scr/sentieon_HAVCR1_exon4_variantCall"

if [ -z $2 ]; then
    if [[ $resume == "resume" ]]; then
        if [ -d $scratchDir ]; then
            echo -e "\nResuming your previous run\n"
            cp ./sentieon.nf ./nextflow.config $scratchDir
            #cp -r ./bin $scratchDir
            cd $scratchDir
            nextflow run -with-report -with-trace -with-timeline -with-dag dag.html sentieon.nf -resume
        else
            echo "There's nothing to resume (scratch directory doesn't exist)"
        fi
    elif [[ $resume == "new" ]]; then
        if [ ! -d $scratchDir ]; then
            echo -e "\nStarting a new project\n"
            mkdir -p $scratchDir/{results/{fastp,fastqc,bqsr,vqsr,bam/{stats,coverage},gvcf,vcf/stats},bin}
            cp ./sentieon.nf ./nextflow.config $scratchDir
            cp -r ./bin $scratchDir
            cd $scratchDir
            nextflow run -with-report -with-trace -with-timeline -with-dag dag.html sentieon.nf
        else
            echo "The directory $scratchDir already exists. Delete it and start over or resume it"
        fi
    else
        echo "Takes string argument 'resume' or 'new'"
    fi
else
    if [[ $resume == "resume" ]]; then
        if [ -d $scratchDir ]; then
            echo -e "\nResuming your previous run\n"
            cp ./sentieon.nf ./nextflow.config $scratchDir
            cp -r ./bin $scratchDir
            cd $scratchDir
            nextflow -C $configFile run -with-report -with-trace -with-timeline -with-dag dag.html sentieon.nf -resume
        else
            echo "There's nothing to resume (scratch directory doesn't exist)"
        fi
    elif [[ $resume == "new" ]]; then
        if [ ! -d $scratchDir ]; then
            echo -e "\nStarting a new project\n"
            mkdir -p $scratchDir/{results/{fastp,fastqc,bqsr,vqsr,bam/{stats,coverage},gvcf,vcf/stats},bin}
            cp ./sentieon.nf ./nextflow.config $scratchDir
            cp -r ./bin $scratchDir
            cd $scratchDir
            nextflow -C $configFile  run -with-report -with-trace -with-timeline -with-dag dag.html sentieon.nf
        else
            echo "A project already exists. Delete it and start over or resume it"
        fi
    else
        echo "Takes string argument 'resume' or 'new'"
    fi
fi
