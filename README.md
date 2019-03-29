# How to Run the Pipeline

TODO: destinationDir need to agree in both runSentieonPL.sh and config
TODO: setting scratch to the same scratch directory

All of the fastq files that you want to push through the pipeline need to be in the same directory. If the fastqs are in different
directories you can create symbolic links that are all in the same directory for the fastqs of interest.

You need a config file to run the pipeline. You can create a new one if you'd like, or you can use one of the ones provided in the `./configs` directory. The config file you end up using must be named `nextflow.config` and must be in the same directory as the `sentieon.nf` file. You'll need to rename and move one of these files in order to run the pipeline. You should execute a command similar to this:

`cp ./bin/demux.config ./nextflow.config`

This `nextflow.config` file and the `runSentieonPL.sh` script should be the only two files you have to edit.

The `nextflow.config` file contains all of the parameters specific to your project that nextflow needs to successfully run your
data through to the end. In the `params` section, you'll want to make sure that the following variables are set to work with your
data: `project`, `dataDir`, `isDemuxNeeded`, `barcodeFile`, `sampleKey`, `isIntervalNeeded`, and `bedFile`. In the `process`
section, you should only change the `clusterOptions` variable. The `clusterOptions` variable sets the SLURM parameters for each
job that is scheduled to a CHPC cluster.

The pipeline is designed to run on the Pezzolesi node (notch026). Once set to run, more jobs will be scheduled on whatever
cluster you want. The pipeline is started through the `runSentieonPL.sh` script with this command:

`sbatch runSentieonPL.sh new`

 or 

`sbatch runSentieonPL.sh resume`

You can choose the cluster you want the pipeline to schedule job to by modifying the `SLURM_CLUSTERS` variable in the
`runSentieonPL.sh` script. Note that your choice of cluster needs to be the same in both this variable and the `clusterOptions`
variable in the `nextflow.config` file. You can change the name of the scratch directory you wish to save your project to by
modifying the `scratchDir` variable.

The pipeline is set to run for 5 days. If you think it'll take longer than that then modify the `--time` flag in the
`runSentieonPL.sh` SLURM header. 

The pipeline stalls out after finishing (still can't get it to recognize when it's done), so you'll have to cancel the job
yourself in order to free up the Pezzolesi node for others to use.
