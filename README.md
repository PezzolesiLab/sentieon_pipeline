# How to Run the Pipeline

## FASTQs

All of the fastq files that you want to push through the pipeline need to be in the same directory. If the fastqs are in different directories you can create symbolic links that are all in the same directory for the fastqs of interest.

## Config File

You need a config file to run the pipeline. You can create a new one if you'd like, or you can use one of the ones provided in the `./configs` directory. The config file you end up using must be named `nextflow.config` and must be in the same directory as the `sentieon.nf` file. You'll need to rename and move one of these files in order to run the pipeline. You should execute a command similar to this:

`cp ./configs/demux.config ./nextflow.config`

This `nextflow.config` file and the `runSentieonPL.sh` script should be the only two files you have to edit.

The `nextflow.config` file contains all of the parameters specific to your project that nextflow needs to successfully run your
data through to the end. In the `params` section, you'll want to make sure that the following variables are set to work with your
data: `project`, `dataDir`, `isDemuxNeeded`, `barcodeFile`, `sampleKey`, `isIntervalNeeded`, and `bedFile`. In the `process`
section, you should only change the `clusterOptions` and the `scratch` variables. The `clusterOptions` variable sets the SLURM parameters for each
job that is scheduled to a CHPC cluster. The `scratch` directory specifies the directory where the heavy lifting will be done and where the end results of the pipeline. Be sure that the `scratch` variable matches the `scratchDir` variable in the `runSentieonPL.sh` script.

The pipeline is designed to run on the Pezzolesi node (notch026). Our node has 32 CPUs on it. By using the "shared" feature, we can schedule jobs on each of these CPUs individually or on groups of them. Once it's running, nextflow should run on a small set of CPUs, the subsequent jobs will be scheduled using the remaining CPUs. The command that controls this behavior is `clusterOptions`.

## runSentieonPL.sh File

You need to make sure that the `scratchDir` variable in this file matches the scratch directory you specify in the config file.

## Starting the Pipeline

The pipeline is started through the `runSentieonPL.sh` script by issuing one of the two commands below from the command line:

`sbatch runSentieonPL.sh new`

 or 

`sbatch runSentieonPL.sh resume`

## A Couple of Final Notes

The pipeline is set to run for 5 days. If you think it'll take longer than that then modify the `--time` flag in the
`runSentieonPL.sh` SLURM header.

The pipeline stalls out after finishing (still can't get it to recognize when it's done), so you'll have to cancel the job
yourself in order to free up the Pezzolesi node for others to use.
