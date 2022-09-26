# Capabilities

The pipeline can handle FASTQ, BAM, and GVCF files as inputs. It can handle single-end and paired-end reads. It can handle
targeted, whole exome, and whole genome data. It can demultiplex fastqs.

# How to Run the Pipeline

## FASTQs

All of the fastq files that you want to push through the pipeline need to be in the same directory. If the fastqs are in different directories you can create symbolic links that are all in the same directory for the fastqs of interest.

## Config File

You need a config file to run the pipeline. You can create a new one if you'd like, or you can use one of the ones provided in the `./configs` directory. The config file you end up using must be named `nextflow.config` and must be in the same directory as the `sentieon.nf` file. You'll need to rename and move one of these files in order to run the pipeline. You should execute a command similar to this:

`cp ./bin/demux.config ./nextflow.config`

This `nextflow.config` file and the `runSentieonPL.sh` script should be the only two files you have to edit.

The `nextflow.config` file contains all of the parameters specific to your project that nextflow needs to successfully run your
data through to the end. In the `params` section, you'll want to make sure that the following variables are set to work with your
data: `project`, `dataDir`, `isDemuxNeeded`, `barcodeFile`, `sampleKey`, `isIntervalNeeded`, and `bedFile`. In the `process`
section, you should only change the `clusterOptions` and the `scratch` variables. The `clusterOptions` variable sets the SLURM parameters for each
job that is scheduled to a CHPC cluster. The `scratch` directory specifies the directory where the heavy lifting will be done and where the end results of the pipeline. Be sure that the `scratch` variable matches the `scratchDir` variable in the `runSentieonPL.sh` script.

The pipeline is designed to run on the Pezzolesi node (notch026). Our node has 32 CPUs on it. By using the "shared" feature, we can schedule jobs on each of these CPUs individually or on groups of them. Once it's running, nextflow should run on a small set of CPUs, the subsequent jobs will be scheduled using the remaining CPUs. The command that controls this behavior is `clusterOptions`.
[//]: #(## runSentieonPL.sh File)
[//]: #(You need to make sure that the `scratchDir` variable in this file matches the `scratch` variable you specify in your `nextflow.config` file.)

## Starting the Pipeline

The pipeline is started through the `runSentieonPL.sh` script by issuing one of the two commands below from the command line:

`sbatch runSentieonPL.sh new`

 or 

`sbatch runSentieonPL.sh resume`

# Version

Release 201808.08
