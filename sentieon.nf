#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

/*
 * Pezzolesi Lab Sentieon Pipeline
 * Author: Scott Frodsham
 * Copyright (c) 2018: none (lol)
*/

// Log info.
log.info "=================================="
log.info "          Pezzolesi Lab           "
log.info "           Sentieon PL            "
log.info "=================================="
log.info "Version: $params.version"
log.info ""
log.info "~~ With the following VQSR knowns ~~"
//log.info "reference: $params. "
//log.info "dbsnp: $params. "
//log.info "hapmap: $params. "
//log.info "1G snps: $params. "
//log.info "1G indels: $params. "
log.info "~~ Set to the following parameters ~~"
//log.info "dbsnp: $params. _par "
//log.info "hapmap: $params. _par "
//log.info "1G snps: $params. _par "
//log.info "1G indels: $params. _par "
log.info "~~ Using the following region bed files ~~"
log.info "$params. beds"
log.info "=================================="

// ---------------------------------------- //

process getFastqs {

    output:
    stdout pooledSampleFqs_ch

    """
    // run script that captures ids relevant to demuxing and filepath (e.g. "flow_cell,index,lane,read,filepath") and outputs to stdout.
    // capture this stout into a channel
    """
}

// do I need to split pooledSamples_ch ?

process fastqThruFastp {
    tag { pooled_fq_id }

    input:
    set val(pooled_fq_id), val(primaryFqLoc) from pooledSampleFqs_ch.splitCsv().map{ id, path -> tuple(id, file(path))}

    output:
    // output of fastp program into channel

    """
    // run each pooled fastq through fastp
    """
}

process demux {
    tag { pooled_fq_id }

    input:
    // output of fastqThruFastp

    output:
    file 'out_dir/*' into ind_fqs

    """
    // run pooled fqs through defq
    """
}

}
