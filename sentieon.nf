#!/usr/bin/nextflow

// Log info
log.info "=================================="
log.info "          Pezzolesi Lab           "
log.info "           Sentieon PL            "
log.info "=================================="
log.info ""
log.info "~~ With the following VQSR knowns ~~"
log.info "reference: $params.reference "
log.info "dbsnp: $params.dbsnp "
log.info "hapmap: $params.hapmap"
log.info "1G snps: $params.snp_1G"
log.info "1G indels: $params.indel_1G "
log.info "mills indels: $params.indel_mills"
log.info "~~ Set to the following parameters ~~"
log.info "hapmap: $params.hapmap_par "
log.info "1G snps: $params.snp_1G_par "
log.info "1G indels: $params.indel_1G_par "
log.info "mills indels: $params.indel_mills_par"
log.info "~~ Using the following region bed files ~~"
log.info "Wuxi: $params.bedFile"
log.info "~~ Software dependencies ~~"
log.info "ANNOVAR"
log.info "fastp/0.19.6"
log.info "multiqc/1.7"
log.info "~~ Modules Needed ~~"
log.info "bcftools/1.7"
log.info "bgzip/1.7"
log.info "samtools/1.7"
log.info "tabix/1.7"
log.info "sentieon/201711.05"
log.info ""
log.info "=================================="

// Get sample name, samples are split into forward and reverse reads

if (params.isDmuxNeeded) {
    
    // read in
    // demultiplex - both
    // trim - both
    //
    // how to delete irrelevant fastqs - parse barcode file
    // how to rename the fastqs. With python script or not?

  libraryReadLaneFq = Channel
    .fromPath( "${params.dataDir}" )
    .map { file ->
        fNameExt = file.baseName
        fName = fNameExt.tokenize('.')[0]
        library = fName.tokenize('_')[0]
        lane = fName.tokenize('_')[6]
        read_num = fName.tokenize('_')[7]
        [ library, read_num, lane, file ]
    }

  libraryReadLaneFq
      .groupTuple()
      .set { demux_in }

  process demuxFq {
      tag { library + "_" + lane }

      input:
      set val(library), val(read_nums), val(lanes), file(fq_files) from demux_in

      output:
      set val(library), val(read1), val(lane), file("${library}_${read1}_${lane}.*.fastq.gz") into demux_out1
      set val(library), val(read2), val(lane), file("${library}_${read2}_${lane}.*.fastq.gz") into demux_out2

      shell:

      fq1 = fq_files[0]
      fq2 = fq_files[1]
      read1 = read_nums[0]
      read2 = read_nums[1]
      lane = lanes[0]

      '''
      fastq-multx \\
      -b \\
      -B !{params.barcodeFile} \\
      !{fq1} \\
      !{fq2} \\
      -o !{library}_!{read1}_!{lane}.%.fastq.gz \\
      -o !{library}_!{read2}_!{lane}.%.fastq.gz
      '''
  }

  demux_out1.join(demux_out2).println()

  // START HERE 3/5/19: trim both w/ cutadapt
}


// specific to UDDCS samples sequenced by WuXi
idReadFq = Channel
    .fromPath( "${params.dataDir}" )
    .map { file ->
        fName = file.baseName
        fName2 = fName.replaceAll("LU01-", "LU01_")
        id = fName2.tokenize('.')[0].tokenize('_')[1]
        read_num = fName2.tokenize('.')[0].tokenize('_')[3]
        [ id, read_num, file ]
    }

idReadFq
    .groupTuple()
    .set { fastp_in }

process runFastp {
    tag { sample_id }

    publishDir "${params.fastp}", mode: 'copy', pattern: '*.html'

    input:
    set val(sample_id), val(read_nums), file(fq_files) from fastp_in

    output:
    set val(sample_id), val(read1), file("${sample_id}_${read1}.trimmed.fastq.gz"), val(sample_id), val(read2), file("${sample_id}_${read2}.trimmed.fastq.gz") into fastp_out
    file("${sample_id}.fastp.report.html")

    shell:

    fq1 = fq_files[0]
    fq2 = fq_files[1]
    read1 = read_nums[0]
    read2 = read_nums[1]

    '''
    fastp \\
    --thread !{params.kp_cpus} \\
    --in1 !{fq1} \\
    --in2 !{fq2} \\
    --out1 "!{sample_id}_!{read1}.trimmed.fastq.gz" \\
    --out2 "!{sample_id}_!{read2}.trimmed.fastq.gz" \\
    --length_required 25 \\
    --low_complexity_filter 5 \\
    --detect_adapter_for_pe \\
    --html "!{sample_id}.fastp.report.html"
    '''
}

fastp_out
    .flatten()
    .collate( 3 )
    .set { collated }

collated.into { fastqc_in; gunzip_fq_in }

process runFastqc {
    tag { sample_id + "_" + read_num }

    input:
    set val(sample_id), val(read_num), file(fq_file) from fastqc_in

    //output:
    //val 'complete' into fastqc_done
    
    shell:
    '''
    fastqc !{fq_file} -o !{params.fastqc} -t !{params.kp_cpus} 
    '''
}

process unzipFastqs {
    tag { sample_id + "_" + read_num }

    input:
    set val(sample_id), val(read_num), file(fq_file) from gunzip_fq_in 

    output:
    set val(sample_id), val(read_num), file("${sample_id}_${read_num}.fastq") into gunzip_fq_out

    """
    gunzip -c ${fq_file} > ${sample_id}_${read_num}.fastq
    """
}

gunzip_fq_out
    .groupTuple()
    .set { bwa_in }
    
process BWA {
    tag { sample_id }

    input:
    set val(sample_id), file(read_nums), file(fq_files) from bwa_in

    output:
    set sample_id, file("${sample_id}.sorted.bam"), file("${sample_id}.sorted.bam.bai") into bwa_out

    shell:

    fq1 = fq_files[0]
    fq2 = fq_files[1]

    '''
    export RG=$(pezzAlign !{fq1})

    ( bwa mem -M \\
    -R $RG \\
    -K 10000000 \\
    -t $SLURM_CPUS_ON_NODE \\
    !{params.reference} \\
    !{fq1} \\
    !{fq2} \\
    || echo -n 'error' ) \\
    | sentieon util sort \\
    -o "!{sample_id}.sorted.bam" \\
    -t $SLURM_CPUS_ON_NODE \\
    --sam2bam -i -
    '''

}

process dedup {
    tag { sample_id }

    input:
    set sample_id, file(bam), file(index) from bwa_out

    output:
    set sample_id, file("${sample_id}.dedup.bam"), file("${sample_id}.dedup.bam.bai") into dedup_out

    shell:
    '''
    sentieon driver \\
    -t $SLURM_CPUS_ON_NODE \\
    -i !{bam} \\
    --algo LocusCollector \\
    --fun score_info "!{sample_id}.score.txt" \\

    sentieon driver \\
    -t $SLURM_CPUS_ON_NODE \\
    -i !{bam} \\
    --algo Dedup \\
    --score_info "!{sample_id}.score.txt" \\
    --metrics "!{sample_id}.metric.txt" \\
    "!{sample_id}.dedup.bam"
    '''
}

process indelRealigner {
    tag { sample_id }

    input:
    set sample_id, file(deduped), file(index) from dedup_out

    output:
    set sample_id, file("${deduped.baseName}.realigned.bam"), file("${deduped.baseName}.realigned.bam.bai") into realign_out

    shell:
    '''
    sentieon driver \\
    -t $SLURM_CPUS_ON_NODE \\
    -r !{params.reference} \\
    -i !{deduped} \\
    --algo Realigner \\
    -k !{params.indel_1G} \\
    -k !{params.indel_mills} \\
    "!{deduped.baseName}.realigned.bam"
    '''
}

process BQSR {
    tag { sample_id }

    publishDir "${params.bam}", mode: 'copy', pattern: '*.{bam,bai}'

    input:
    set sample_id, file(realigned), file(index) from realign_out

    output:
    set sample_id, file("${realigned.baseName}.realigned.bqsr.bam"), file("${realigned.baseName}.realigned.bqsr.bam.bai"), file("${realigned.baseName}.recal_data_table") into bqsr_out
    set sample_id, file("${realigned.baseName}.recal_data_table"), file("${realigned.baseName}.recal_data.table.post") into bqsr_graph_out

    shell:
    '''
    sentieon driver \\
    -t $SLURM_CPUS_ON_NODE \\
    -r !{params.reference} \\
    -i !{realigned} \\
    --algo QualCal \\
    -k !{params.dbsnp} \\
    "!{realigned.baseName}.recal_data_table"

    sentieon driver \\
    -t $SLURM_CPUS_ON_NODE \\
    -r !{params.reference} \\
    -i !{realigned} \\
    -q "!{realigned.baseName}.recal_data_table" \\
    --algo QualCal \\
    -k !{params.dbsnp} \\
    "!{realigned.baseName}.recal_data.table.post" \\
    --algo ReadWriter "!{realigned.baseName}.realigned.bqsr.bam"
    '''
}

//bqsr_out.into { stats_in; flagstat_in; coverageMetrics_in; haplotyper_in }
bqsr_out.into { stats_in; flagstat_in; coverageMetrics_in }

process graphBQSR {
    tag { sample_id }

    publishDir "${params.bqsr}", mode: 'copy'

    input:
    set sample_id, file(table), file(post) from bqsr_graph_out

    output:
    file("${sample_id}.bqsr.pdf")
    //val 'complete' into graph_done

    shell:
    '''
    sentieon driver \\
    -t $SLURM_CPUS_ON_NODE \\
    --algo QualCal \\
    --plot \\
    --before !{table} \\
    --after !{post} \\
    "!{sample_id}.recal.result.csv"

    sentieon plot bqsr \\
    -o "!{sample_id}.bqsr.pdf" \\
    "!{sample_id}.recal.result.csv"
    '''
}

process samStats {
    tag { sample_id }

    publishDir "${params.stats}", mode: 'copy'

    input:
    set sample_id, file(bam), file(index), file(recal_table) from stats_in

    output:
    file("${sample_id}.stats")
    //val 'complete' into samStats_done

    """
    samtools stats ${bam} > "${sample_id}.stats"
    """
}

process samFlagstat {
    tag { sample_id }

    publishDir "${params.stats}", mode: 'copy'

    input:
    set sample_id, file(bam), file(index), file(recal_table) from flagstat_in

    output:
    file("${sample_id}.flagstat")
    //val 'complete' into samFlagstat_done

    """
    samtools flagstat ${bam} > "${sample_id}.flagstat"
    """
}

process coverageMetrics {
    tag { sample_id }

    publishDir "${params.coverage}", mode: 'copy'

    input:
    set sample_id, file(bam), file(index), file(recal_table) from coverageMetrics_in

    output:
    file("${sample_id}.sample_summary")
    //val 'complete' into coverage_done
    
    shell:
    '''
    sentieon driver \\
    -t $SLURM_CPUS_ON_NODE \\
    -i !{bam} \\
    --interval !{params.bedFile} \\
    -r !{params.reference} \\
    --algo CoverageMetrics \\
    --partition sample \\
    --omit_base_output \\
    --omit_interval_stat \\
    --omit_locus_stat \\
    "!{sample_id}"
    '''
}

/*
process haplotyper {
    tag { sample_id }

    publishDir "${params.gvcf}", mode: 'copy'

    input:
    set sample_id, file(recalbam), file(index), file(recalTable) from haplotyper_in

    output:
    set file("${sample_id}.g.vcf.gz"), file("${sample_id}.g.vcf.gz.tbi") into haplotyper_out

    shell:
    '''
    sentieon driver \\
    -t $SLURM_CPUS_ON_NODE \\
    -r !{params.reference} \\
    -i !{recalbam} \\
    -q !{recalTable} \\
    --algo Haplotyper \\
    --emit_mode gvcf \\
    -d !{params.dbsnp} \\
    "!{sample_id}.g.vcf"

    bgzip -c "!{sample_id}.g.vcf" > "!{sample_id}.g.vcf.gz"
    tabix -p vcf "!{sample_id}.g.vcf.gz"
    '''
} 

haplotyper_out
    .toSortedList()
    .transpose()
    .first()
    .set { gvcfTyper_in }

process gvcfTyper {
    tag { "$params.project" }

    input:
    val gvcfs from gvcfTyper_in

    output:
    file ("${params.project}.g.vcf.gz") into gvcfTyper_out

    shell:

    inputGVCFs = gvcfs.join(' -v ')

    '''
    sentieon driver \\
    -t !{params.np_cpus} \\
    -r !{params.reference} \\
    --interval !{params.bedFile} \\
    --algo GVCFtyper \\
    "!{params.project}.g.vcf.gz" \\
    -v !{inputGVCFs}
    '''
}

gvcfTyper_out
    .unique()
    .toList()
    .set { mergeGVCFs_in }

process mergeGVCFs {
    tag { "$params.project" }

    input:
    val chrFiles from mergeGVCFs_in

    output:
    set file("${params.project}_merged.vcf.gz"), file("${params.project}_merged.vcf.gz.tbi") into mergeGVCFs_out

    shell:

    inputVCFs = chrFiles.join(' ')

    '''
    bcftools concat --thread !{params.np_cpus} -O z !{inputVCFs} -o "!{params.project}_merged.vcf.gz"
    tabix -p vcf "!{params.project}_merged.vcf.gz"
    '''
}

process varCalSnp {
    tag { "$params.project" }

    input:
    set file(merged_vcf), file(index) from mergeGVCFs_out

    output:
    set file(merged_vcf), file(index), file("${params.project}_recal.tranches.snp"), file("${params.project}_recal.snp.vcf.gz"), file("${params.project}_recal.snp.vcf.gz.tbi") into varCalSNP_out

    """
    sentieon driver \\
    --thread_count $params.np_cpus \\
    -r $params.reference \\
    --algo VarCal \\
    --vcf $merged_vcf \\
    --var_type SNP \\
    --tranches_file "${params.project}_recal.tranches.snp" \\
    --resource $params.hapmap \\
    --resource_param $params.hapmap_par \\
    --resource $params.omni \\
    --resource_param $params.omni_par \\
    --resource $params.snp_1G \\
    --resource_param $params.snp_1G_par \\
    --annotation DP \\
    --annotation QD \\
    --annotation MQRankSum \\
    --annotation ReadPosRankSum \\
    --annotation FS \\
    "${params.project}_recal.snp.vcf.gz"
    """
}

process applyVarCalSnp {
    tag {"$params.project" }

    input:
    set file(merged_vcf), file(merged_index), file(snpTranchFile), file(recalVCF), file(index) from varCalSNP_out

    output:
    set file("${params.project}_applyRecal.snp.vcf.gz"), file("${params.project}_applyRecal.snp.vcf.gz.tbi") into appliedSNP_out

    """
    sentieon driver \\
    --thread_count $params.np_cpus \\
    -r $params.reference \\
    --algo ApplyVarCal \\
    --vcf $merged_vcf \\
    --var_type SNP \\
    --recal $recalVCF \\
    --tranches_file $snpTranchFile \\
    --sensitivity 90 \\
    "${params.project}_applyRecal.snp.vcf.gz"
    """
}

process varCalIndel {
    tag { "$params.project" }

    input:
    set file(snp_recal_file), file(index) from appliedSNP_out

    output:
    set file(snp_recal_file), file(index), file("${params.project}_recal.tranches.indel"), file("${params.project}_recal.indel.vcf.gz"), file("${params.project}_recal.indel.vcf.gz.tbi") into varCalIndel_out

    """
    sentieon driver \\
    --thread_count $params.np_cpus \\
    -r $params.reference \\
    --algo VarCal \\
    --vcf $snp_recal_file \\
    --var_type INDEL \\
    --tranches_file "${params.project}_recal.tranches.indel" \\
    --resource $params.indel_mills \\
    --resource_param $params.indel_mills_par \\
    --resource $params.indel_1G \\
    --resource_param $params.indel_1G_par \\
    --annotation DP \\
    --annotation MQRankSum \\
    --annotation ReadPosRankSum \\
    --annotation FS \\
    "${params.project}_recal.indel.vcf.gz"
    """
}

process applyVarCalIndel {
    tag { "$params.project" }

    publishDir "${params.vcf}", mode: 'copy'

    input:
    set file(snp_recal_file), file(recal_index), file(IndelTranch), file(recalVCF), file(index) from varCalIndel_out

    output:
    set file("${params.project}_complete.vcf.gz"), file("${params.project}_complete.vcf.gz.tbi") into appliedIndel_out

    """
    sentieon driver \\
    --thread_count $params.np_cpus \\
    -r $params.reference \\
    --algo ApplyVarCal \\
    --vcf $snp_recal_file \\
    --var_type INDEL \\
    --recal $recalVCF \\
    --tranches_file $IndelTranch \\
    --sensitivity 90 \\
    "${params.project}_complete.vcf.gz"
    """
}
appliedIndel_out.into { vcfStats_in; annotateVCF_in }

process finalStats {
    tag { "${params.project}" }

    publishDir "${params.vcfstats}", mode: 'copy'

    input:
    set file(vcf), file(index) from vcfStats_in

    output:
    file("${params.project}_FinalVCF.stats" )
    val 'complete' into finalStats

    shell:
    '''
    bcftools stats --threads !{params.np_cpus} !{vcf} > "!{params.project}_FinalVCF.stats"
    '''
}
process annotateFinalVCF {
    tag { "$params.project" }

    publishDir "${params.annotatedVCF}", mode: 'copy'

    input:
    set file(vcf), file(index) from annotateVCF_in

    output:
    file("${params.project}_complete.hg19_multianno.vcf")

    shell:
    '''
    table_annovar.pl \\
        !{vcf} \\
        $HOME/modules/annovar/humandb/ \\
        --buildver hg19 \\
        --out !{params.project}_complete \\
        --remove \\
        --protocol refGene,ensGene,knownGene,dbnsfp30a,gnomad_genome,exac03 \\
        --operation g,g,g,f,f,f \\
        --nastring . \\
        --vcfinput
    '''
}

fastqc_done
    .unique()
    .toList()
    .set { allFastqc }

samStats_done
    .unique()
    .toList()
    .set { allSamStats }

samFlagstat_done
    .unique()
    .toList()
    .set { allFlagStats }

coverage_done
    .unique()
    .toList()
    .set { allCoverage }

graph_done
    .unique()
    .toList()
    .set { allGraph }

process multiqc {
    tag { "$params.project" }

    publishDir "${params.multiqc}", mode: 'copy'

    input:
    val 'complete' from finalStats
    val 'complete' from allFastqc
    val 'complete' from allSamStats
    val 'complete' from allFlagStats
    val 'complete' from allCoverage
    val 'complete' from allGraph

    output:
    file("${params.project}.multiqc.report.html")

    """
    multiqc ${params.complete} --force --no-data-dir --filename "${params.project}.multiqc.report"
    """
}

*/
workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"

}
