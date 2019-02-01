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
log.info "=================================="

gzippedFastqs_ch = Channel.fromPath( '/uufs/chpc.utah.edu/common/home/u6013142/projects/eGFR/nextflow_variant_discovery/practiceData/double/*.fastq.gz' )
//gzippedFastqs_ch = Channel.fromPath( '/uufs/chpc.utah.edu/common/home/u6013142/projects/eGFR/nextflow_variant_discovery/realData/*.fastq.gz' )

process unzipFastqs {
    tag { sample_id + "_" + read_num }

    input:
    file gzfq_file from gzippedFastqs_ch

    output:
    set val(sample_id), val(read_num), file("${sample_id}_${read_num}.fastq") into fastqs_ch

    script:
    // Get sample name, samples are split into forward and reverse reads
    fileName = gzfq_file.toString()
    fileName_v2 = fileName.replaceAll("LU01-", "LU01_")
    sample_id = fileName_v2.tokenize('.')[0].tokenize('_')[1]
    read_num = fileName_v2.tokenize('.')[0].tokenize('_')[3]
    """
    gunzip -c ${gzfq_file} > ${sample_id}_${read_num}.fastq
    """
}

fastqs_ch
    .groupTuple()
    .set { fastp_in_ch }

process runFastp {
    tag { sample_id }

    publishDir "${params.fastp}", mode: 'copy', pattern: '*.html'

    input:
    set val(sample_id), val(read_nums), file(fq_files) from fastp_in_ch

    output:
    set val(sample_id), val(read1), file("${sample_id}_${read1}.trimmed.fastq.gz"), val(sample_id), val(read2), file("${sample_id}_${read2}.trimmed.fastq.gz") into fastp_out_ch
    file("${sample_id}.fastp.report.html")

    shell:
    if (fq_files[0].toString().contains("_R1.")) {
        fq1 = fq_files[0]
        fq2 = fq_files[1]
    } else {
        fq2 = fq_files[0]
        fq1 = fq_files[1]
    } 

    if (read_nums[0] == "R1") {
        read1 = read_nums[0]
        read2 = read_nums[1]
    } else {
        read1 = read_nums[1]
        read2 = read_nums[0]
    }
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
    --html "!{sample_id}.fastp.report.html" \\
    '''
}

fastp_out_ch
    .flatten()
    .collate( 3 )
    .set { collated_ch }

collated_ch.into { fastqc_in_ch; group_fq_ch }

process runFastqc {
    tag { sample_id + "_" + read_num }

    input:
    set val(sample_id), val(read_num), file(fq_file) from fastqc_in_ch

    output:
    val 'complete' into fastqc_out_ch
    
    shell:
    '''
    fastqc !{fq_file} -o !{params.fastqc} -t !{params.kp_cpus} 
    '''
}

group_fq_ch
    .groupTuple()
    .set { bwa_in_ch }
    
process BWA {
    tag { sample_id }

    input:
    set val(sample_id), file(read_nums), file(fq_files) from bwa_in_ch

    output:
    set sample_id, file("${sample_id}.sorted.bam"), file("${sample_id}.sorted.bam.bai") into bwa_out_ch

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
    --sam2bam -i - \\
    '''
}

process dedup {
    tag { sample_id }

    input:
    set sample_id, file(bam), file(index) from bwa_out_ch

    output:
    set sample_id, file("${sample_id}.dedup.bam"), file("${sample_id}.dedup.bam.bai") into dedup_out_ch

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
    set sample_id, file(deduped), file(index) from dedup_out_ch

    output:
    set sample_id, file("${deduped.baseName}.realigned.bam"), file("${deduped.baseName}.realigned.bam.bai") into realigner_out_ch

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
    set sample_id, file(realigned), file(index) from realigner_out_ch

    output:
    set sample_id, file("${realigned.baseName}.realigned.bqsr.bam"), file("${realigned.baseName}.realigned.bqsr.bam.bai"), file("${realigned.baseName}.recal_data_table") into bqsr_out_ch
    set sample_id, file("${realigned.baseName}.recal_data_table"), file("${realigned.baseName}.recal_data.table.post") into bqsr_graph_out_ch

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

bqsr_out_ch.into{ haplotyper_in_ch; stats_in_ch; flagstat_in_ch ; coverageMetrics_in_ch }

process samStats {
    tag { sample_id }

    publishDir "${params.stats}", mode: 'copy'

    input:
    set sample_id, file(bam), file(index), file(recal_table) from stats_in_ch

    output:
    file("${sample_id}.stats")
    val 'complete' into samStats_done_ch

    """
    samtools stats ${bam} > "${sample_id}.stats"
    """
}

process samFlagstat {
    tag { sample_id }

    publishDir "${params.stats}", mode: 'copy'

    input:
    set sample_id, file(bam), file(index), file(recal_table) from flagstat_in_ch

    output:
    file("${sample_id}.flagstat")
    val 'complete' into samFlagstat_done_ch

    """
    samtools flagstat ${bam} > "${sample_id}.flagstat"
    """
}

process coverageMetrics {
    tag { sample_id }

    publishDir "${params.coverage}", mode: 'copy'

    input:
    set sample_id, file(bam), file(index), file(recal_table) from coverageMetrics_in_ch

    output:
    file("${sample_id}.sample_summary")
    val 'complete' into coverage_out_ch
    
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

process graphBQSR {
    tag { sample_id }

    publishDir "${params.bqsr}", mode: 'copy'

    input:
    set sample_id, file(table), file(post) from bqsr_graph_out_ch

    output:
    file("${sample_id}.bqsr.pdf")
    val 'complete' into graph_out_ch

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

process haplotyper {
    tag { sample_id }

    publishDir "${params.gvcf}", mode: 'copy'

    input:
    set sample_id, file(recalbam), file(index), file(recalTable) from haplotyper_in_ch

    output:
    set file("${sample_id}.g.vcf.gz"), file("${sample_id}.g.vcf.gz.tbi") into haplotyper_out_ch

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
    "!{sample_id}.g.vcf" \\

    bgzip -c "!{sample_id}.g.vcf" > "!{sample_id}.g.vcf.gz"
    tabix -p vcf "!{sample_id}.g.vcf.gz"
    '''
} 

haplotyper_out_ch
    .toSortedList()
    .transpose()
    .first()
    .set { gvcfCollect_ch }

process gvcfTyper {
    tag { "$params.project" }

    input:
    val gvcfs from gvcfCollect_ch

    output:
    file ("${params.project}.g.vcf.gz") into typed_ch

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

typed_ch
    .unique()
    .toList()
    .set { chrTyped_ch }

process mergeGVCFs {
    tag { "$params.project" }

    input:
    val chrFiles from chrTyped_ch

    output:
    set file("${params.project}_merged.vcf.gz"), file("${params.project}_merged.vcf.gz.tbi") into merged_ch

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
    set file(merged_vcf), file(index) from merged_ch

    output:
    set file(merged_vcf), file(index), file("${params.project}_recal.tranches.snp"), file("${params.project}_recal.snp.vcf.gz"), file("${params.project}_recal.snp.vcf.gz.tbi") into VarCalSNP_ch

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
    set file(merged_vcf), file(merged_index), file(snpTranchFile), file(recalVCF), file(index) from VarCalSNP_ch

    output:
    set file("${params.project}_applyRecal.snp.vcf.gz"), file("${params.project}_applyRecal.snp.vcf.gz.tbi") into applyedSNP_ch

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
    set file(snp_recal_file), file(index) from applyedSNP_ch

    output:
    set file(snp_recal_file), file(index), file("${params.project}_recal.tranches.indel"), file("${params.project}_recal.indel.vcf.gz"), file("${params.project}_recal.indel.vcf.gz.tbi") into VarCalINDEL_ch

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
    set file(snp_recal_file), file(recal_index), file(IndelTranch), file(recalVCF), file(index) from VarCalINDEL_ch

    output:
    set file("${params.project}_complete.vcf.gz"), file("${params.project}_complete.vcf.gz.tbi") into applyedINDEL_ch

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
applyedINDEL_ch.into { finalVCF_ch; vcfStats_ch }

fastqc_out_ch
    .unique()
    .toList()
    .set { allFastqc_ch }

samStats_done_ch
    .unique()
    .toList()
    .set { allSamStats_ch }

samFlagstat_done_ch
    .unique()
    .toList()
    .set { allFlagStats_ch }

coverage_out_ch
    .unique()
    .toList()
    .set { allCoverage_ch }

graph_out_ch
    .unique()
    .toList()
    .set { allGraph_ch }

process finalStats {
    tag { "${params.project}" }

    publishDir "${params.vcfstats}", mode: 'copy'

    input:
    set file(vcf), file(index) from vcfStats_ch

    output:
    file("${params.project}_FinalVCF.stats" )
    val 'complete' into finalStats_ch

    shell:
    '''
    bcftools stats --threads !{params.np_cpus} !{vcf} > "!{params.project}_FinalVCF.stats"
    '''
}

process multiqc {
    tag { "$params.project" }

    publishDir "${params.multiqc}", mode: 'copy'

    input:
    val 'complete' from finalStats_ch
    val 'complete' from allFastqc_ch
    val 'complete' from allSamStats_ch
    val 'complete' from allFlagStats_ch
    val 'complete' from allCoverage_ch
    val 'complete' from allGraph_ch

    output:
    file("${params.project}.multiqc.report.html")

    """
    multiqc ${params.complete} --force --no-data-dir --filename "${params.project}.multiqc.report"
    """
}

process annotateFinalVCF {
    tag { "$params.project" }

    publishDir "${params.annotatedVCF}", mode: 'copy'

    input:
    set file(vcf), file(index) from finalVCF_ch

    output:
    file("${params.project}_complete.hg19_multianno.vcf")

    shell:
    //~/modules/annovar/table_annovar.pl \\
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

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}
