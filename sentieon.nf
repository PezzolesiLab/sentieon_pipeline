#!/usr/bin/nextflow

log.info """\
    ==================================
              Pezzolesi Lab           
            Sentieon Pipeline            
    ==================================
    
    -- With the following VQSR knowns --
       reference: $params.reference 
       dbsnp: $params.dbsnp 
       hapmap: $params.hapmap
       omni: $params.omni
       1G snps: $params.snp_1G
       1G indels: $params.indel_1G 
       mills indels: $params.indel_mills

    -- Set to the following parameters --
       hapmap: $params.hapmap_par 
       omni: $params.omni_par
       1G snps: $params.snp_1G_par 
       1G indels: $params.indel_1G_par 
       mills indels: $params.indel_mills_par

    -- Using the following region bed files --
       Targeted BED: $params.targetedBedFile
       Tiled BED: $params.tiledBedFile

    -- Software dependencies --
       March 21 2019: Not sure if including these executables in the ./bin will work
       annovar (currently pointing to my copy of the .pl script but the pezzolesi-group1 copy of the annovar DB)
       fastp/0.19.6
       multiqc/1.7
       cutadapt/1.6 or higher
       fastq-multx/1.3.1

    -- Modules Needed --
       bcftools/1.7
       bgzip/1.7
       samtools/1.7
       tabix/1.7
       sentieon/201711.05
    
    ==================================
"""

jointCalling = params.jointCallWithBackground
demuxing     = params.isDemuxNeeded
interval     = params.isIntervalNeeded
bamStart     = params.startingFromBams
isTargSeq    = params.isTargSeq
gvcfStart    = params.startingFromGvcfs
twoFqs       = params.inputTwoFastqs

if ( !gvcfStart ) {
    if ( !bamStart ) {
        if ( twoFqs ) {
            if ( demuxing ) {
            
                // specific to file names generated by the University of Utah DNA Sequencing Core Facility
                libraryReadLaneFq = Channel
                  .fromPath( "${params.dataDir}" )
                  .map { file ->
                      fNameExt = file.baseName
                      fName    = fNameExt.tokenize('.')[0]
                      library  = fName.tokenize('_')[0]
                      lane     = fName.tokenize('_')[6]
                      read_num = fName.tokenize('_')[7]
                      [ library, read_num, lane, file ]
                  }

                libraryReadLaneFq
                    .groupTuple()
                    .set { demux_in }

                process demuxFq {
                    tag { library + "_" + lane }

                    input:
                    set val(library), val(reads), val(lanes), val(fqs) from demux_in

                    output:
                    set val("${library}_${r1}_${lane}"), file("${library}_${r1}_${lane}.*.fastq.gz"), val("${library}_${r2}_${lane}"), file("${library}_${r2}_${lane}.*.fastq.gz") into demux_out

                    shell:

                    fq1  = fqs[0]
                    fq2  = fqs[1]
                    r1   = reads[0]
                    r2   = reads[1]
                    lane = lanes[0]
                    
                    '''
                    fastq-multx \\
                    -x \\
                    -b \\
                    -B !{params.barcodeFile} \\
                    !{fq1} \\
                    !{fq2} \\
                    -o !{library}_!{r1}_!{lane}.%.fastq.gz \\
                    -o !{library}_!{r2}_!{lane}.%.fastq.gz
                    '''
                }

                //.filter { ( it[1] =~ /bc[0-3][0-9]|bc4[0-8]/ & it[1] =~ /15686/ ) | ( it[1] =~ /bc49|bc[5-8][0-9]|bc9[0-6]/ & it[1] =~ /15887/ ) }
                //.filter { it[0] =~ /bc[0][1-2]/ }
                demux_out
                  .transpose()
                  .filter { ( it.toString() =~ /15686\w*\.(bc[0-3][0-9]|bc4[0-8])|15887\w*\.(bc49|bc[5-9][0-9]|bc9[0-6])|15942X1\w*\.(bc[0-3][0-9]|bc4[0-8])|15942X2\w*\.(bc49|bc[5-9][0-9]|bc9[0-6])/ ) }
                  .set{ trim_in }

                process trimReads {
                  tag { library1 + "_" + barcode1 }

                  input:
                  set val(library_read_lane1), file(fq_file1), val(library_read_lane2), file(fq_file2) from trim_in

                  output:
                  set val("${library1}_${r1}_${lane1}_${barcode1}"), file("${library1}_${r1}_${lane1}_${barcode1}.nobc.fastq.gz"), val(r1), file("${library2}_${r2}_${lane2}_${barcode2}.nobc.fastq.gz"), val(r2) into trim_out

                  shell:

                  fName1   = fq_file1.baseName
                  barcode1 = fName1.tokenize('.')[1]
                  library1 = library_read_lane1.tokenize('_')[0]
                  r1       = library_read_lane1.tokenize('_')[1]
                  lane1    = library_read_lane1.tokenize('_')[2]

                  fName2   = fq_file2.baseName
                  barcode2 = fName2.tokenize('.')[1]
                  library2 = library_read_lane2.tokenize('_')[0]
                  r2       = library_read_lane2.tokenize('_')[1]
                  lane2    = library_read_lane2.tokenize('_')[2]

                  '''
                  cutadapt \\
                  -u 6 \\
                  -U 6 \\
                  -j !{params.cpus_left} \\
                  -o "!{library1}_!{r1}_!{lane1}_!{barcode1}.nobc.fastq.gz" \\
                  -p "!{library2}_!{r2}_!{lane2}_!{barcode2}.nobc.fastq.gz" \\
                  !{fq_file1} \\
                  !{fq_file2}
                  '''
                }

                  //.filter { it[1] =~ /bc[0][1-2]/ }
                sampleNamesBarcodeKey = Channel
                  .fromPath( "${params.sampleKey}" )
                  .splitCsv(header: true)
                  .map { row -> tuple(row.Library + "_" + row.Read + "_L" + row.Lane + "_" + row.Barcode_ID, row.Sample_ID) }
                  //.filter { ( it.toString() =~ /15686\w*_(bc[0-3][0-9]|bc4[0-8])|15887\w*_(bc49|bc[5-9][0-9]|bc9[0-6])/ ) }
                  //.filter { (it[1] =~ /bc[0-3][0-9]|bc4[0-8]/ & it[1] =~ /15686/) | (it[1] =~ /bc49|bc[5-8][0-9]|bc9[0-6]/ & it[1] =~ /15887/ ) }

                trim_out.join(sampleNamesBarcodeKey).flatten().buffer(size:5, skip:1).set{ rename_in }

                rename_in
                  .map { renameList -> 
                      myfq1    = file(renameList[0])
                      r1       = renameList[1]
                      newid    = renameList[4]
                      dir1     = myfq1.getParent()
                      newName1 = "${dir1}/${newid}_${r1}.fastq.gz"
                      myfq1.copyTo(newName1)

                      myfq2    = file(renameList[2])
                      r2       = renameList[3]
                      dir2     = myfq2.getParent()
                      newName2 = "${dir2}/${newid}_${r2}.fastq.gz"
                      myfq2.copyTo(newName2)
                      [ newid, tuple(r1, r2), tuple(newName1, newName2) ]
                  }
                  .set { fastp_in }

            } else {

                // specific to UDDCS samples sequenced by WuXi
                idReadFq = Channel
                    .fromPath( "${params.dataDir}" )
                    .map { file ->
                        fName    = file.baseName
                        fName2   = fName.replaceAll("LU01-", "LU01_")
                        id       = fName2.tokenize('.')[0].tokenize('_')[1]
                        read_num = fName2.tokenize('.')[0].tokenize('_')[3]
                        [ id, read_num, file ]
                    }
                
                idReadFq
                    .groupTuple()
                    .set { fastp_in }
              
            }

            process runFastp {
                tag { sample_id }

                publishDir "${params.fastp}", mode: 'copy', pattern: '*.html'

                input:
                set val(sample_id), val(reads), val(fqs) from fastp_in

                output:
                set val(sample_id), val(r1), file("${sample_id}_${r1}.trimmed.fastq.gz"), val(sample_id), val(r2), file("${sample_id}_${r2}.trimmed.fastq.gz") into fastp_out
                file("${sample_id}.fastp.report.html")

                script:

                if (fqs[0].contains("_R1")) {
                    fq1 = fqs[0]
                    fq2 = fqs[1]
                    r1  = reads[0]
                    r2  = reads[1]
                } else {
                    fq1 = fqs[1]
                    fq2 = fqs[0]
                    r1  = reads[1]
                    r2  = reads[0]
                }

                if (demuxing) {

                    """
                    fastp \\
                    --thread "${params.single_cpus}" \\
                    --in1 "${fq1}" \\
                    --in2 "${fq2}" \\
                    --out1 "${sample_id}_${r1}.trimmed.fastq.gz" \\
                    --out2 "${sample_id}_${r2}.trimmed.fastq.gz" \\
                    --length_required "${params.minReadLength}" \\
                    --low_complexity_filter 5 \\
                    --html "${sample_id}.fastp.report.html"
                    """

                } else {

                    """
                    fastp \\
                    --thread "${params.single_cpus}" \\
                    --in1 ${fq1} \\
                    --in2 ${fq2} \\
                    --out1 "${sample_id}_${r1}.trimmed.fastq.gz" \\
                    --out2 "${sample_id}_${r2}.trimmed.fastq.gz" \\
                    --length_required "${params.minReadLength}" \\
                    --low_complexity_filter 5 \\
                    --detect_adapter_for_pe \\
                    --html "${sample_id}.fastp.report.html"
                    """

                }
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

                shell:
                '''
                fastqc !{fq_file} -o !{params.fastqc} -t !{params.cpus_left} 
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
                sync
                sleep 5
                """
            }

            gunzip_fq_out
                .groupTuple()
                .set { bwa_in }
                
            process BWA {
                tag { sample_id }
                echo true 
                beforeScript 'export MODULEPATH=$MODULEPATH:/scratch/ucgd/serial/common/modulefiles/notchpeak.peaks'
                cache 'deep'

                input:
                set val(sample_id), file(reads), file(fqs) from bwa_in

                output:
                set sample_id, file("${sample_id}.sorted.bam"), file("${sample_id}.sorted.bam.bai") into bwa_out

                shell:

                fq1 = fqs[0]
                fq2 = fqs[1]

                '''
                export RG=$(pezzAlign !{fq1})
                echo $MODULEPATH
                module load sentieon/201711.05

                ( bwa mem -M \\
                -R $RG \\
                -K 10000000 \\
                -t !{params.cpus_left} \\
                !{params.reference} \\
                !{fq1} \\
                !{fq2} \\
                || echo -n 'error' ) \\
                | sentieon util sort \\
                -o "!{sample_id}.sorted.bam" \\
                -t !{params.cpus_left} \\
                --sam2bam -i -
                '''
            }

        } else {
            // specific to WGS UDDCS samples
            idReadFq = Channel
                .fromPath( "${params.dataDir}" )
                .map { file ->
                    fName = file.baseName
                    id = fName.tokenize('.')[0].tokenize('_')[0]
                    [ id, file ]
                }

            idReadFq.into { fastp_in }
            
            process runFastpOneFq {
                tag { sample_id }

                publishDir "${params.fastp}", mode: 'copy', pattern: '*.html'

                input:
                set val(sample_id), val(fq) from fastp_in

                output:
                set val(sample_id), file("${sample_id}.trimmed.fastq.gz") into fastp_out
                file("${sample_id}.fastp.report.html")

                script:

                """
                fastp \\
                --thread "${params.single_cpus}" \\
                --in ${fq} \\
                --out "${sample_id}.trimmed.fastq.gz" \\
                --length_required 25 \\
                --low_complexity_filter 5 \\
                --detect_adapter_for_pe \\
                --html "${sample_id}.fastp.report.html"
                """
            }

            fastp_out.into { fastqc_in; gunzip_fq_in }

            process runFastqcOneFq {
                tag { sample_id }

                input:
                set val(sample_id), file(fq_file) from fastqc_in

                //output:
                //val 'fqc_complete' into fastqc_done
                
                shell:
                '''
                fastqc !{fq_file} -o !{params.fastqc} -t !{params.cpus_left} 
                '''
            }

            process unzipFastqsOneFq {
                tag { sample_id }

                input:
                set val(sample_id), file(fq_file) from gunzip_fq_in 

                output:
                set val(sample_id), file("${sample_id}.fastq") into gunzip_fq_out

                """
                gunzip -c ${fq_file} > ${sample_id}.fastq
                sync
                sleep 5
                """
            }

            process BWAOneFq {
                tag { sample_id }
                echo true 
                beforeScript 'export MODULEPATH=$MODULEPATH:/scratch/ucgd/serial/common/modulefiles/notchpeak.peaks'
                cache 'deep'

                input:
                set val(sample_id), file(fq) from bwa_in

                output:
                set sample_id, file("${sample_id}.sorted.bam"), file("${sample_id}.sorted.bam.bai") into bwa_out

                shell:

                fq1 = fqs[0]

                '''
                export RG=$(pezzAlign !{fq})
                echo $MODULEPATH
                module load sentieon/201711.05

                ( bwa mem -M \\
                -R $RG \\
                -K 10000000 \\
                -t !{params.cpus_left} \\
                !{params.reference} \\
                !{fq1} \\
                || echo -n 'error' ) \\
                | sentieon util sort \\
                -o "!{sample_id}.sorted.bam" \\
                -t !{params.cpus_left} \\
                --sam2bam -i -
                '''
            }
        }

        process dedup {
            tag { sample_id }
            beforeScript 'export MODULEPATH=$MODULEPATH:/scratch/ucgd/serial/common/modulefiles/notchpeak.peaks'
            
            input:
            set sample_id, file(bam), file(index) from bwa_out

            output:
            set sample_id, file("${sample_id}.dedup.bam"), file("${sample_id}.dedup.bam.bai") into dedup_out

            shell:
            '''
            module load sentieon/201711.05
            sentieon driver \\
            -t !{params.cpus_left} \\
            -i !{bam} \\
            --algo LocusCollector \\
            --fun score_info "!{sample_id}.score.txt" \\

            sentieon driver \\
            -t !{params.cpus_left} \\
            -i !{bam} \\
            --algo Dedup \\
            --score_info "!{sample_id}.score.txt" \\
            --metrics "!{sample_id}.metric.txt" \\
            "!{sample_id}.dedup.bam"
            '''
        }

        process indelRealigner {
            tag { sample_id }
            beforeScript 'export MODULEPATH=$MODULEPATH:/scratch/ucgd/serial/common/modulefiles/notchpeak.peaks'

            input:
            set sample_id, file(deduped), file(index) from dedup_out

            output:
            set sample_id, file("${deduped.baseName}.realigned.bam"), file("${deduped.baseName}.realigned.bam.bai") into realign_out

            shell:
            '''
            module load sentieon/201711.05
            sentieon driver \\
            -t !{params.cpus_left} \\
            -r !{params.reference} \\
            -i !{deduped} \\
            --algo Realigner \\
            -k !{params.indel_mills} \\
            -k !{params.indel_1G} \\
            "!{deduped.baseName}.realigned.bam"
            '''
        }

        process BQSR {
            tag { sample_id }
            beforeScript 'export MODULEPATH=$MODULEPATH:/scratch/ucgd/serial/common/modulefiles/notchpeak.peaks'

            publishDir "${params.bam}", mode: 'copy', pattern: '*.{bam,bai}'

            input:
            set sample_id, file(realigned), file(index) from realign_out

            output:
            set sample_id, file("${realigned.baseName}.realigned.bqsr.bam"), file("${realigned.baseName}.realigned.bqsr.bam.bai"), file("${realigned.baseName}.recal_data_table") into bqsr_out
            set sample_id, file("${realigned.baseName}.recal_data_table"), file("${realigned.baseName}.recal_data.table.post") into bqsr_graph_out

            shell:
            '''
            module load sentieon/201711.05
            sentieon driver \\
            -t !{params.cpus_left} \\
            -r !{params.reference} \\
            -i !{realigned} \\
            --algo QualCal \\
            -k !{params.dbsnp} \\
            "!{realigned.baseName}.recal_data_table"

            sentieon driver \\
            -t !{params.cpus_left} \\
            -r !{params.reference} \\
            -i !{realigned} \\
            -q "!{realigned.baseName}.recal_data_table" \\
            --algo QualCal \\
            -k !{params.dbsnp} \\
            "!{realigned.baseName}.recal_data.table.post" \\
            --algo ReadWriter "!{realigned.baseName}.realigned.bqsr.bam"
            '''
        }

        bqsr_out.into { stats_in; flagstat_in; coverageMetrics_in; haplotyper_in }

        process graphBQSR {
            tag { sample_id }
            beforeScript 'export MODULEPATH=$MODULEPATH:/scratch/ucgd/serial/common/modulefiles/notchpeak.peaks'

            publishDir "${params.bqsr}", mode: 'copy'

            input:
            set sample_id, file(table), file(post) from bqsr_graph_out

            output:
            file("${sample_id}.bqsr.pdf")
            //val 'graphbqsr_complete' into bqsrGraph_done

            shell:
            '''
            module load sentieon/201711.05
            sentieon driver \\
            -t !{params.cpus_left} \\
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
            beforeScript 'export MODULEPATH=$MODULEPATH:/scratch/ucgd/serial/common/modulefiles/notchpeak.peaks'

            input:
            set sample_id, file(bam), file(index), file(recal_table) from stats_in

            output:
            file("${sample_id}.stats")
            //val 'ss_complete' into samStats_done

            """
            module load samtools/1.7
            samtools stats ${bam} > "${sample_id}.stats"
            """
        }

        process samFlagstat {
            tag { sample_id }

            publishDir "${params.stats}", mode: 'copy'
            beforeScript 'export MODULEPATH=$MODULEPATH:/scratch/ucgd/serial/common/modulefiles/notchpeak.peaks'

            input:
            set sample_id, file(bam), file(index), file(recal_table) from flagstat_in

            output:
            file("${sample_id}.flagstat")
            val 'sfs_complete' into samFlagstat_done

            """
            module load samtools/1.7
            samtools flagstat ${bam} > "${sample_id}.flagstat"
            """
        }

        process coverageMetrics {
            tag { sample_id }
            beforeScript 'export MODULEPATH=$MODULEPATH:/scratch/ucgd/serial/common/modulefiles/notchpeak.peaks'

            publishDir "${params.coverage}", mode: 'copy'

            input:
            set sample_id, file(bam), file(index), file(recal_table) from coverageMetrics_in

            output:
            file("${sample_id}.sample_summary")
            
            shell:
            if (interval) {
                '''
                module load sentieon/201711.05
                sentieon driver \\
                -t !{params.cpus_left} \\
                -i !{bam} \\
                --interval !{params.targetedBedFile} \\
                -r !{params.reference} \\
                --algo CoverageMetrics \\
                --partition sample \\
                --omit_base_output \\
                --omit_interval_stat \\
                --omit_locus_stat \\
                "!{sample_id}"
                '''
            } else {
                '''
                module load sentieon/201711.05
                sentieon driver \\
                -t !{params.cpus_left} \\
                -i !{bam} \\
                -r !{params.reference} \\
                --algo CoverageMetrics \\
                --partition sample \\
                --omit_base_output \\
                --omit_locus_stat \\
                "!{sample_id}"
                '''
            }
        }

        process haplotyper {
            tag { sample_id }
            beforeScript 'export MODULEPATH=$MODULEPATH:/scratch/ucgd/serial/common/modulefiles/notchpeak.peaks'

            publishDir "${params.gvcf}", mode: 'copy'

            input:
            set sample_id, file(recalbam), file(index), file(recalTable) from haplotyper_in

            output:
            set file("${sample_id}.g.vcf.gz"), file("${sample_id}.g.vcf.gz.tbi") into haplotyper_out

            shell:
            '''
            module load sentieon/201711.05
            sentieon driver \\
            -t !{params.cpus_left} \\
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

    } else {

        sampleid_bam = Channel
          .fromPath( "${params.dataDir}" )
          .map { bam ->
              sampleidExt = bam.baseName
              sample_id = sampleidExt.tokenize('.')[0]
              [ sample_id, bam ]
          }
        sampleid_bam.into { stats_in; flagstat_in; coverageMetrics_in; haplotyper_in }

        process samStatsFromBam {
            tag { sample_id }

            publishDir "${params.stats}", mode: 'copy'
            beforeScript 'export MODULEPATH=$MODULEPATH:/scratch/ucgd/serial/common/modulefiles/notchpeak.peaks'

            //input:
            //set sample_id, file(bam), file(index), file(recal_table) from stats_in
            input:
            set sample_id, bam from stats_in

            output:
            file("${sample_id}.stats")
            //val 'ss_complete' into samStats_done

            """
            module load samtools/1.7
            samtools stats ${bam} > "${sample_id}.stats"
            """
        }

        process samFlagstatFromBam {
            tag { sample_id }

            publishDir "${params.stats}", mode: 'copy'
            beforeScript 'export MODULEPATH=$MODULEPATH:/scratch/ucgd/serial/common/modulefiles/notchpeak.peaks'

            //input:
            //set sample_id, file(bam), file(index), file(recal_table) from flagstat_in
            input:
            set sample_id, bam from flagstat_in

            output:
            file("${sample_id}.flagstat")
            //val 'sfs_complete' into samFlagstat_done

            """
            module load samtools/1.7
            samtools flagstat ${bam} > "${sample_id}.flagstat"
            """
        }

        process coverageMetricsFromBam {
            tag { sample_id }
            beforeScript 'export MODULEPATH=$MODULEPATH:/scratch/ucgd/serial/common/modulefiles/notchpeak.peaks'

            publishDir "${params.coverage}", mode: 'copy'

            //input:
            //set sample_id, file(bam), file(index), file(recal_table) from coverageMetrics_in
            input:
            set sample_id, bam from coverageMetrics_in

            output:
            file("${sample_id}.sample_summary")
            
            shell:

            if (interval) {
                '''
                module load sentieon/201711.05
                sentieon driver \\
                -t !{params.cpus_left} \\
                -i !{bam} \\
                --interval !{params.targetedBedFile} \\
                -r !{params.reference} \\
                --algo CoverageMetrics \\
                --partition sample \\
                --omit_base_output \\
                --omit_interval_stat \\
                --omit_locus_stat \\
                "!{sample_id}"
                '''
            } else {
                '''
                module load sentieon/201711.05
                sentieon driver \\
                -t !{params.cpus_left} \\
                -i !{bam} \\
                -r !{params.reference} \\
                --algo CoverageMetrics \\
                --partition sample \\
                --omit_base_output \\
                --omit_locus_stat \\
                "!{sample_id}"
                '''
            }
        }

        process haplotyperFromBam {
            tag { sample_id }
            beforeScript 'export MODULEPATH=$MODULEPATH:/scratch/ucgd/serial/common/modulefiles/notchpeak.peaks'

            publishDir "${params.gvcf}", mode: 'copy'

            input:
            set sample_id, bam from haplotyper_in

            output:
            set file("${sample_id}.g.vcf.gz"), file("${sample_id}.g.vcf.gz.tbi") into haplotyper_out

            shell:
            '''
            module load sentieon/201711.05
            sentieon driver \\
            -t !{params.cpus_left} \\
            -r !{params.reference} \\
            -i !{bam} \\
            --algo Haplotyper \\
            --emit_mode gvcf \\
            -d !{params.dbsnp} \\
            "!{sample_id}.g.vcf"

            bgzip -c "!{sample_id}.g.vcf" > "!{sample_id}.g.vcf.gz"
            tabix -p vcf "!{sample_id}.g.vcf.gz"
            '''
        } 
    }
}

if (jointCalling) {
    if (gvcfStart) {
        gVCF_channel = Channel
            .fromPath( "${params.dataDir}" )

        backgroundGVCF_channel = Channel
            .fromPath( "${params.pathToBackgroundGVCFs}" )

        gVCF_channel
            .concat(backgroundGVCF_channel)
            .toSortedList()
            .set { gvcfTyper_in }
    } else {
        //toCombineWithBackground = haplotyper_out
            //.set { toCombineWithBackground }
        toCombineWithBackground = haplotyper_out
            .filter{ ( it =~ /gz$/ ) }

        backgroundGVCF_channel = Channel
            .fromPath( "${params.pathToBackgroundGVCFs}" )

        toCombineWithBackground
            .concat(backgroundGVCF_channel)
            .toSortedList()
            .set { gvcfTyper_in }
    }
} else {
    if (gvcfStart) {
        gVCF_channel = Channel
                .fromPath( "${params.dataDir}" )
        gVCF_channel
            .toSortedList()
            .set { gvcfTyper_in }
    } else {
        haplotyper_out
            .toSortedList()
            .transpose()
            .first()
            .set { gvcfTyper_in }
    }
}

process gvcfTyper {
    tag { "$params.project" }
    beforeScript 'export MODULEPATH=$MODULEPATH:/scratch/ucgd/serial/common/modulefiles/notchpeak.peaks'

    publishDir "${params.vcf}", mode: 'copy'

    input:
    val gvcfs from gvcfTyper_in

    output:
    set file ("${params.project}_noVQSR.vcf.gz"), file("${params.project}_noVQSR.vcf.gz.tbi") into gvcfTyper_out

    shell:
    if (interval) {

        inputGVCFs = gvcfs.join(' -v ')

        '''
        module load sentieon/201711.05
        sentieon driver \\
        -t !{params.cpus_left} \\
        -r !{params.reference} \\
        --interval !{params.tiledBedFile} \\
        --algo GVCFtyper \\
        "!{params.project}_noVQSR.vcf.gz" \\
        -v !{inputGVCFs}
        '''

    } else {

        inputGVCFs = gvcfs.join(' -v ')

        '''
        module load sentieon/201711.05
        sentieon driver \\
        -t !{params.cpus_left} \\
        -r !{params.reference} \\
        --algo GVCFtyper \\
        "!{params.project}_noVQSR.vcf.gz" \\
        -v !{inputGVCFs}
        '''
    }
}

if (!isTargSeq) {
    process varCalSnp {
        tag { "$params.project" }
        beforeScript 'export MODULEPATH=$MODULEPATH:/scratch/ucgd/serial/common/modulefiles/notchpeak.peaks'

        input:
        set file(jointCalled_vcf), file(index) from gvcfTyper_out

        output:
        set file(jointCalled_vcf), file(index), file("${params.project}_recal.tranches.snp"), file("${params.project}_recal.snp.vcf.gz"), file("${params.project}_recal.snp.vcf.gz.tbi") into varCalSNP_out

        """
        module load sentieon/201711.05
        sentieon driver \\
        --thread_count $params.cpus_left \\
        -r $params.reference \\
        --algo VarCal \\
        --vcf $jointCalled_vcf \\
        --var_type SNP \\
        --tranche 90 \\
        --tranche 95 \\
        --tranche 99 \\
        --tranche 99.5 \\
        --tranche 99.6 \\
        --tranche 99.7 \\
        --tranche 99.8 \\
        --tranche 99.9 \\
        --tranche 99.95 \\
        --tranche 100 \\
        --tranches_file "${params.project}_recal.tranches.snp" \\
        --plot_file "${params.vqsr}/${params.project}_vqsrPlot.snp.txt" \\
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
        beforeScript 'export MODULEPATH=$MODULEPATH:/scratch/ucgd/serial/common/modulefiles/notchpeak.peaks'

        input:
        set file(jointCalled_vcf), file(merged_index), file(snpTranchFile), file(recalVCF), file(index) from varCalSNP_out

        output:
        set file("${params.project}_applyRecal.snp.vcf.gz"), file("${params.project}_applyRecal.snp.vcf.gz.tbi") into appliedSNP_out

        """
        module load sentieon/201711.05
        sentieon driver \\
        --thread_count $params.cpus_left \\
        -r $params.reference \\
        --algo ApplyVarCal \\
        --vcf $jointCalled_vcf \\
        --var_type SNP \\
        --recal $recalVCF \\
        --tranches_file $snpTranchFile \\
        --sensitivity "${params.truthSensitivitySnp}" \\
        "${params.project}_applyRecal.snp.vcf.gz"
        """
    }

    process varCalIndel {
        tag { "$params.project" }
        beforeScript 'export MODULEPATH=$MODULEPATH:/scratch/ucgd/serial/common/modulefiles/notchpeak.peaks'

        input:
        set file(snp_recal_file), file(index) from appliedSNP_out

        output:
        set file(snp_recal_file), file(index), file("${params.project}_recal.tranches.indel"), file("${params.project}_recal.indel.vcf.gz"), file("${params.project}_recal.indel.vcf.gz.tbi") into varCalIndel_out

        """
        module load sentieon/201711.05
        sentieon driver \\
        --thread_count $params.cpus_left \\
        -r $params.reference \\
        --algo VarCal \\
        --vcf $snp_recal_file \\
        --var_type INDEL \\
        --tranche 90 \\
        --tranche 95 \\
        --tranche 96 \\
        --tranche 97 \\
        --tranche 98 \\
        --tranche 99 \\
        --tranche 99.5 \\
        --tranche 100 \\
        --tranches_file "${params.project}_recal.tranches.indel" \\
        --plot_file "${params.vqsr}/${params.project}_vqsrPlot.indel.txt" \\
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
        beforeScript 'export MODULEPATH=$MODULEPATH:/scratch/ucgd/serial/common/modulefiles/notchpeak.peaks'

        publishDir "${params.vcf}", mode: 'copy'

        input:
        set file(snp_recal_file), file(recal_index), file(IndelTranch), file(recalVCF), file(index) from varCalIndel_out

        output:
        set file("${params.project}_VQSR.vcf.gz"), file("${params.project}_VQSR.vcf.gz.tbi") into appliedIndel_out

        """
        module load sentieon/201711.05
        sentieon driver \\
        --thread_count $params.cpus_left \\
        -r $params.reference \\
        --algo ApplyVarCal \\
        --vcf $snp_recal_file \\
        --var_type INDEL \\
        --recal $recalVCF \\
        --tranches_file $IndelTranch \\
        --sensitivity "${params.truthSensitivityIndel}" \\
        "${params.project}_VQSR.vcf.gz"
        """
    }
    appliedIndel_out.set { vcfStats_in }
} else {
    gvcfTyper_out.set { vcfStats_in }
}

process finalStats {
    tag { "${params.project}" }

    publishDir "${params.vcfstats}", mode: 'copy'
    beforeScript 'export MODULEPATH=$MODULEPATH:/scratch/ucgd/serial/common/modulefiles/notchpeak.peaks'

    input:
    set file(vcf), file(index) from vcfStats_in

    output:
    file("${params.project}_FinalVCF.stats" )
    val 'bcfstats_complete' into multiqc_greenlight

    shell:
    '''
    module load bcftools/1.7
    bcftools stats --threads !{params.cpus_left} !{vcf} > "!{params.project}_FinalVCF.stats"
    '''
}

process multiqc {
    tag { "$params.project" }

    publishDir "${params.multiqc}", mode: 'copy'

    input:
    val greenlight from multiqc_greenlight

    output:
    file("${params.project}.multiqc.report.html")

    script:
    if (bamStart) {
        """
        rm -r ${params.fastqc} ${params.fastp} ${params.bqsr} ${params.bam}
        multiqc ${params.complete} --force --no-data-dir --filename "${params.project}.multiqc.report"
        """
    } else if (gvcfStart) {
        """
        rm -r ${params.fastqc} ${params.fastp} ${params.bqsr} ${params.bam} ${params.gvcf}
        multiqc ${params.complete} --force --data-dir --filename "${params.project}.multiqc.report"
        """
    } else {
        """
        multiqc ${params.complete} --force --data-dir --filename "${params.project}.multiqc.report"
        """
    }
}

workflow.onComplete {
    log.info "Pipeline completed at: $workflow.complete"
    log.info "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}
