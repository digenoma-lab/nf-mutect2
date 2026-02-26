#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// ------------------------------------------------------------------
// Parameters (WGS hg38 defaults; SLURM profile in nextflow.config)
// ------------------------------------------------------------------
params.reads            = null          // Glob for paired FASTQs: "data/*_{R1,R2}.fastq.gz"
params.crams            = null          // Glob for existing CRAMs/BAMs
params.reference        = null          // hg38 FASTA
params.known_sites      = null          // gnomAD/common VCF for BQSR & contamination
params.intervals        = null          // BED or interval_list to seed SplitIntervals
params.scatter_count    = 50            // WGS shard count
params.interval_padding = 100
params.outdir           = "results"
params.tumor_sample     = null
params.normal_sample    = null
params.panel_of_normals = null          // Pre-built PON VCF
params.pon_crams        = null          // Glob of normal CRAMs to build PON
params.wgs              = true
params.sample_sheet     = null          // CSV with columns: subject,type,tumor_sample(optional),sample_id,fastq1,fastq2,cram
params.canonical_only   = true          // default restrict to autosomes + sex chromosomes
params.extra_intervals  = null          // BED/interval_list to add non-canonical contigs

// ------------------------------------------------------------------
// Processes
// ------------------------------------------------------------------

process GATK_CREATE_DICT {
    tag "$reference.baseName"
    publishDir "${params.outdir}/ref", mode: 'copy', saveAs: { "genome.dict" }

    input:
    path reference

    output:
    path "*.dict", emit: dict

    script:
    """
    gatk CreateSequenceDictionary -R ${reference} -O genome.dict
    """
}

process SAMTOOLS_FAIDX {
    tag "$reference.baseName"
    publishDir "${params.outdir}/ref", mode: 'copy'

    input:
    path reference

    output:
    path "*.fai", emit: fai

    script:
    """
    samtools faidx ${reference}
    """
}

process MAKE_CANONICAL_BED {
    tag "$dict.baseName"
    publishDir "${params.outdir}/intervals", mode: 'copy', saveAs: { "canonical_autosomes_sex.bed" }

    input:
    path dict

    output:
    path "canonical_autosomes_sex.bed", emit: bed

    script:
    """
    awk 'BEGIN{OFS="\\t"} /^@SQ/{
        sn=""; ln="";
        for(i=1;i<=NF;i++){
            if(\$i~/^SN:/){sn=substr(\$i,4)}
            if(\$i~/^LN:/){ln=substr(\$i,4)}
        }
        if(sn ~ /^(chr)?([1-9]|1[0-9]|2[0-2]|X|Y)$/){print sn,0,ln}
    }' ${dict} > canonical_autosomes_sex.bed
    """
}

process MERGE_INTERVALS {
    tag "merge_intervals"
    publishDir "${params.outdir}/intervals", mode: 'copy', saveAs: { "merged_intervals.bed" }

    input:
    path primary
    path extra optional: true

    output:
    path "merged_intervals.bed", emit: bed

    script:
    """
    cat ${primary} > merged_intervals.bed
    if [ -n "${extra}" ]; then
        cat ${extra} >> merged_intervals.bed
    fi
    sort -k1,1 -k2,2n merged_intervals.bed | uniq > merged_intervals.tmp && mv merged_intervals.tmp merged_intervals.bed
    """
}

process FASTQC {
    tag "$sampleId-mem"
    //label 'process_medium'
    publishDir "$params.outdir/QC/FASTQC", mode: "copy"



    input:
    tuple val(sampleId), val(part), file(read1), file(read2)

    output:
    path("${sampleId}-${part}.fastqc"), emit: fqc

    script:
    if(params.debug == true){
    """
    echo fastqc -o ${sampleId}-${part}.fastqc $read1 $read2
    mkdir -p ${sampleId}-${part}.fastqc
    touch ${sampleId}-${part}.fastqc/report.fastqc

    """
    } else{
    """
    mkdir -p ${sampleId}-${part}.fastqc
    fastqc -t $task.cpus -o ${sampleId}-${part}.fastqc $read1 $read2
    """
    }

}

//we do run bwa-mem2 or bwa mem
process BWAMEM {

    tag "$sampleId-mem"
    //label 'process_high'
    publishDir "$params.outdir/BWA", mode: "copy", pattern: '*.log.*'
    publishDir "$params.outdir/BWA/HLA", mode: "copy", pattern: '*.hla.all'

    input:
    tuple val(sampleId), val(part), file(read1), file(read2)

    output:
    tuple val("${sampleId}"), val("${part}"), file("${sampleId}-${part}.mkdup.cram"),file("${sampleId}-${part}.mkdup.cram.crai"), emit: bams
    path("${sampleId}-${part}.log.bwamem")
    path("${sampleId}-${part}.hla.all") , optional: true
    path("${sampleId}-${part}.log.hla") , optional: true

   script:
    def aln="bwa-mem2"
    //we define the aln tool
    if(params.aligner=="bwa"){
    aln="bwa"
    }
    if(params.debug == true){
    """
    echo "seqtk mergepe $read1 $read2 | ${aln} mem -p -t $task.cpus -R'@RG\\tID:${sampleId}-${part}\\tSM:${sampleId}\\tPL:ill' ${params.ref} - 2> ${sampleId}-${part}.log.bwamem | k8 ${params.alt_js} -p ${sampleId}-${part}.hla ${params.ref}.alt | samtools view -1 - > ${sampleId}-${part}.aln.bam"
    echo "run-HLA ${sampleId}-${part}.hla > ${sampleId}-${part}.hla.top 2> ${sampleId}-${part}.log.hla;"
    echo "touch ${sampleId}-${part}.hla.HLA-dummy.gt; cat ${sampleId}-${part}.hla.HLA*.gt | grep ^GT | cut -f2- > ${sampleId}-${part}.hla.all"
    echo "rm -f ${sampleId}-${part}.hla.HLA*;"
    touch ${sampleId}-${part}.mkdup.cram
    touch ${sampleId}-${part}.mkdup.cram.crai
    touch ${sampleId}-${part}.log.bwamem
    touch ${sampleId}-${part}.hla.all
    """
    }else{
    if(params.hla == true){
    """
    seqtk mergepe $read1 $read2 \\
        | ${aln} mem -p -t $task.cpus -R'@RG\\tID:${sampleId}-${part}\\tSM:${sampleId}\\tPL:ill' ${params.ref} - 2> ${sampleId}-${part}.log.bwamem \\
        | k8 ${params.alt_js} -p ${sampleId}-${part}.hla ${params.ref}.alt \
        | samtools view -Sb -  \
        | samtools fixmate -m - -  \
        | samtools sort -m 1G -@ ${task.cpus - 8 } -  \
        | samtools markdup -O CRAM  --write-index --reference ${params.ref} -@ ${task.cpus - 8} - ${sampleId}-${part}.mkdup.cram

   run-HLA ${sampleId}-${part}.hla > ${sampleId}-${part}.hla.top 2> ${sampleId}-${part}.log.hla;
     touch ${sampleId}-${part}.hla.HLA-dummy.gt; cat ${sampleId}-${part}.hla.HLA*.gt | grep ^GT | cut -f2- > ${sampleId}-${part}.hla.all;
     rm -f ${sampleId}-${part}.hla.HLA*;
    """
    }
    else if (params.alt == true){
     """
    seqtk mergepe $read1 $read2  \\
    | ${aln} mem -p -t $task.cpus  -R'@RG\\tID:${sampleId}-${part}\\tSM:${sampleId}\\tPL:ill' ${params.ref} - 2> ${sampleId}-${part}.log.bwamem \\
    | k8 ${params.alt_js} -p ${sampleId}-${part}.hla hs38DH.fa.alt \
    | samtools view -Sb -  \
    | samtools fixmate -m - -  \
    | samtools sort -m 1G -@ ${task.cpus - 8 } -  \
    | samtools markdup -O CRAM  --write-index --reference ${params.ref} -@ ${task.cpus - 8} - ${sampleId}-${part}.mkdup.cram

     """
    }else{
    //normal mapping mode
     """
      seqtk mergepe $read1 $read2 \\
    | ${aln} mem -p -t $task.cpus  -R'@RG\\tID:${sampleId}-${part}\\tSM:${sampleId}\\tPL:ill' ${params.ref} - 2> ${sampleId}-${part}.log.bwamem \\
    | samtools view -Sb -  \
    | samtools fixmate -m - -  \
    | samtools sort -@ ${task.cpus} -  \
    | samtools markdup -O CRAM  --write-index --reference ${params.ref} -@ ${task.cpus} - ${sampleId}-${part}.mkdup.cram

     """

    }
  }

}

//merge bams by sample
process MERGEB{

  tag "$sampleId-merge"
  publishDir "$params.outdir/CRAM", mode: "copy"
  container "oras://community.wave.seqera.io/library/bwa-mem2_instrain_multiqc_qualimap_samtools:850f96dbd001458f"

  input:
  tuple val(sampleId), val(parts), file(cramFiles), file(craiFiles)
  path reference

  output:
  tuple val(sampleId), file("${sampleId}.merged.cram"), file("${sampleId}.merged.cram.crai") ,emit: mbams

  script:
  def filesb = cramFiles instanceof List ? cramFiles : [cramFiles]
  if ( filesb.size() == 1 ) {
    if ( params.debug == true ) {
                """
                echo ln -s ${cramFiles[0]} ${sampleId}.merged.cram
                touch ${sampleId}.merged.cram.crai
                touch ${sampleId}.merged.cram
                """
            } else {
                """
                ln -s ${cramFiles[0]} ${sampleId}.merged.cram
                ln -s ${craiFiles[0]} ${sampleId}.merged.cram.crai
                """
            }
  }else{
  if(params.debug == true){
  """
    echo samtools merge --write-index --reference  ${reference} -O CRAM -@ $task.cpus -f ${sampleId}.merged.cram ${cramFiles}
    touch ${sampleId}.merged.cram
    touch ${sampleId}.merged.cram.crai
  """
  }else{
  """
  samtools merge --write-index --reference ${reference} -O CRAM -@ $task.cpus -f ${sampleId}.merged.cram ${cramFiles}
  """
  }
 }
}

process QUALIMAP {
    tag "${sampleId}"
    publishDir "${params.outdir}/qualimap", mode: 'copy'


    container "oras://community.wave.seqera.io/library/bwa-mem2_instrain_multiqc_qualimap_samtools:850f96dbd001458f"

    input:
    tuple val(sampleId), path(cram), path(crai)
    path reference
    path fai

    output:
    path "${sampleId}", emit: qualimap_results

    script:
    if(params.debug){
    """
    echo qualimap bamqc \
        -bam ${cram} \
        -outdir qualimap_results/${sampleId} \
        -nt ${task.cpus} \
        --java-mem-size=${task.memory.toGiga()}G
        mkdir ${sampleId}
    """
    }else{
    """
    # reference FASTA (necesario para descomprimir CRAM)
    REF=${reference}

     # ---------- 1)  Crear FIFO ----------
    fifo="${sampleId}.bam.pipe"
     mkfifo "\$fifo"

    # ---------- 2)  Convertir CRAM → BAM en segundo plano ----------
     samtools view -@ ${task.cpus - 2 } -T ${reference} -b ${cram} > "\$fifo" &
     sam_pid=\$!

     # ---------- 3)  Lanzar Qualimap leyendo del FIFO ----------
      qualimap bamqc \
        -bam    "\$fifo" \
        -outdir ${sampleId} \
        -nt     ${task.cpus} \
        --java-mem-size=${task.memory.toGiga()}G
       rc=\$?

     # ---------- 4)  Limpiar ----------
     kill "\$sam_pid" 2>/dev/null || true   # por si qualimap terminó antes
     rm -f "\$fifo"

    exit "\$rc"
    """
    }
}


process VALIDATE_CRAM {
    tag "$meta.id"

    input:
    tuple val(meta), path(cram)
    path reference

    output:
    tuple val(meta), path(cram), path("${cram}.crai"), emit: cram_crai

    script:
    """
    [ -f ${cram}.crai ] || samtools index ${cram}
    samtools quickcheck ${cram}
    ln -s ${cram}.crai ${cram}.crai || true
    """
}

process GATK_MARKDUPLICATES {
    tag "$meta.id"
    publishDir "${params.outdir}/align", mode: 'copy', pattern: "*_markdup.cram"

    input:
    tuple val(meta), path(cram), path(crai)
    path reference

    output:
    tuple val(meta), path("${meta.id}_markdup.cram"), path("${meta.id}_markdup.cram.crai"), emit: cram_crai
    path("${meta.id}_markdup_metrics.txt"), emit: metrics

    script:
    """
    gatk MarkDuplicates \\
      -I ${cram} \\
      -O ${meta.id}_markdup.cram \\
      -M ${meta.id}_markdup_metrics.txt \\
      --CREATE_INDEX true \\
      --VALIDATION_STRINGENCY LENIENT \\
      --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500
    """
}

process GATK_BASERECALIBRATOR {
    tag "$meta.id"

    input:
    tuple val(meta), path(cram), path(crai)
    path reference
    path known_sites

    output:
    tuple val(meta), path("${meta.id}_recal.table"), emit: table

    script:
    """
    gatk BaseRecalibrator \\
      -I ${cram} \\
      -R ${reference} \\
      --known-sites ${known_sites} \\
      -O ${meta.id}_recal.table
    """
}

process GATK_APPLYBQSR {
    tag "$meta.id"
    publishDir "${params.outdir}/align", mode: 'copy', pattern: "*_recal.cram"

    input:
    tuple val(meta), path(cram), path(crai)
    tuple val(meta2), path(table)
    path reference

    output:
    tuple val(meta), path("${meta.id}_recal.cram"), path("${meta.id}_recal.cram.crai"), emit: cram_crai

    script:
    """
    gatk ApplyBQSR \\
      -I ${cram} \\
      -R ${reference} \\
      --bqsr-recal-file ${table} \\
      -O ${meta.id}_recal.cram \\
      --create-output-bam-index true
    """
}

process GATK_SPLITINTERVALS {
    tag "$reference.baseName"
    publishDir "${params.outdir}/intervals", mode: 'copy'

    input:
    path reference
    path dict
    path base_intervals optional: true
    val scatter_count

    output:
    path("intervals/*"), emit: intervals

    script:
    def baseArg = base_intervals ? "-L ${base_intervals}" : ""
    """
    mkdir -p intervals
    gatk SplitIntervals \\
      -R ${reference} \\
      ${baseArg} \\
      --sequence-dictionary ${dict} \\
      --scatter-count ${scatter_count} \\
      --interval-padding ${params.interval_padding} \\
      --subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION \\
      -O intervals
    """
}

process GATK_MUTECT2_SCATTER {
    tag "$tumor_meta.id:${interval.baseName}"
    publishDir "${params.outdir}/variants/shards", mode: 'copy'

    input:
    tuple path(interval), val(tumor_meta), path(tumor_cram), path(tumor_crai), val(normal_meta), path(normal_cram), path(normal_crai)
    path reference
    path pon optional: true

    output:
    tuple val(tumor_meta), path("${tumor_meta.id}_${interval.baseName}.vcf.gz"), path("${tumor_meta.id}_${interval.baseName}.vcf.gz.tbi"), emit: vcf
    tuple val(tumor_meta), path("${tumor_meta.id}_${interval.baseName}_f1r2.tar.gz"), emit: f1r2
    path("${tumor_meta.id}_${interval.baseName}.stats"), emit: stats

    script:
    def ponArg = pon ? "--panel-of-normals ${pon}" : ""
    """
    gatk Mutect2 \\
      -R ${reference} \\
      -I ${tumor_cram} \\
      -I ${normal_cram} \\
      -tumor ${tumor_meta.sample} \\
      -normal ${normal_meta.sample} \\
      -L ${interval} \\
      ${ponArg} \\
      --f1r2-tar-gz ${tumor_meta.id}_${interval.baseName}_f1r2.tar.gz \\
      -O ${tumor_meta.id}_${interval.baseName}.vcf.gz
    """
}

process GATK_MUTECT2_SCATTER_TONLY {
    tag "$tumor_meta.id:${interval.baseName}"
    publishDir "${params.outdir}/variants/shards", mode: 'copy'

    input:
    tuple path(interval), val(tumor_meta), path(tumor_cram), path(tumor_crai)
    path reference
    path pon optional: true

    output:
    tuple val(tumor_meta), path("${tumor_meta.id}_${interval.baseName}.vcf.gz"), path("${tumor_meta.id}_${interval.baseName}.vcf.gz.tbi"), emit: vcf
    tuple val(tumor_meta), path("${tumor_meta.id}_${interval.baseName}_f1r2.tar.gz"), emit: f1r2
    path("${tumor_meta.id}_${interval.baseName}.stats"), emit: stats

    script:
    def ponArg = pon ? "--panel-of-normals ${pon}" : ""
    """
    gatk Mutect2 \\
      -R ${reference} \\
      -I ${tumor_cram} \\
      -tumor ${tumor_meta.sample} \\
      -L ${interval} \\
      ${ponArg} \\
      --f1r2-tar-gz ${tumor_meta.id}_${interval.baseName}_f1r2.tar.gz \\
      -O ${tumor_meta.id}_${interval.baseName}.vcf.gz
    """
}

process GATK_MERGEVCFS {
    tag "$meta.id"
    publishDir "${params.outdir}/variants", mode: 'copy', pattern: "*_merged.vcf.gz"

    input:
    tuple val(meta), path(vcf_list)
    path reference

    output:
    tuple val(meta), path("${meta.id}_merged.vcf.gz"), path("${meta.id}_merged.vcf.gz.tbi"), emit: vcf

    script:
    def inputs = vcf_list.collect { "-I ${it}" }.join(" ")
    """
    gatk MergeVcfs \\
      ${inputs} \\
      -R ${reference} \\
      -O ${meta.id}_merged.vcf.gz
    gatk IndexFeatureFile -I ${meta.id}_merged.vcf.gz
    """
}

process GATK_LEARNREADORIENTATIONMODEL {
    tag "$meta.id"

    input:
    tuple val(meta), path(f1r2_list)

    output:
    tuple val(meta), path("${meta.id}_read_orientation_model.tar.gz"), emit: model

    script:
    def inputs = f1r2_list.collect { "-I ${it}" }.join(" ")
    """
    gatk LearnReadOrientationModel \\
      ${inputs} \\
      -O ${meta.id}_read_orientation_model.tar.gz
    """
}

process GATK_GETPILEUPSUMMARIES {
    tag "$meta.id"

    input:
    tuple val(meta), path(cram), path(crai)
    path known_sites

    output:
    tuple val(meta), path("${meta.id}_pileups.table"), emit: table

    script:
    """
    gatk GetPileupSummaries \\
      -I ${cram} \\
      -V ${known_sites} \\
      -L ${known_sites} \\
      -O ${meta.id}_pileups.table
    """
}

process GATK_CALCULATECONTAMINATION {
    tag "$tumor_meta.id"

    input:
    tuple val(tumor_meta), path(tumor_pileups)
    path normal_pileups optional: true

    output:
    tuple val(tumor_meta), path("${tumor_meta.id}_contamination.table"), emit: contamination
    tuple val(tumor_meta), path("${tumor_meta.id}_segments.table"), emit: segments

    script:
    def matchedArg = normal_pileups ? "-matched ${normal_pileups}" : ""
    """
    gatk CalculateContamination \\
      -I ${tumor_pileups} \\
      ${matchedArg} \\
      -O ${tumor_meta.id}_contamination.table \\
      --tumor-segmentation ${tumor_meta.id}_segments.table
    """
}

process GATK_FILTERMUTECTCALLS {
    tag "$meta.id"
    publishDir "${params.outdir}/variants/filtered", mode: 'copy'

    input:
    tuple val(meta), path(vcf), path(tbi), path(contamination), path(segments) optional: true, path(orientation_model)
    path reference

    output:
    tuple val(meta), path("${meta.id}_filtered.vcf.gz"), path("${meta.id}_filtered.vcf.gz.tbi"), emit: vcf
    path("${meta.id}_filtering_stats.tsv"), emit: stats

    script:
    """
    gatk FilterMutectCalls \\
      -R ${reference} \\
      -V ${vcf} \\
      --contamination-table ${contamination} \\
      ${segments ? "--tumor-segmentation ${segments}" : ""} \\
      --ob-priors ${orientation_model} \\
      -O ${meta.id}_filtered.vcf.gz \\
      --filtering-stats ${meta.id}_filtering_stats.tsv
    """
}

process MULTIQC {
    publishDir "${params.outdir}/qc/multiqc", mode: 'copy'

    input:
    path('*')

    output:
    path("multiqc_report.html"), emit: report

    script:
    """
    multiqc .
    """
}

process GATK_CONTAMINATION {
    tag "$tumor_meta.id"

    input:
    tuple val(tumor_meta), path(tumor_cram), path(tumor_crai), val(normal_meta), path(normal_cram), path(normal_crai)
    path known_sites

    output:
    tuple val(tumor_meta), path("${tumor_meta.id}_contamination.table"), emit: contamination
    tuple val(tumor_meta), path("${tumor_meta.id}_segments.table"), emit: segments

    script:
    """
    gatk GetPileupSummaries -I ${tumor_cram} -V ${known_sites} -L ${known_sites} -O tumor_pileups.table
    gatk GetPileupSummaries -I ${normal_cram} -V ${known_sites} -L ${known_sites} -O normal_pileups.table
    gatk CalculateContamination \\
      -I tumor_pileups.table \\
      -matched normal_pileups.table \\
      -O ${tumor_meta.id}_contamination.table \\
      --tumor-segmentation ${tumor_meta.id}_segments.table
    """
}

process GATK_CONTAMINATION_TONLY {
    tag "$tumor_meta.id"

    input:
    tuple val(tumor_meta), path(tumor_cram), path(tumor_crai)
    path known_sites

    output:
    tuple val(tumor_meta), path("${tumor_meta.id}_contamination.table"), emit: contamination
    path("${tumor_meta.id}_segments.table"), emit: segments

    script:
    """
    gatk GetPileupSummaries -I ${tumor_cram} -V ${known_sites} -L ${known_sites} -O tumor_pileups.table
    gatk CalculateContamination \\
      -I tumor_pileups.table \\
      -O ${tumor_meta.id}_contamination.table \\
      --tumor-segmentation ${tumor_meta.id}_segments.table
    """
}

// PON creation
process GATK_MUTECT2_NORMAL {
    tag "$meta.id:${interval.baseName}"
    publishDir "${params.outdir}/pon/shards", mode: 'copy'

    input:
    tuple path(interval), val(meta), path(cram), path(crai)
    path reference

    output:
    tuple val(meta), path("${meta.id}_${interval.baseName}_pon.vcf.gz"), emit: vcf

    script:
    """
    gatk Mutect2 \\
      -R ${reference} \\
      -I ${cram} \\
      -tumor ${meta.sample} \\
      -L ${interval} \\
      --max-mnp-distance 0 \\
      -O ${meta.id}_${interval.baseName}_pon.vcf.gz
    """
}

process GATK_MERGE_PON_NORMAL {
    tag "$meta.id"
    publishDir "${params.outdir}/pon", mode: 'copy'

    input:
    tuple val(meta), path(vcf_list)
    path reference

    output:
    tuple val(meta), path("${meta.id}_pon_merged.vcf.gz"), emit: vcf

    script:
    def inputs = vcf_list.collect { "-I ${it}" }.join(" ")
    """
    gatk MergeVcfs ${inputs} -R ${reference} -O ${meta.id}_pon_merged.vcf.gz
    gatk IndexFeatureFile -I ${meta.id}_pon_merged.vcf.gz
    """
}

process GATK_CREATE_SOMATIC_PON {
    publishDir "${params.outdir}/pon", mode: 'copy'

    input:
    path vcf_list

    output:
    path "panel_of_normals.vcf.gz", emit: pon

    script:
    def inputs = vcf_list.collect { "-vcfs ${it}" }.join(" ")
    """
    gatk CreateSomaticPanelOfNormals \\
      ${inputs} \\
      -O panel_of_normals.vcf.gz
    """
}

// ------------------------------------------------------------------
// Main workflow
// ------------------------------------------------------------------

workflow {
    if (!params.sample_sheet && !params.reads && !params.crams)
        error "Provide --sample_sheet (preferred) or --reads/--crams"
    if (!params.reference) error "Provide --reference (hg38 FASTA)"
    if (!params.known_sites) error "Provide --known_sites VCF for BQSR/contamination"
    if (!params.sample_sheet && (!params.tumor_sample))
        error "Provide --tumor_sample (and --normal_sample if paired) when not using --sample_sheet"

    reference_ch    = Channel.fromPath(params.reference, checkIfExists: true).first()
    reference_ch.into { ref_dict; ref_fai; ref_bwa; ref_split; ref_mutect; ref_merge; ref_filter; ref_pon }

    known_sites_ch  = Channel.fromPath(params.known_sites, checkIfExists: true).first()
    known_sites_ch.into { ks_bqsr; ks_contam_t; ks_contam_n; ks_pon }
    ks_contam_t.into { ks_contam_paired; ks_contam_tonly }

    dict_ch         = GATK_CREATE_DICT(ref_dict).out.dict
    canonical_bed   = MAKE_CANONICAL_BED(dict_ch).out.bed
    fai_ch          = SAMTOOLS_FAIDX(ref_fai).out.fai
    bwa_idx_ch      = BWA_MEM2_INDEX(ref_bwa).out.index
    base_interval_seed = params.intervals ? Channel.fromPath(params.intervals, checkIfExists: true).first()
                                          : (params.canonical_only ? canonical_bed.first() : Channel.empty())
    if (params.extra_intervals) {
        extra_intervals_ch = Channel.fromPath(params.extra_intervals, checkIfExists: true).first()
        base_interval_seed = MERGE_INTERVALS(base_interval_seed, extra_intervals_ch).out.bed
    }
    base_interval_seed.into { intervals_seed_main; intervals_seed_pon }

    // Input handling (sample sheet preferred)
    Channel crams_valid
    Channel fastqc_zip = Channel.empty()

    if (params.sample_sheet) {
        sample_rows = Channel.fromPath(params.sample_sheet, checkIfExists: true)
            .splitCsv(header:true)
            .map { row ->
                def meta = [
                    id     : row.sample_id ?: row.subject ?: row.type,
                    sample : row.sample_id ?: row.subject ?: row.type,
                    subject: row.subject ?: row.sample_id ?: row.type,
                    type   : row.type?.toLowerCase()
                ]
                def fastqs = (row.fastq1 && row.fastq2) ? [row.fastq1, row.fastq2] : null
                def cram = row.cram ?: null
                [meta, fastqs, cram]
            }
            .filter { meta, fastqs, cram -> meta.type in ['tumor','normal'] }

        fastq_samples = sample_rows.filter { meta, fastqs, cram -> fastqs }.map { meta, fastqs, cram -> [meta, fastqs] }
        cram_samples  = sample_rows.filter { meta, fastqs, cram -> cram }.map { meta, fastqs, cram -> [meta, cram] }

        if (!fastq_samples.isEmpty()) {
            FASTQC(fastq_samples)
            fastqc_zip = FASTQC.out.zip
            crams_from_fastq = BWA_MEM2_ALIGN(fastq_samples, ref_mutect, bwa_idx_ch, fai_ch).out.cram_crai
        } else {
            crams_from_fastq = Channel.empty()
        }

        crams_from_cram = cram_samples.isEmpty() ? Channel.empty() : VALIDATE_CRAM(cram_samples, ref_mutect).out.cram_crai
        crams_valid = crams_from_fastq.mix(crams_from_cram)

    } else if (params.crams) {
        crams_in = Channel.fromPath(params.crams, checkIfExists: true).map { cram ->
            def meta = [id: cram.baseName.replaceAll(/\\.cram$|\\.bam$/, ''), sample: cram.baseName.replaceAll(/\\.cram$|\\.bam$/, ''), subject: cram.baseName.replaceAll(/\\.cram$|\\.bam$/, ''), type: 'tumor']
            [meta, cram]
        }
        crams_valid = VALIDATE_CRAM(crams_in, ref_mutect).out.cram_crai
    } else {
        reads_ch = Channel.fromFilePairs(params.reads, flat: true, checkIfExists: true).map { sample, reads ->
            def meta = [id: sample, sample: sample, subject: sample, type: sample == params.normal_sample ? 'normal' : 'tumor']
            [meta, reads]
        }
        FASTQC(reads_ch)
        fastqc_zip = FASTQC.out.zip
        crams_valid = BWA_MEM2_ALIGN(reads_ch, ref_mutect, bwa_idx_ch, fai_ch).out.cram_crai
    }

    // Preprocess
    md = GATK_MARKDUPLICATES(crams_valid, ref_mutect)
    bqsr = GATK_BASERECALIBRATOR(md.out.cram_crai, ref_mutect, ks_bqsr)
    recals = GATK_APPLYBQSR(md.out.cram_crai, bqsr.out.table, ref_mutect)
    final_crams = recals.out.cram_crai

    // Pairing
    Channel subject_pairs
    if (params.sample_sheet) {
        grouped = final_crams.groupBy { meta, cram, crai -> meta.subject }
        subject_pairs = grouped.map { subject, items ->
            def tumor = items.find { it[0].type == 'tumor' }
            if (!tumor) return null
            def normal = items.find { it[0].type == 'normal' }
            [subject, tumor, normal]
        }.filter { it != null }
    } else {
        tumor_cram = final_crams.filter { meta, cram, crai -> meta.sample == params.tumor_sample }
        normal_cram = final_crams.filter { meta, cram, crai -> params.normal_sample ? meta.sample == params.normal_sample : false }
        subject_pairs = tumor_cram.map { t -> [t[0].sample, t, null] }
        if (!normal_cram.isEmpty()) {
            subject_pairs = subject_pairs
                .combine(normal_cram.first()) { pair, n -> [pair[0], pair[1], n] }
        }
    }

    paired_pairs = subject_pairs.filter { subject, tumor, normal -> normal != null }
    tumor_only_pairs = subject_pairs.filter { subject, tumor, normal -> normal == null }

    paired_pairs.into { pairs_for_mutect; pairs_for_contam }
    tumor_only_pairs.into { to_for_mutect; to_for_contam }

    // Split intervals
    intervals = GATK_SPLITINTERVALS(ref_split, dict_ch, intervals_seed_main, params.scatter_count).out.intervals

    // Panel of normals
    pon_ch = Channel.empty()
    if (!params.panel_of_normals && params.pon_crams) {
        pon_inputs = Channel.fromPath(params.pon_crams, checkIfExists: true).map { cram ->
            def meta = [id: cram.baseName.replaceAll(/\\.cram$|\\.bam$/, ''), sample: cram.baseName.replaceAll(/\\.cram$|\\.bam$/, '')]
            [meta, cram]
        }
        pon_valid = VALIDATE_CRAM(pon_inputs, ref_pon).out.cram_crai
        pon_md = GATK_MARKDUPLICATES(pon_valid, ref_pon)
        pon_bqsr = GATK_BASERECALIBRATOR(pon_md.out.cram_crai, ref_pon, ks_pon)
        pon_recals = GATK_APPLYBQSR(pon_md.out.cram_crai, pon_bqsr.out.table, ref_pon)

        pon_intervals = GATK_SPLITINTERVALS(ref_split, dict_ch, intervals_seed_pon, params.scatter_count).out.intervals
        pon_shards = pon_intervals.cross(pon_recals.out.cram_crai).map { interval, tupleVal ->
            def (meta, cram, crai) = tupleVal
            [interval, meta, cram, crai]
        }
        GATK_MUTECT2_NORMAL(pon_shards, ref_pon)
        pon_grouped = GATK_MUTECT2_NORMAL.out.vcf.groupTuple()
        GATK_MERGE_PON_NORMAL(pon_grouped, ref_pon)
        pon_final = GATK_MERGE_PON_NORMAL.out.vcf.collect()
        GATK_CREATE_SOMATIC_PON(pon_final)
        pon_ch = GATK_CREATE_SOMATIC_PON.out.pon
    } else if (params.panel_of_normals) {
        pon_ch = Channel.fromPath(params.panel_of_normals, checkIfExists: true).first()
    }

    // Scatter Mutect2 across intervals
    pon_ready = pon_ch.ifEmpty { Channel.value(null) }

    paired_scatter = intervals.cross(pairs_for_mutect).map { interval, triple ->
        def (subject, t, n) = triple
        def (t_meta, t_cram, t_crai) = t
        def (n_meta, n_cram, n_crai) = n
        [interval, t_meta, t_cram, t_crai, n_meta, n_cram, n_crai]
    }
    to_scatter = intervals.cross(to_for_mutect).map { interval, triple ->
        def (subject, t, _) = triple
        def (t_meta, t_cram, t_crai) = t
        [interval, t_meta, t_cram, t_crai]
    }

    mutect_paired = GATK_MUTECT2_SCATTER(paired_scatter, ref_mutect, pon_ready)
    mutect_tonly  = GATK_MUTECT2_SCATTER_TONLY(to_scatter, ref_mutect, pon_ready)

    vcf_shards = Channel.empty()
        .mix(mutect_paired.out.vcf)
        .mix(mutect_tonly.out.vcf)
    f1r2_shards = Channel.empty()
        .mix(mutect_paired.out.f1r2)
        .mix(mutect_tonly.out.f1r2)

    vcf_group = vcf_shards.groupTuple()
    merged_vcf = GATK_MERGEVCFS(vcf_group, ref_merge)

    f1r2_group = f1r2_shards.groupTuple()
    lrom = GATK_LEARNREADORIENTATIONMODEL(f1r2_group)

    // Contamination
    contam_inputs_paired = pairs_for_contam.map { subject, t, n ->
        def (t_meta, t_cram, t_crai) = t
        def (n_meta, n_cram, n_crai) = n
        [t_meta, t_cram, t_crai, n_meta, n_cram, n_crai]
    }
    contam_inputs_tonly = to_for_contam.map { subject, t, _ ->
        def (t_meta, t_cram, t_crai) = t
        [t_meta, t_cram, t_crai]
    }

    contam_paired = GATK_CONTAMINATION(contam_inputs_paired, ks_contam_paired)
    contam_tonly  = GATK_CONTAMINATION_TONLY(contam_inputs_tonly, ks_contam_tonly)
    contam = Channel.empty()
        .mix(contam_paired.out.contamination)
        .mix(contam_tonly.out.contamination)
    contam_segments = Channel.empty()
        .mix(contam_paired.out.segments)
        .mix(contam_tonly.out.segments)

    // Align channels by sample for filtering
    vcf_k       = merged_vcf.out.vcf.map { meta, vcf, tbi -> [meta.id, meta, vcf, tbi] }
    contam_k    = contam.map { meta, table -> [meta.id, table] }
    segments_k  = contam_segments.map { meta, seg -> [meta.id, seg] }
    lrom_k      = lrom.out.model.map { meta, model -> [meta.id, model] }

    vcf_contam = vcf_k.join(contam_k) { it[0] } { it[0] }
        .map { key, v, c -> [key, v[1], v[2], v[3], c[1]] }
    vcf_contam_seg = vcf_contam.join(segments_k) { it[0] } { it[0] }
        .map { key, vc, s -> [key, vc[1], vc[2], vc[3], vc[4], s[1]] }
    filter_input = vcf_contam_seg.join(lrom_k) { it[0] } { it[0] }
        .map { key, vcs, m -> [vcs[1], vcs[2], vcs[3], vcs[4], m[1]] } // meta, vcf, tbi, contam, seg, model

    GATK_FILTERMUTECTCALLS(
        filter_input,
        ref_filter
    )

    // QC
    qc_inputs = fastqc_zip.mix(md.out.metrics)
    MULTIQC(qc_inputs)
}
