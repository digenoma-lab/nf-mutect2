#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// ------------------------------------------------------------------
// Parameters (WGS hg38 defaults; SLURM profile in nextflow.config)
// ------------------------------------------------------------------
params.reads            = null          // Glob for paired FASTQs: "data/*_{R1,R2}.fastq.gz"
params.crams            = null          // Glob for existing CRAMs/BAMs
params.reference        = null          // hg38 FASTA
params.known_sites      = null          // gnomAD/common VCF for BQSR & contamination
params.germline_resource = null         // Mutect2 germline resource (e.g. af-only-gnomad.vcf.gz)
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
params.run_bqsr         = false         // optional BQSR; disabled by default
params.debug            = false

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

    stub:
    """
    touch genome.dict
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

    stub:
    """
    touch ${reference}.fai
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
        if(sn ~ /^(chr)?([1-9]|1[0-9]|2[0-2]|X|Y)\$/){print sn,0,ln}
    }' ${dict} > canonical_autosomes_sex.bed
    """

    stub:
    """
    printf "chr1\t0\t1000000\n" > canonical_autosomes_sex.bed
    """
}

process MERGE_INTERVALS {
    tag "merge_intervals"
    publishDir "${params.outdir}/intervals", mode: 'copy', saveAs: { "merged_intervals.bed" }

    input:
    path primary
    path extra

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

    stub:
    """
    cp ${primary} merged_intervals.bed
    """
}

process FASTQC {
    tag "$meta.id"
    publishDir "${params.outdir}/qc/fastqc", mode: "copy"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.zip"), emit: zip

    script:
    """
    fastqc -q ${reads}
    """

    stub:
    """
    touch ${meta.id}_fastqc.zip
    """
}

process BWAMEM {
    tag "${meta.id}-${part}"
    publishDir "${params.outdir}/align", mode: "copy", pattern: "*.mkdup.cram*"

    input:
    tuple val(meta), val(part), file(read1), file(read2)
    path reference

    output:
    tuple val(meta), val(part), path("${meta.id}-${part}.mkdup.cram"), path("${meta.id}-${part}.mkdup.cram.crai"), emit: bams

    script:
    def rg = "@RG\\tID:${meta.id}-${part}\\tSM:${meta.sample}\\tPL:ILLUMINA\\tLB:${meta.sample}"
    """
    bwa-mem2 mem -t ${task.cpus} -R "${rg}" ${reference} ${read1} ${read2} | \\
      samtools view -Sb - | \\
      samtools fixmate -m - - | \\
      samtools sort -@ ${task.cpus} - | \\
      samtools markdup -@ ${task.cpus} --write-index -O CRAM --reference ${reference} - ${meta.id}-${part}.mkdup.cram
    """

    stub:
    """
    touch ${meta.id}-${part}.mkdup.cram ${meta.id}-${part}.mkdup.cram.crai
    """
}

process MERGEB {
    tag "${meta.id}-merge"
    publishDir "${params.outdir}/align", mode: "copy", pattern: "*.cram*"

    input:
    tuple val(meta), val(parts), path(cramFiles), path(craiFiles)
    path reference

    output:
    tuple val(meta), path("${meta.id}.cram"), path("${meta.id}.cram.crai"), emit: cram_crai

    script:
    def filesb = cramFiles instanceof List ? cramFiles : [cramFiles]
    if (filesb.size() == 1) {
        """
        ln -s ${cramFiles[0]} ${meta.id}.cram
        ln -s ${craiFiles[0]} ${meta.id}.cram.crai
        """
    } else {
        """
        samtools merge --write-index --reference ${reference} -O CRAM -@ ${task.cpus} -f ${meta.id}.cram ${cramFiles}
        """
    }

    stub:
    """
    touch ${meta.id}.cram ${meta.id}.cram.crai
    """
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

    stub:
    """
    touch ${cram}.crai
    """
}

process VALIDATE_CRAM_PON {
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

    stub:
    """
    touch ${cram}.crai
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

    stub:
    """
    touch ${meta.id}_recal.table
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

    stub:
    """
    touch ${meta.id}_recal.cram ${meta.id}_recal.cram.crai
    """
}

process GATK_SPLITINTERVALS {
    tag "$reference.baseName"
    publishDir "${params.outdir}/intervals", mode: 'copy'

    input:
    path reference
    path dict
    path base_intervals
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

    stub:
    """
    mkdir -p intervals
    touch intervals/0000-scattered.interval_list intervals/0001-scattered.interval_list
    """
}

process GATK_MUTECT2_SCATTER {
    tag "$tumor_meta.id:${interval.baseName}"
    publishDir "${params.outdir}/variants/shards", mode: 'copy'

    input:
    tuple path(interval), val(tumor_meta), path(tumor_cram), path(tumor_crai), val(normal_meta), path(normal_cram), path(normal_crai)
    path reference
    path germline_resource
    val pon

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
      --germline-resource ${germline_resource} \\
      ${ponArg} \\
      --f1r2-tar-gz ${tumor_meta.id}_${interval.baseName}_f1r2.tar.gz \\
      -O ${tumor_meta.id}_${interval.baseName}.vcf.gz
    """

    stub:
    """
    touch ${tumor_meta.id}_${interval.baseName}.vcf.gz
    touch ${tumor_meta.id}_${interval.baseName}.vcf.gz.tbi
    touch ${tumor_meta.id}_${interval.baseName}_f1r2.tar.gz
    touch ${tumor_meta.id}_${interval.baseName}.stats
    """
}

process GATK_MUTECT2_SCATTER_TONLY {
    tag "$tumor_meta.id:${interval.baseName}"
    publishDir "${params.outdir}/variants/shards", mode: 'copy'

    input:
    tuple path(interval), val(tumor_meta), path(tumor_cram), path(tumor_crai)
    path reference
    path germline_resource
    val pon

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
      --germline-resource ${germline_resource} \\
      ${ponArg} \\
      --f1r2-tar-gz ${tumor_meta.id}_${interval.baseName}_f1r2.tar.gz \\
      -O ${tumor_meta.id}_${interval.baseName}.vcf.gz
    """

    stub:
    """
    touch ${tumor_meta.id}_${interval.baseName}.vcf.gz
    touch ${tumor_meta.id}_${interval.baseName}.vcf.gz.tbi
    touch ${tumor_meta.id}_${interval.baseName}_f1r2.tar.gz
    touch ${tumor_meta.id}_${interval.baseName}.stats
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

    stub:
    """
    touch ${meta.id}_merged.vcf.gz ${meta.id}_merged.vcf.gz.tbi
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

    stub:
    """
    touch ${meta.id}_read_orientation_model.tar.gz
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
    path normal_pileups

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
    tuple val(meta), path(vcf), path(tbi), path(contamination), path(segments), path(orientation_model)
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
      --tumor-segmentation ${segments} \\
      --ob-priors ${orientation_model} \\
      -O ${meta.id}_filtered.vcf.gz \\
      --filtering-stats ${meta.id}_filtering_stats.tsv
    """

    stub:
    """
    touch ${meta.id}_filtered.vcf.gz ${meta.id}_filtered.vcf.gz.tbi ${meta.id}_filtering_stats.tsv
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

    stub:
    """
    mkdir -p multiqc_data
    touch multiqc_report.html
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

    stub:
    """
    touch ${tumor_meta.id}_contamination.table ${tumor_meta.id}_segments.table
    """
}

process GATK_CONTAMINATION_TONLY {
    tag "$tumor_meta.id"

    input:
    tuple val(tumor_meta), path(tumor_cram), path(tumor_crai)
    path known_sites

    output:
    tuple val(tumor_meta), path("${tumor_meta.id}_contamination.table"), emit: contamination
    tuple val(tumor_meta), path("${tumor_meta.id}_segments.table"), emit: segments

    script:
    """
    gatk GetPileupSummaries -I ${tumor_cram} -V ${known_sites} -L ${known_sites} -O tumor_pileups.table
    gatk CalculateContamination \\
      -I tumor_pileups.table \\
      -O ${tumor_meta.id}_contamination.table \\
      --tumor-segmentation ${tumor_meta.id}_segments.table
    """

    stub:
    """
    touch ${tumor_meta.id}_contamination.table ${tumor_meta.id}_segments.table
    """
}

// PON creation
process GATK_MUTECT2_NORMAL {
    tag "$meta.id:${interval.baseName}"
    publishDir "${params.outdir}/pon/shards", mode: 'copy'

    input:
    tuple path(interval), val(meta), path(cram), path(crai)
    path reference
    path germline_resource

    output:
    tuple val(meta), path("${meta.id}_${interval.baseName}_pon.vcf.gz"), emit: vcf

    script:
    """
    gatk Mutect2 \\
      -R ${reference} \\
      -I ${cram} \\
      -tumor ${meta.sample} \\
      -L ${interval} \\
      --germline-resource ${germline_resource} \\
      --max-mnp-distance 0 \\
      -O ${meta.id}_${interval.baseName}_pon.vcf.gz
    """

    stub:
    """
    touch ${meta.id}_${interval.baseName}_pon.vcf.gz
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

    stub:
    """
    touch ${meta.id}_pon_merged.vcf.gz
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

    stub:
    """
    touch panel_of_normals.vcf.gz
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
    if (!params.germline_resource) error "Provide --germline_resource (e.g. af-only-gnomad.vcf.gz) for Mutect2"
    if (!params.sample_sheet && (!params.tumor_sample))
        error "Provide --tumor_sample (and --normal_sample if paired) when not using --sample_sheet"

    reference_ch    = Channel.value(file(params.reference))
    known_sites_ch  = Channel.value(file(params.known_sites))
    germline_ch     = Channel.value(file(params.germline_resource))

    GATK_CREATE_DICT(reference_ch)
    dict_ch = GATK_CREATE_DICT.out.dict
    MAKE_CANONICAL_BED(dict_ch)
    canonical_bed = MAKE_CANONICAL_BED.out.bed
    SAMTOOLS_FAIDX(reference_ch)
    base_interval_seed = params.intervals ? Channel.value(file(params.intervals))
                                          : (params.canonical_only ? canonical_bed : Channel.empty())
    if (params.extra_intervals) {
        extra_intervals_ch = Channel.value(file(params.extra_intervals))
        MERGE_INTERVALS(base_interval_seed, extra_intervals_ch)
        base_interval_seed = MERGE_INTERVALS.out.bed
    }
    intervals_seed_main = base_interval_seed

    // Input handling (sample sheet preferred)
    crams_valid = null
    fastqc_zip = Channel.empty()

    if (params.sample_sheet) {
        sample_sheet_file = file(params.sample_sheet)
        sample_sheet_dir = sample_sheet_file.parent
        resolve_sheet_path = { p ->
            if (!p) return null
            return p.startsWith('/') ? file(p) : file(sample_sheet_dir.resolve(p).toString())
        }

        sample_rows = Channel.fromPath(params.sample_sheet, checkIfExists: true)
            .splitCsv(header:true)
            .map { row ->
                def meta = [
                    id     : row.sample_id ?: row.subject ?: row.type,
                    sample : row.sample_id ?: row.subject ?: row.type,
                    subject: row.subject ?: row.sample_id ?: row.type,
                    type   : row.type?.toLowerCase()
                ]
                def fastqs = (row.fastq1 && row.fastq2) ? [resolve_sheet_path.call(row.fastq1), resolve_sheet_path.call(row.fastq2)] : null
                def cram = resolve_sheet_path.call(row.cram)
                [meta, fastqs, cram]
            }
            .filter { meta, fastqs, cram -> meta.type in ['tumor','normal'] }

        fastq_samples = sample_rows.filter { meta, fastqs, cram -> fastqs }.map { meta, fastqs, cram -> [meta, fastqs] }
        cram_samples  = sample_rows.filter { meta, fastqs, cram -> cram }.map { meta, fastqs, cram -> [meta, cram] }

        FASTQC(fastq_samples)
        fastqc_zip = FASTQC.out.zip
        bwamem_inputs = fastq_samples.map { meta, fastqs ->
            def part = fastqs[0].baseName.replaceAll(/[^A-Za-z0-9._-]/, '_')
            [meta, part, fastqs[0], fastqs[1]]
        }
        BWAMEM(bwamem_inputs, reference_ch)
        merged_parts = BWAMEM.out.bams.groupTuple()
        MERGEB(merged_parts, reference_ch)
        crams_from_fastq = MERGEB.out.cram_crai
        VALIDATE_CRAM(cram_samples, reference_ch)
        crams_from_cram = VALIDATE_CRAM.out.cram_crai
        crams_valid = crams_from_fastq.mix(crams_from_cram)

    } else if (params.crams) {
        crams_in = Channel.fromPath(params.crams, checkIfExists: true).map { cram ->
            def meta = [id: cram.baseName.replaceAll(/\\.cram$|\\.bam$/, ''), sample: cram.baseName.replaceAll(/\\.cram$|\\.bam$/, ''), subject: cram.baseName.replaceAll(/\\.cram$|\\.bam$/, ''), type: 'tumor']
            [meta, cram]
        }
        VALIDATE_CRAM(crams_in, reference_ch)
        crams_valid = VALIDATE_CRAM.out.cram_crai
    } else {
        reads_ch = Channel.fromFilePairs(params.reads, checkIfExists: true).map { sample, reads ->
            def meta = [id: sample, sample: sample, subject: sample, type: sample == params.normal_sample ? 'normal' : 'tumor']
            [meta, reads]
        }
        FASTQC(reads_ch)
        fastqc_zip = FASTQC.out.zip
        bwamem_inputs = reads_ch.map { meta, reads ->
            [meta, 'part1', reads[0], reads[1]]
        }
        BWAMEM(bwamem_inputs, reference_ch)
        merged_parts = BWAMEM.out.bams.groupTuple()
        MERGEB(merged_parts, reference_ch)
        crams_valid = MERGEB.out.cram_crai
    }

    // Preprocess: duplicate marking is already done during alignment.
    if (params.run_bqsr) {
        GATK_BASERECALIBRATOR(crams_valid, reference_ch, known_sites_ch)
        bqsr_table = GATK_BASERECALIBRATOR.out.table
        GATK_APPLYBQSR(crams_valid, bqsr_table, reference_ch)
        final_crams = GATK_APPLYBQSR.out.cram_crai
    } else {
        final_crams = crams_valid
    }

    // Pairing
    subject_pairs = null
    if (params.sample_sheet) {
        grouped = final_crams
            .map { meta, cram, crai -> [meta.subject, meta, cram, crai] }
            .groupTuple()
        subject_pairs = grouped.map { subject, metas, crams, crais ->
            def tuples = (0..<metas.size()).collect { i -> [metas[i], crams[i], crais[i]] }
            def tumor = tuples.find { it[0].type == 'tumor' }
            def normal = tuples.find { it[0].type == 'normal' }
            tumor ? [subject, tumor, normal] : null
        }.filter { it != null }
    } else {
        tumor_cram = final_crams.filter { meta, cram, crai -> meta.sample == params.tumor_sample }
        normal_cram = final_crams.filter { meta, cram, crai -> params.normal_sample ? meta.sample == params.normal_sample : false }
        if (params.normal_sample) {
            subject_pairs = tumor_cram
                .combine(normal_cram.first()) { t, n -> [t[0].sample, t, n] }
        } else {
            subject_pairs = tumor_cram.map { t -> [t[0].sample, t, null] }
        }
    }

    pairs_for_mutect = subject_pairs.filter { rec -> rec[2] != null }
    to_for_mutect = subject_pairs.filter { rec -> rec[2] == null }


    // Split intervals
    GATK_SPLITINTERVALS(reference_ch, dict_ch, intervals_seed_main, params.scatter_count)
    intervals = GATK_SPLITINTERVALS.out.intervals
    intervals_all = intervals.collect()
    intervals_for_paired = intervals_all.flatMap { it }
    intervals_for_tonly = intervals_all.flatMap { it }


    // Panel of normals
    pon_ch = Channel.empty()
    if (!params.panel_of_normals && params.pon_crams) {
        pon_inputs = Channel.fromPath(params.pon_crams, checkIfExists: true).map { cram ->
            def meta = [id: cram.baseName.replaceAll(/\\.cram$|\\.bam$/, ''), sample: cram.baseName.replaceAll(/\\.cram$|\\.bam$/, '')]
            [meta, cram]
        }
        VALIDATE_CRAM_PON(pon_inputs, reference_ch)
        pon_valid = VALIDATE_CRAM_PON.out.cram_crai
        // Keep PON path independent from primary BQSR/splitting invocations to avoid process re-use conflicts.
        pon_prepped = pon_valid
        pon_intervals = intervals_all.flatMap { it }
        pon_shards = pon_intervals.combine(pon_prepped)
        GATK_MUTECT2_NORMAL(pon_shards, reference_ch, germline_ch)
        pon_grouped = GATK_MUTECT2_NORMAL.out.vcf.groupTuple()
        GATK_MERGE_PON_NORMAL(pon_grouped, reference_ch)
        pon_final = GATK_MERGE_PON_NORMAL.out.vcf.map { meta, vcf -> vcf }.collect()
        GATK_CREATE_SOMATIC_PON(pon_final)
        pon_ch = GATK_CREATE_SOMATIC_PON.out.pon
    } else if (params.panel_of_normals) {
        pon_ch = Channel.value(file(params.panel_of_normals))
    }

    // Scatter Mutect2 across intervals
    pon_ready = pon_ch.ifEmpty(null)

    paired_scatter = intervals_for_paired.combine(pairs_for_mutect).map { interval, subject, t, n ->
        def (t_meta, t_cram, t_crai) = t
        def (n_meta, n_cram, n_crai) = n
        [interval, t_meta, file(t_cram.toString()), file(t_crai.toString()), n_meta, file(n_cram.toString()), file(n_crai.toString())]
    }
    to_scatter = intervals_for_tonly.combine(to_for_mutect).map { interval, subject, t, _ ->
        def (t_meta, t_cram, t_crai) = t
        [interval, t_meta, file(t_cram.toString()), file(t_crai.toString())]
    }


    GATK_MUTECT2_SCATTER(paired_scatter, reference_ch, germline_ch, pon_ready)
    GATK_MUTECT2_SCATTER_TONLY(to_scatter, reference_ch, germline_ch, pon_ready)

    vcf_shards = Channel.empty()
        .mix(GATK_MUTECT2_SCATTER.out.vcf)
        .mix(GATK_MUTECT2_SCATTER_TONLY.out.vcf)
    f1r2_shards = Channel.empty()
        .mix(GATK_MUTECT2_SCATTER.out.f1r2)
        .mix(GATK_MUTECT2_SCATTER_TONLY.out.f1r2)

    vcf_group = vcf_shards
        .map { meta, vcf, tbi -> [meta, vcf] }
        .groupTuple()
    GATK_MERGEVCFS(vcf_group, reference_ch)
    merged_vcf = GATK_MERGEVCFS.out.vcf

    f1r2_group = f1r2_shards.groupTuple()
    GATK_LEARNREADORIENTATIONMODEL(f1r2_group)
    lrom = GATK_LEARNREADORIENTATIONMODEL.out.model

    // Contamination
    contam_inputs_paired = paired_scatter
        .map { interval, t_meta, t_cram, t_crai, n_meta, n_cram, n_crai ->
            [t_meta, t_cram, t_crai, n_meta, n_cram, n_crai]
        }
        .unique { row -> row[0].id }
    contam_inputs_tonly = to_scatter
        .map { interval, t_meta, t_cram, t_crai -> [t_meta, t_cram, t_crai] }
        .unique { row -> row[0].id }

    GATK_CONTAMINATION(contam_inputs_paired, known_sites_ch)
    GATK_CONTAMINATION_TONLY(contam_inputs_tonly, known_sites_ch)
    contam = Channel.empty()
        .mix(GATK_CONTAMINATION.out.contamination)
        .mix(GATK_CONTAMINATION_TONLY.out.contamination)
    contam_segments = Channel.empty()
        .mix(GATK_CONTAMINATION.out.segments)
        .mix(GATK_CONTAMINATION_TONLY.out.segments)

    // Align channels by sample for filtering
    vcf_k       = merged_vcf.map { meta, vcf, tbi -> [meta.id, meta, vcf, tbi] }
    contam_k    = contam.map { meta, table -> [meta.id, table] }
    segments_k  = contam_segments.map { meta, seg -> [meta.id, seg] }
    lrom_k      = lrom.map { meta, model -> [meta.id, model] }

    vcf_contam = vcf_k.join(contam_k)
        .map { key, meta, vcf, tbi, table -> [key, meta, vcf, tbi, table] }
    vcf_contam_seg = vcf_contam.join(segments_k)
        .map { key, meta, vcf, tbi, table, seg -> [key, meta, vcf, tbi, table, seg] }
    filter_input = vcf_contam_seg.join(lrom_k)
        .map { key, meta, vcf, tbi, table, seg, model -> [meta, vcf, tbi, table, seg, model] }

    GATK_FILTERMUTECTCALLS(
        filter_input,
        reference_ch
    )

    // QC
    qc_inputs = fastqc_zip.map { meta, zip -> zip }
    MULTIQC(qc_inputs)
}
