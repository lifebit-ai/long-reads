#!/usr/bin/env nextflow

int threads = Runtime.getRuntime().availableProcessors()

Channel
    .fromPath(params.fasta)
    .ifEmpty { exit 1, "fasta annotation file not found: ${params.fasta}" }
    .into { fasta; fasta_minimap2 }

Channel
      .fromPath(params.reads)
      .map { file -> tuple(file.baseName, file) }
      .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes" }
      .into { reads; reads_minimap2 }

reads.subscribe { println "value: $it"}

minimap2 = reads_minimap2.combine(fasta_minimap2)

process minimap2 {
    tag "$reads"
    publishDir "${params.outdir}", mode: 'copy'
    container 'evolbioinfo/minimap2:v2.14'
    cpus threads

    input:
    set val(name), file(reads), file(fasta) from minimap2

    output:
    set val(name), file("${name}.sam") into mapped_reads

    script:
    """
    minimap2 -ax map-ont -t ${task.cpus} $fasta $reads > ${name}.sam
    """
}

process bwa_sort {
    tag "$sam"
    container 'lifebitai/samtools:latest'

    input:
    set val(name), file(sam) from mapped_reads

    output:
    set val(name), file("${name}-sorted.bam") into sorted_bam

    """
    samtools sort -o ${name}-sorted.bam -O BAM $sam
    """
}

process mark_duplicates {
    tag "$bam"
    container 'broadinstitute/gatk:latest'

    input:
    set val(name), file(bam) from sorted_bam

    output:
    set val(name), file("${name}.bam"), file("${name}.bai") into marked_dup
    file ("${name}.bam.metrics") into mark_dup_report

    """
    gatk MarkDuplicates \
    -I ${bam} \
    --CREATE_INDEX true \
    -M ${name}.bam.metrics \
    -O ${name}.bam
    """
}


