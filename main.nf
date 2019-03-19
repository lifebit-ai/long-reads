#!/usr/bin/env nextflow

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

    input:
    set val(name), file(reads), file(fasta) from minimap2

    output:
    file("${name}.sam") into results

    script:
    """
    minimap2 -ax map-ont -t ${task.cpus} $fasta $reads > ${name}.sam
    """
}


